/*
 * fasttrop C core. Refer to trop_core.h for the interface and the conventions.
 *
 * The weighted low-rank fit has two solvers. The default solver (fit_ws_prox)
 * is an accelerated proximal gradient method on the low-rank term. At each
 * step, it calculates the unit and time effects again, and it uses an adaptive
 * restart.
 *
 * The block solver (fit_ws) is the block coordinate method of the reference
 * implementations. It calculates the exact weighted means for the unit and
 * time effects. Then it does FISTA steps with soft thresholds on the singular
 * values. The objective is convex. Thus, both solvers find the same minimum in
 * the limits of the tolerance.
 */
#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "trop_core.h"

#define SVT_EPS 1e-10

static int nthreads(const trop_opts *o) { return (o->threads < 1) ? 1 : o->threads; }

/* ------------------------------------------------------------------------ */
/* Workspace                                                                */
/* ------------------------------------------------------------------------ */

typedef struct {
  int N, T, K, lwork;
  double *alpha, *beta, *L, *aold, *bold, *Lold, *Lin, *Lprev, *G;
  double *A, *U, *S, *VT, *Us, *work, *w, *rs, *cs;
  int *iwork, *queue;
  unsigned char *seen;
} ws_t;

static void ws_free(ws_t *ws) {
  if (!ws) return;
  free(ws->alpha); free(ws->beta); free(ws->L); free(ws->aold); free(ws->bold);
  free(ws->Lold); free(ws->Lin); free(ws->Lprev); free(ws->G); free(ws->A);
  free(ws->U); free(ws->S); free(ws->VT); free(ws->Us); free(ws->work);
  free(ws->w); free(ws->rs); free(ws->cs); free(ws->iwork); free(ws->queue);
  free(ws->seen);
  free(ws);
}

static ws_t *ws_new(int N, int T) {
  ws_t *ws = (ws_t *) calloc(1, sizeof(ws_t));
  if (!ws) return NULL;
  size_t NT = (size_t) N * (size_t) T;
  int K = (N < T) ? N : T;
  ws->N = N; ws->T = T; ws->K = K;
  ws->alpha = (double *) calloc((size_t) N, sizeof(double));
  ws->beta = (double *) calloc((size_t) T, sizeof(double));
  ws->aold = (double *) malloc((size_t) N * sizeof(double));
  ws->bold = (double *) malloc((size_t) T * sizeof(double));
  ws->rs = (double *) malloc((size_t) N * sizeof(double));
  ws->cs = (double *) malloc((size_t) T * sizeof(double));
  ws->L = (double *) calloc(NT, sizeof(double));
  ws->Lold = (double *) malloc(NT * sizeof(double));
  ws->Lin = (double *) malloc(NT * sizeof(double));
  ws->Lprev = (double *) malloc(NT * sizeof(double));
  ws->G = (double *) malloc(NT * sizeof(double));
  ws->A = (double *) malloc(NT * sizeof(double));
  ws->w = (double *) malloc(NT * sizeof(double));
  ws->U = (double *) malloc((size_t) N * (size_t) K * sizeof(double));
  ws->Us = (double *) malloc((size_t) N * (size_t) K * sizeof(double));
  ws->VT = (double *) malloc((size_t) K * (size_t) T * sizeof(double));
  ws->S = (double *) malloc((size_t) K * sizeof(double));
  ws->iwork = (int *) malloc((size_t) 8 * (size_t) K * sizeof(int));
  ws->queue = (int *) malloc((size_t) (N + T) * sizeof(int));
  ws->seen = (unsigned char *) malloc((size_t) (N + T));
  if (!ws->alpha || !ws->beta || !ws->aold || !ws->bold || !ws->rs || !ws->cs ||
      !ws->L || !ws->Lold || !ws->Lin || !ws->Lprev || !ws->G || !ws->A || !ws->w ||
      !ws->U || !ws->Us || !ws->VT || !ws->S || !ws->iwork || !ws->queue || !ws->seen) {
    ws_free(ws);
    return NULL;
  }
  /* LAPACK workspace query for the thin SVD. */
  int m = N, n = T, lda = N, ldu = N, ldvt = K, info = 0, lw = -1;
  double wq = 0.0;
  F77_CALL(dgesdd)("S", &m, &n, ws->A, &lda, ws->S, ws->U, &ldu, ws->VT, &ldvt,
                   &wq, &lw, ws->iwork, &info FCONE);
  ws->lwork = (info == 0 && wq > 0) ? (int) wq + 1 : 7 * K * K + 4 * (N + T) + 64;
  ws->work = (double *) malloc((size_t) ws->lwork * sizeof(double));
  if (!ws->work) { ws_free(ws); return NULL; }
  return ws;
}

static void ws_cold(ws_t *ws) {
  memset(ws->alpha, 0, (size_t) ws->N * sizeof(double));
  memset(ws->beta, 0, (size_t) ws->T * sizeof(double));
  memset(ws->L, 0, (size_t) ws->N * (size_t) ws->T * sizeof(double));
}

/* ------------------------------------------------------------------------ */
/* Weighted low-rank fit                                                     */
/* ------------------------------------------------------------------------ */

/* Singular-value soft-thresholding: Lout = U diag((s - thr)_+) V'. */
static void svt(ws_t *ws, const double *Gin, double thr, double *Lout) {
  int N = ws->N, T = ws->T, K = ws->K;
  size_t NT = (size_t) N * (size_t) T;
  if (!(thr > 0.0)) { memcpy(Lout, Gin, NT * sizeof(double)); return; }
  if (!isfinite(thr)) { memset(Lout, 0, NT * sizeof(double)); return; }
  memcpy(ws->A, Gin, NT * sizeof(double));
  int m = N, n = T, lda = N, ldu = N, ldvt = K, info = 0, lw = ws->lwork;
  F77_CALL(dgesdd)("S", &m, &n, ws->A, &lda, ws->S, ws->U, &ldu, ws->VT, &ldvt,
                   ws->work, &lw, ws->iwork, &info FCONE);
  if (info != 0) { memset(Lout, 0, NT * sizeof(double)); return; }
  int r = 0;
  while (r < K && ws->S[r] - thr > SVT_EPS) r++;
  if (r == 0) { memset(Lout, 0, NT * sizeof(double)); return; }
  for (int k = 0; k < r; k++) {
    double sk = ws->S[k] - thr;
    const double *uk = ws->U + (size_t) N * k;
    double *usk = ws->Us + (size_t) N * k;
    for (int i = 0; i < N; i++) usk[i] = uk[i] * sk;
  }
  double one = 1.0, zero = 0.0;
  int ldvt2 = K;
  F77_CALL(dgemm)("N", "N", &m, &n, &r, &one, ws->Us, &m, ws->VT, &ldvt2, &zero,
                  Lout, &m FCONE FCONE);
}

static double maxabsdiff(const double *a, const double *b, size_t n) {
  double d = 0.0;
  for (size_t i = 0; i < n; i++) {
    double v = fabs(a[i] - b[i]);
    if (v > d) d = v;
  }
  return d;
}

/* Core solver on the workspace. Uses ws->alpha, ws->beta, ws->L as the start
 * point and leaves the solution there. Returns 1 if converged. */
static int fit_ws(ws_t *ws, const double *Y, const double *w, double lam,
                  const trop_opts *o, int *iters) {
  int N = ws->N, T = ws->T;
  size_t NT = (size_t) N * (size_t) T;
  double *a = ws->alpha, *b = ws->beta, *L = ws->L;
  int lowrank = isfinite(lam) ? 1 : 0;
  double wmax = 0.0;
  for (int i = 0; i < N; i++) ws->rs[i] = 0.0;
  for (int t = 0; t < T; t++) {
    double c = 0.0;
    const double *wt = w + (size_t) N * t;
    for (int i = 0; i < N; i++) {
      double v = wt[i];
      c += v;
      ws->rs[i] += v;
      if (v > wmax) wmax = v;
    }
    ws->cs[t] = c;
  }
  if (!lowrank) memset(L, 0, NT * sizeof(double));
  double thr = (wmax > 0.0) ? lam / (2.0 * wmax) : lam / 2.0;
  double winv = (wmax > 0.0) ? 1.0 / wmax : 1.0;
  int converged = 0, it;
  for (it = 0; it < o->max_iter; it++) {
    memcpy(ws->aold, a, (size_t) N * sizeof(double));
    memcpy(ws->bold, b, (size_t) T * sizeof(double));
    if (lowrank) memcpy(ws->Lold, L, NT * sizeof(double));
    /* unit effects */
    for (int i = 0; i < N; i++) a[i] = 0.0;
    for (int t = 0; t < T; t++) {
      const double *wt = w + (size_t) N * t, *yt = Y + (size_t) N * t, *lt = L + (size_t) N * t;
      double bt = b[t];
      for (int i = 0; i < N; i++) if (wt[i] > 0.0) a[i] += wt[i] * (yt[i] - lt[i] - bt);
    }
    for (int i = 0; i < N; i++) a[i] = (ws->rs[i] > 0.0) ? a[i] / ws->rs[i] : 0.0;
    /* time effects */
    for (int t = 0; t < T; t++) {
      const double *wt = w + (size_t) N * t, *yt = Y + (size_t) N * t, *lt = L + (size_t) N * t;
      double num = 0.0;
      for (int i = 0; i < N; i++) if (wt[i] > 0.0) num += wt[i] * (yt[i] - lt[i] - a[i]);
      b[t] = (ws->cs[t] > 0.0) ? num / ws->cs[t] : 0.0;
    }
    /* low-rank component */
    if (lowrank) {
      if (!(lam > 0.0)) {
        for (int t = 0; t < T; t++) {
          const double *wt = w + (size_t) N * t, *yt = Y + (size_t) N * t;
          double *lt = L + (size_t) N * t;
          for (int i = 0; i < N; i++) if (wt[i] > 0.0) lt[i] = yt[i] - a[i] - b[t];
        }
      } else {
        double tf = 1.0;
        memcpy(ws->Lprev, L, NT * sizeof(double));
        for (int k = 0; k < o->max_inner; k++) {
          memcpy(ws->Lin, L, NT * sizeof(double));
          double tfn = 0.5 * (1.0 + sqrt(1.0 + 4.0 * tf * tf));
          double mom = (tf - 1.0) / tfn;
          for (int t = 0; t < T; t++) {
            const double *wt = w + (size_t) N * t, *yt = Y + (size_t) N * t;
            const double *lt = L + (size_t) N * t, *lp = ws->Lprev + (size_t) N * t;
            double *gt = ws->G + (size_t) N * t;
            double bt = b[t];
            for (int i = 0; i < N; i++) {
              double lm = lt[i] + mom * (lt[i] - lp[i]);
              gt[i] = (wt[i] > 0.0) ? lm + wt[i] * winv * ((yt[i] - a[i] - bt) - lm) : lm;
            }
          }
          memcpy(ws->Lprev, L, NT * sizeof(double));
          svt(ws, ws->G, thr, L);
          tf = tfn;
          if (maxabsdiff(L, ws->Lin, NT) < o->tol) break;
        }
      }
    }
    double d = maxabsdiff(a, ws->aold, (size_t) N);
    double db = maxabsdiff(b, ws->bold, (size_t) T);
    if (db > d) d = db;
    if (lowrank) {
      double dl = maxabsdiff(L, ws->Lold, NT);
      if (dl > d) d = dl;
    }
    if (d < o->tol) { converged = 1; it++; break; }
  }
  if (iters) *iters = it;
  return converged;
}

/* Weighted two-way fixed effects for R = Y - Z. The function calculates the
 * weighted means for the units and for the periods in turn. It starts from
 * ws->alpha and ws->beta, and it returns the number of sweeps. */
static int twfe_sweeps(ws_t *ws, const double *Y, const double *Z, const double *w,
                       double tol, int max_sweeps) {
  int N = ws->N, T = ws->T, s;
  double *a = ws->alpha, *b = ws->beta;
  for (s = 0; s < max_sweeps; s++) {
    double d = 0.0;
    for (int i = 0; i < N; i++) ws->aold[i] = a[i];
    for (int i = 0; i < N; i++) a[i] = 0.0;
    for (int t = 0; t < T; t++) {
      const double *wt = w + (size_t) N * t, *yt = Y + (size_t) N * t;
      const double *zt = Z ? Z + (size_t) N * t : NULL;
      double bt = b[t];
      for (int i = 0; i < N; i++) {
        if (wt[i] > 0.0) a[i] += wt[i] * (yt[i] - (zt ? zt[i] : 0.0) - bt);
      }
    }
    for (int i = 0; i < N; i++) {
      a[i] = (ws->rs[i] > 0.0) ? a[i] / ws->rs[i] : 0.0;
      double v = fabs(a[i] - ws->aold[i]);
      if (v > d) d = v;
    }
    for (int t = 0; t < T; t++) {
      const double *wt = w + (size_t) N * t, *yt = Y + (size_t) N * t;
      const double *zt = Z ? Z + (size_t) N * t : NULL;
      double num = 0.0;
      for (int i = 0; i < N; i++) {
        if (wt[i] > 0.0) num += wt[i] * (yt[i] - (zt ? zt[i] : 0.0) - a[i]);
      }
      double nb = (ws->cs[t] > 0.0) ? num / ws->cs[t] : 0.0;
      double v = fabs(nb - b[t]);
      if (v > d) d = v;
      b[t] = nb;
    }
    if (d < tol) { s++; break; }
  }
  return s;
}

/* Accelerated proximal gradient on the profile objective. The method does
 * FISTA steps on L for
 * h(L) = min_{a,b} sum w (Y - a - b - L)^2 + lam ||L||_*.
 * At each extrapolated point, it fits the fixed effects again to the
 * tolerance. It uses an adaptive restart from the gradient. The objective and
 * the minimum are the same as for fit_ws. */
static int fit_ws_prox(ws_t *ws, const double *Y, const double *w, double lam,
                       const trop_opts *o, int *iters) {
  int N = ws->N, T = ws->T;
  size_t NT = (size_t) N * (size_t) T;
  double *L = ws->L;
  double wmax = 0.0;
  for (int i = 0; i < N; i++) ws->rs[i] = 0.0;
  for (int t = 0; t < T; t++) {
    double c = 0.0;
    const double *wt = w + (size_t) N * t;
    for (int i = 0; i < N; i++) {
      double v = wt[i];
      c += v;
      ws->rs[i] += v;
      if (v > wmax) wmax = v;
    }
    ws->cs[t] = c;
  }
  double sweep_tol = 0.1 * o->tol;
  int max_steps = o->max_iter * o->max_inner;
  int max_sweeps = 1000;
  if (!isfinite(lam)) {
    memset(L, 0, NT * sizeof(double));
    int sw = twfe_sweeps(ws, Y, NULL, w, sweep_tol, max_steps > max_sweeps ? max_steps : max_sweeps);
    if (iters) *iters = sw;
    return sw < (max_steps > max_sweeps ? max_steps : max_sweeps);
  }
  if (!(lam > 0.0)) {
    /* Unpenalized degenerate case, kept for compatibility with the block scheme. */
    return fit_ws(ws, Y, w, lam, o, iters);
  }
  double winv = (wmax > 0.0) ? 1.0 / wmax : 1.0;
  double thr = (wmax > 0.0) ? lam / (2.0 * wmax) : lam / 2.0;
  double tf = 1.0;
  memcpy(ws->Lprev, L, NT * sizeof(double));
  int converged = 0, k;
  for (k = 0; k < max_steps; k++) {
    double tfn = 0.5 * (1.0 + sqrt(1.0 + 4.0 * tf * tf));
    double mom = (tf - 1.0) / tfn;
    /* Extrapolated point Z (stored in Lin). */
    for (size_t c = 0; c < NT; c++) ws->Lin[c] = L[c] + mom * (L[c] - ws->Lprev[c]);
    twfe_sweeps(ws, Y, ws->Lin, w, sweep_tol, max_sweeps);
    for (int t = 0; t < T; t++) {
      const double *wt = w + (size_t) N * t, *yt = Y + (size_t) N * t;
      const double *zt = ws->Lin + (size_t) N * t;
      double *gt = ws->G + (size_t) N * t;
      double bt = ws->beta[t];
      for (int i = 0; i < N; i++) {
        gt[i] = (wt[i] > 0.0) ? zt[i] + wt[i] * winv * (yt[i] - ws->alpha[i] - bt - zt[i]) : zt[i];
      }
    }
    memcpy(ws->Lprev, L, NT * sizeof(double));
    svt(ws, ws->G, thr, L);
    /* Gradient-based adaptive restart (O'Donoghue and Candes). */
    double ip = 0.0;
    for (size_t c = 0; c < NT; c++) ip += (ws->Lin[c] - L[c]) * (L[c] - ws->Lprev[c]);
    tf = (ip > 0.0) ? 1.0 : tfn;
    if (maxabsdiff(L, ws->Lprev, NT) < o->tol) { converged = 1; k++; break; }
  }
  /* Fixed effects consistent with the final low-rank component. */
  twfe_sweeps(ws, Y, L, w, sweep_tol, max_sweeps);
  if (iters) *iters = k;
  return converged;
}

static int fit_any(ws_t *ws, const double *Y, const double *w, double lam,
                   const trop_opts *o, int *iters) {
  return (o->solver == 1) ? fit_ws(ws, Y, w, lam, o, iters) : fit_ws_prox(ws, Y, w, lam, o, iters);
}

int trop_fit_weighted(int N, int T, const double *Y, const double *w,
                      double lambda_nn, const trop_opts *opt,
                      double *alpha, double *beta, double *L,
                      int *iters, int *converged) {
  ws_t *ws = ws_new(N, T);
  if (!ws) return TROP_ERR_ALLOC;
  int it = 0;
  int conv = fit_any(ws, Y, w, lambda_nn, opt, &it);
  memcpy(alpha, ws->alpha, (size_t) N * sizeof(double));
  memcpy(beta, ws->beta, (size_t) T * sizeof(double));
  memcpy(L, ws->L, (size_t) N * (size_t) T * sizeof(double));
  if (iters) *iters = it;
  if (converged) *converged = conv;
  ws_free(ws);
  return TROP_OK;
}

/* ------------------------------------------------------------------------ */
/* Panel helpers                                                             */
/* ------------------------------------------------------------------------ */

typedef struct {
  int N, T;
  double *Y;            /* sanitized outcomes (0 where unobserved) */
  unsigned char *ctrl;  /* observed control cell */
  double *S;            /* pairwise sums of squared gaps over shared control cells */
  int *n;               /* pairwise counts */
} prep_t;

static void prep_free(prep_t *q) {
  if (!q) return;
  free(q->Y); free(q->ctrl); free(q->S); free(q->n);
  free(q);
}

static prep_t *prep_new(int N, int T) {
  prep_t *q = (prep_t *) calloc(1, sizeof(prep_t));
  if (!q) return NULL;
  size_t NT = (size_t) N * (size_t) T, NN = (size_t) N * (size_t) N;
  q->N = N; q->T = T;
  q->Y = (double *) malloc(NT * sizeof(double));
  q->ctrl = (unsigned char *) malloc(NT);
  q->S = (double *) malloc(NN * sizeof(double));
  q->n = (int *) malloc(NN * sizeof(int));
  if (!q->Y || !q->ctrl || !q->S || !q->n) { prep_free(q); return NULL; }
  return q;
}

/* Fill the outcomes and the control mask. If necessary, calculate the sums
 * for the pairwise distances. */
static void prep_fill(prep_t *q, const double *Y, const int *W, const int *obs,
                      const int *rows, int need_dist) {
  int N = q->N, T = q->T;
  for (int t = 0; t < T; t++) {
    for (int i = 0; i < N; i++) {
      int src = rows ? rows[i] : i;
      size_t s = (size_t) src + (size_t) N * t, d = (size_t) i + (size_t) N * t;
      int ob = obs[s] && isfinite(Y[s]);
      q->Y[d] = ob ? Y[s] : 0.0;
      q->ctrl[d] = (unsigned char) (ob && W[s] == 0);
    }
  }
  if (!need_dist) return;
  for (int i = 0; i < N; i++) {
    q->S[i + (size_t) N * i] = 0.0;
    q->n[i + (size_t) N * i] = 0;
    for (int j = i + 1; j < N; j++) {
      double s = 0.0;
      int c = 0;
      for (int t = 0; t < T; t++) {
        size_t a = (size_t) i + (size_t) N * t, b = (size_t) j + (size_t) N * t;
        if (q->ctrl[a] && q->ctrl[b]) {
          double d = q->Y[a] - q->Y[b];
          s += d * d;
          c++;
        }
      }
      q->S[i + (size_t) N * j] = q->S[j + (size_t) N * i] = s;
      q->n[i + (size_t) N * j] = q->n[j + (size_t) N * i] = c;
    }
  }
}

/* Weights for the target cell (i, t), from Eq. 3 of the paper. The time
 * weight is exp(-lt |t - s|). The unit weight is exp(-lu d_{-t}(j, i)). The
 * RMS gap uses the periods in which both units are observed controls, but not
 * period t. */
static void local_weights(const prep_t *q, int i, int t, double lt, double lu,
                          int exclude_target, double *w, double *theta, double *omega) {
  int N = q->N, T = q->T;
  for (int s = 0; s < T; s++) theta[s] = (lt == 0.0) ? 1.0 : exp(-lt * fabs((double) (t - s)));
  if (lu == 0.0) {
    for (int j = 0; j < N; j++) omega[j] = 1.0;
  } else {
    size_t it = (size_t) i + (size_t) N * t;
    for (int j = 0; j < N; j++) {
      if (j == i) { omega[j] = 1.0; continue; }
      double S = q->S[i + (size_t) N * j];
      int c = q->n[i + (size_t) N * j];
      size_t jt = (size_t) j + (size_t) N * t;
      if (q->ctrl[it] && q->ctrl[jt]) {
        double d = q->Y[it] - q->Y[jt];
        S -= d * d;
        c -= 1;
        if (S < 0.0) S = 0.0;
      }
      omega[j] = (c > 0) ? exp(-lu * sqrt(S / (double) c)) : 0.0;
    }
  }
  for (int s = 0; s < T; s++) {
    double th = theta[s];
    double *ws_ = w + (size_t) N * s;
    const unsigned char *cs = q->ctrl + (size_t) N * s;
    for (int j = 0; j < N; j++) ws_[j] = cs[j] ? omega[j] * th : 0.0;
  }
  if (exclude_target) w[(size_t) i + (size_t) N * t] = 0.0;
}

/* alpha_i + beta_t is identified iff unit i and period t share a connected
 * component of the bipartite graph of positively weighted control cells. */
static int estimable(ws_t *ws, const double *w, int i, int t) {
  int N = ws->N, T = ws->T;
  int has_row = 0, has_col = 0;
  for (int s = 0; s < T && !has_row; s++) if (w[(size_t) i + (size_t) N * s] > 0.0) has_row = 1;
  for (int j = 0; j < N && !has_col; j++) if (w[(size_t) j + (size_t) N * t] > 0.0) has_col = 1;
  if (!has_row || !has_col) return 0;
  memset(ws->seen, 0, (size_t) (N + T));
  int head = 0, tail = 0;
  ws->queue[tail++] = N + t;
  ws->seen[N + t] = 1;
  while (head < tail) {
    int v = ws->queue[head++];
    if (v >= N) {
      int s = v - N;
      const double *wc = w + (size_t) N * s;
      for (int j = 0; j < N; j++) {
        if (!ws->seen[j] && wc[j] > 0.0) {
          if (j == i) return 1;
          ws->seen[j] = 1;
          ws->queue[tail++] = j;
        }
      }
    } else {
      for (int s = 0; s < T; s++) {
        if (!ws->seen[N + s] && w[(size_t) v + (size_t) N * s] > 0.0) {
          ws->seen[N + s] = 1;
          ws->queue[tail++] = N + s;
        }
      }
    }
  }
  return ws->seen[i];
}

/* ------------------------------------------------------------------------ */
/* Local estimator                                                           */
/* ------------------------------------------------------------------------ */

/* Serial or parallel local estimate on a prepared panel. The function
 * processes the cells of the same unit in time order. Thus, the warm starts
 * are deterministic. */
static int local_estimate_prep(const prep_t *q, const int *W, const int *obs,
                               const int *rows, double lt, double lu, double ln,
                               const trop_opts *o, int parallel, double *tau,
                               double *att, int *n_used, int *n_skipped, int *n_nonconv) {
  int N = q->N, T = q->T;
  size_t NT = (size_t) N * (size_t) T;
  for (size_t c = 0; c < NT; c++) tau[c] = NAN;
  int err = 0, skipped = 0, nonconv = 0;
  int nth = parallel ? nthreads(o) : 1;
#ifdef _OPENMP
#pragma omp parallel num_threads(nth) reduction(+:skipped, nonconv)
#endif
  {
    trop_qos_enter(o);
    ws_t *ws = ws_new(N, T);
    double *theta = (double *) malloc((size_t) T * sizeof(double));
    double *omega = (double *) malloc((size_t) N * sizeof(double));
    if (!ws || !theta || !omega) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
      err = 1;
    }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (int i = 0; i < N; i++) {
      if (err) continue;
      int started = 0;
      int src = rows ? rows[i] : i;
      for (int t = 0; t < T; t++) {
        size_t s = (size_t) src + (size_t) N * t;
        if (W[s] != 1 || !obs[s]) continue;
        size_t c = (size_t) i + (size_t) N * t;
        if (!isfinite(q->Y[c]) || !(obs[s])) continue;
        local_weights(q, i, t, lt, lu, 0, ws->w, theta, omega);
        if (!estimable(ws, ws->w, i, t)) { skipped++; continue; }
        if (!started || !o->warm) { ws_cold(ws); started = 1; }
        int it = 0;
        if (!fit_any(ws, q->Y, ws->w, ln, o, &it)) nonconv++;
        tau[c] = q->Y[c] - ws->alpha[i] - ws->beta[t] - ws->L[c];
      }
    }
    ws_free(ws);
    free(theta);
    free(omega);
    trop_qos_leave(o);
  }
  if (err) return TROP_ERR_ALLOC;
  double sum = 0.0;
  int used = 0;
  for (size_t c = 0; c < NT; c++) if (isfinite(tau[c])) { sum += tau[c]; used++; }
  *att = (used > 0) ? sum / (double) used : NAN;
  if (n_used) *n_used = used;
  if (n_skipped) *n_skipped = skipped;
  if (n_nonconv) *n_nonconv = nonconv;
  return TROP_OK;
}

int trop_local_estimate(const trop_panel *p, double lam_time, double lam_unit,
                        double lam_nn, const trop_opts *opt, double *tau,
                        double *att, int *n_used, int *n_skipped, int *n_nonconv) {
  prep_t *q = prep_new(p->N, p->T);
  if (!q) return TROP_ERR_ALLOC;
  prep_fill(q, p->Y, p->W, p->obs, NULL, lam_unit != 0.0);
  int rc = local_estimate_prep(q, p->W, p->obs, NULL, lam_time, lam_unit, lam_nn,
                               opt, 1, tau, att, n_used, n_skipped, n_nonconv);
  prep_free(q);
  return rc;
}

static int cmp_int(const void *a, const void *b) {
  int x = *(const int *) a, y = *(const int *) b;
  return (x > y) - (x < y);
}

static int local_loocv_prep(const prep_t *q, double lt, double lu, double ln,
                            const int *cells_in, int n_cells_in, const trop_opts *o,
                            double *Q, int *n_valid, int *n_nonconv) {
  int N = q->N, T = q->T;
  size_t NT = (size_t) N * (size_t) T;
  int n_cells = 0;
  int *cells = NULL;
  if (n_cells_in > 0) {
    cells = (int *) malloc((size_t) n_cells_in * sizeof(int));
    if (!cells) return TROP_ERR_ALLOC;
    for (int k = 0; k < n_cells_in; k++) {
      int c = cells_in[k];
      if (c >= 0 && (size_t) c < NT && q->ctrl[c]) cells[n_cells++] = c;
    }
  } else {
    cells = (int *) malloc(NT * sizeof(int));
    if (!cells) return TROP_ERR_ALLOC;
    for (size_t c = 0; c < NT; c++) if (q->ctrl[c]) cells[n_cells++] = (int) c;
  }
  /* Order by unit, then time, so each unit's cells chain warm starts. */
  int *key = (int *) malloc((size_t) (n_cells > 0 ? n_cells : 1) * sizeof(int));
  double *tau = (double *) malloc((size_t) (n_cells > 0 ? n_cells : 1) * sizeof(double));
  int *ustart = (int *) calloc((size_t) N + 1, sizeof(int));
  if (!key || !tau || !ustart) { free(cells); free(key); free(tau); free(ustart); return TROP_ERR_ALLOC; }
  for (int k = 0; k < n_cells; k++) {
    int c = cells[k];
    int i = c % N, t = c / N;
    key[k] = i * T + t;
  }
  qsort(key, (size_t) n_cells, sizeof(int), cmp_int);
  for (int k = 0; k < n_cells; k++) ustart[key[k] / T + 1]++;
  for (int i = 0; i < N; i++) ustart[i + 1] += ustart[i];
  int err = 0, nonconv = 0;
  int nth = nthreads(o);
#ifdef _OPENMP
#pragma omp parallel num_threads(nth) reduction(+:nonconv)
#endif
  {
    trop_qos_enter(o);
    ws_t *ws = ws_new(N, T);
    double *theta = (double *) malloc((size_t) T * sizeof(double));
    double *omega = (double *) malloc((size_t) N * sizeof(double));
    if (!ws || !theta || !omega) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
      err = 1;
    }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (int i = 0; i < N; i++) {
      if (err) continue;
      int started = 0;
      for (int k = ustart[i]; k < ustart[i + 1]; k++) {
        int t = key[k] % T;
        size_t c = (size_t) i + (size_t) N * t;
        local_weights(q, i, t, lt, lu, 1, ws->w, theta, omega);
        if (!estimable(ws, ws->w, i, t)) { tau[k] = NAN; continue; }
        if (!started || !o->warm) { ws_cold(ws); started = 1; }
        int it = 0;
        if (!fit_any(ws, q->Y, ws->w, ln, o, &it)) nonconv++;
        tau[k] = q->Y[c] - ws->alpha[i] - ws->beta[t] - ws->L[c];
      }
    }
    ws_free(ws);
    free(theta);
    free(omega);
    trop_qos_leave(o);
  }
  int rc = TROP_OK;
  if (err) {
    rc = TROP_ERR_ALLOC;
  } else {
    double sum = 0.0;
    int valid = 0;
    for (int k = 0; k < n_cells; k++) if (isfinite(tau[k])) { sum += tau[k] * tau[k]; valid++; }
    *Q = (valid > 0) ? sum : INFINITY;
    if (n_valid) *n_valid = valid;
    if (n_nonconv) *n_nonconv = nonconv;
  }
  free(cells); free(key); free(tau); free(ustart);
  return rc;
}

int trop_local_loocv(const trop_panel *p, double lam_time, double lam_unit,
                     double lam_nn, const int *cells, int n_cells,
                     const trop_opts *opt, double *Q, int *n_valid, int *n_nonconv) {
  prep_t *q = prep_new(p->N, p->T);
  if (!q) return TROP_ERR_ALLOC;
  prep_fill(q, p->Y, p->W, p->obs, NULL, lam_unit != 0.0);
  int rc = local_loocv_prep(q, lam_time, lam_unit, lam_nn, cells, n_cells, opt,
                            Q, n_valid, n_nonconv);
  prep_free(q);
  return rc;
}

/* ------------------------------------------------------------------------ */
/* Search helpers                                                            */
/* ------------------------------------------------------------------------ */

typedef struct {
  double *v; /* rows of (time, unit, nn, score) */
  int n, cap;
} cache_t;

static int same(double a, double b) {
  if (isinf(a) || isinf(b)) return isinf(a) && isinf(b) && ((a > 0) == (b > 0));
  return a == b;
}

static int cache_get(const cache_t *c, double lt, double lu, double ln, double *score) {
  for (int k = 0; k < c->n; k++) {
    const double *r = c->v + 4 * (size_t) k;
    if (same(r[0], lt) && same(r[1], lu) && same(r[2], ln)) { *score = r[3]; return 1; }
  }
  return 0;
}

static int cache_put(cache_t *c, double lt, double lu, double ln, double score) {
  if (c->n == c->cap) {
    int cap = c->cap ? 2 * c->cap : 64;
    double *v = (double *) realloc(c->v, (size_t) cap * 4 * sizeof(double));
    if (!v) return TROP_ERR_ALLOC;
    c->v = v;
    c->cap = cap;
  }
  double *r = c->v + 4 * (size_t) c->n;
  r[0] = lt; r[1] = lu; r[2] = ln; r[3] = score;
  c->n++;
  return TROP_OK;
}

static void cache_trace(const cache_t *c, double *trace, int max_trace, int *n_evals) {
  if (n_evals) *n_evals = c->n;
  if (!trace) return;
  int m = c->n < max_trace ? c->n : max_trace;
  for (int k = 0; k < m; k++)
    for (int j = 0; j < 4; j++) trace[k + (size_t) max_trace * j] = c->v[4 * (size_t) k + j];
}

/* ------------------------------------------------------------------------ */
/* Local search (paper footnote 2)                                           */
/* ------------------------------------------------------------------------ */

typedef struct {
  const prep_t *q;
  const int *cells;
  int n_cells;
  const trop_opts *o;
  cache_t cache;
} lsearch_t;

static int lsearch_eval(lsearch_t *s, double lt, double lu, double ln, double *score) {
  if (cache_get(&s->cache, lt, lu, ln, score)) return TROP_OK;
  double Q = INFINITY;
  int nv = 0, nc = 0;
  int rc = local_loocv_prep(s->q, lt, lu, ln, s->cells, s->n_cells, s->o, &Q, &nv, &nc);
  if (rc != TROP_OK) return rc;
  *score = Q;
  return cache_put(&s->cache, lt, lu, ln, Q);
}

/* Minimum along one coordinate. The value of which is 0 for time, 1 for unit,
 * and 2 for nn. For equal scores, the function keeps the first grid value. If
 * all scores are infinite, the function returns the first grid value. */
static int lsearch_line(lsearch_t *s, int which, const double *grid, int ng,
                        double lt, double lu, double ln, trop_interrupt_fn chk,
                        double *best_v, double *best_s) {
  double bv = grid[0], bs = INFINITY;
  for (int g = 0; g < ng; g++) {
    if (chk && chk()) return TROP_ERR_INTERRUPT;
    double a = lt, b = lu, c = ln, sc = INFINITY;
    if (which == 0) a = grid[g]; else if (which == 1) b = grid[g]; else c = grid[g];
    int rc = lsearch_eval(s, a, b, c, &sc);
    if (rc != TROP_OK) return rc;
    if (sc < bs) { bs = sc; bv = grid[g]; }
  }
  *best_v = bv;
  *best_s = bs;
  return TROP_OK;
}

int trop_local_search(const trop_panel *p, const double *g_time, int n_time,
                      const double *g_unit, int n_unit, const double *g_nn, int n_nn,
                      const int *cells, int n_cells, int max_cycles,
                      const trop_opts *opt, trop_interrupt_fn chk,
                      double *best, double *best_score, int *n_evals,
                      double *trace, int max_trace) {
  if (n_time < 1 || n_unit < 1 || n_nn < 1) return TROP_ERR_INPUT;
  int need_dist = 0;
  for (int g = 0; g < n_unit; g++) if (g_unit[g] != 0.0) need_dist = 1;
  prep_t *q = prep_new(p->N, p->T);
  if (!q) return TROP_ERR_ALLOC;
  prep_fill(q, p->Y, p->W, p->obs, NULL, need_dist);
  lsearch_t s = { q, cells, n_cells, opt, { NULL, 0, 0 } };
  double sc, lt0, lu0, ln0, lt, lu, ln, score = INFINITY, prev = INFINITY;
  int rc;
  /* Stage 1: univariate searches with the other components disabled. */
  if ((rc = lsearch_line(&s, 0, g_time, n_time, 0.0, 0.0, INFINITY, chk, &lt0, &sc)) != TROP_OK) goto done;
  if ((rc = lsearch_line(&s, 2, g_nn, n_nn, 0.0, 0.0, INFINITY, chk, &ln0, &sc)) != TROP_OK) goto done;
  if ((rc = lsearch_line(&s, 1, g_unit, n_unit, 0.0, 0.0, INFINITY, chk, &lu0, &sc)) != TROP_OK) goto done;
  lt = lt0; lu = lu0; ln = ln0;
  /* Stage 2: coordinate cycles until the score is stable. */
  for (int cyc = 0; cyc < max_cycles; cyc++) {
    if ((rc = lsearch_line(&s, 1, g_unit, n_unit, lt, lu, ln, chk, &lu, &sc)) != TROP_OK) goto done;
    if ((rc = lsearch_line(&s, 0, g_time, n_time, lt, lu, ln, chk, &lt, &sc)) != TROP_OK) goto done;
    if ((rc = lsearch_line(&s, 2, g_nn, n_nn, lt, lu, ln, chk, &ln, &score)) != TROP_OK) goto done;
    if (fabs(score - prev) < 1e-6) break;
    prev = score;
  }
  if ((rc = lsearch_eval(&s, lt, lu, ln, &score)) != TROP_OK) goto done;
  best[0] = lt; best[1] = lu; best[2] = ln;
  *best_score = score;
  cache_trace(&s.cache, trace, max_trace, n_evals);
done:
  free(s.cache.v);
  prep_free(q);
  return rc;
}

/* ------------------------------------------------------------------------ */
/* Global estimator (official block implementation)                          */
/* ------------------------------------------------------------------------ */

static int global_estimate_ws(ws_t *ws, const prep_t *q, const int *W, const int *obs,
                              const int *rows, const int *treated, int n_treated, int k,
                              double lt, double lu, double ln, const trop_opts *o,
                              double *theta, double *omega, double *avg,
                              double *tau, double *att, int *converged) {
  int N = q->N, T = q->T, pre = T - k;
  if (n_treated < 1 || k < 1 || pre < 1) return TROP_ERR_INPUT;
  double center = (double) T - (double) k / 2.0;
  for (int s = 0; s < T; s++) theta[s] = exp(-lt * fabs((double) s - center));
  /* In a block design, the cells before the treatment are untreated. Thus,
   * ctrl shows which of these cells have an outcome. As in diff-diff, the
   * calculation does not use missing outcomes. A unit without an outcome
   * before the treatment gets a weight of zero. */
  for (int s = 0; s < pre; s++) {
    double a = 0.0;
    int c = 0;
    for (int m = 0; m < n_treated; m++) {
      size_t idx = (size_t) treated[m] + (size_t) N * s;
      if (q->ctrl[idx]) { a += q->Y[idx]; c++; }
    }
    avg[s] = (c > 0) ? a / (double) c : NAN;
  }
  for (int j = 0; j < N; j++) {
    double S = 0.0;
    int c = 0;
    for (int s = 0; s < pre; s++) {
      size_t idx = (size_t) j + (size_t) N * s;
      if (!q->ctrl[idx] || !isfinite(avg[s])) continue;
      double d = avg[s] - q->Y[idx];
      S += d * d;
      c++;
    }
    omega[j] = (c > 0) ? exp(-lu * sqrt(S / (double) c)) : 0.0;
  }
  for (int s = 0; s < T; s++) {
    for (int j = 0; j < N; j++) {
      size_t c = (size_t) j + (size_t) N * s;
      ws->w[c] = q->ctrl[c] ? omega[j] * theta[s] : 0.0;
    }
  }
  ws_cold(ws);
  int it = 0;
  int conv = fit_any(ws, q->Y, ws->w, ln, o, &it);
  if (converged) *converged = conv;
  double sum = 0.0;
  int used = 0;
  for (int s = 0; s < T; s++) {
    for (int j = 0; j < N; j++) {
      int src = rows ? rows[j] : j;
      size_t c = (size_t) j + (size_t) N * s, cs = (size_t) src + (size_t) N * s;
      double v = NAN;
      if (W[cs] == 1 && obs[cs]) {
        v = q->Y[c] - ws->alpha[j] - ws->beta[s] - ws->L[c];
        sum += v;
        used++;
      }
      if (tau) tau[c] = v;
    }
  }
  *att = (used > 0) ? sum / (double) used : NAN;
  return TROP_OK;
}

int trop_global_estimate(const trop_panel *p, const int *treated, int n_treated,
                         int k, double lam_time, double lam_unit, double lam_nn,
                         const trop_opts *opt, double *tau, double *att, int *converged) {
  int N = p->N, T = p->T;
  prep_t *q = prep_new(N, T);
  ws_t *ws = ws_new(N, T);
  double *theta = (double *) malloc((size_t) T * sizeof(double));
  double *omega = (double *) malloc((size_t) N * sizeof(double));
  double *avg = (double *) malloc((size_t) T * sizeof(double));
  int rc = TROP_ERR_ALLOC;
  if (q && ws && theta && omega && avg) {
    prep_fill(q, p->Y, p->W, p->obs, NULL, 0);
    rc = global_estimate_ws(ws, q, p->W, p->obs, NULL, treated, n_treated, k, lam_time,
                            lam_unit, lam_nn, opt, theta, omega, avg, tau, att, converged);
  }
  prep_free(q); ws_free(ws); free(theta); free(omega); free(avg);
  return rc;
}

int trop_global_placebo(const trop_panel *pc, int k, const int *sets, int n_sets,
                        int set_size, double lam_time, double lam_unit, double lam_nn,
                        const trop_opts *opt, double *rmse, int *n_valid) {
  int N = pc->N, T = pc->T;
  size_t NT = (size_t) N * (size_t) T;
  double *ate = (double *) malloc((size_t) (n_sets > 0 ? n_sets : 1) * sizeof(double));
  if (!ate) return TROP_ERR_ALLOC;
  int err = 0;
  int nth = nthreads(opt);
#ifdef _OPENMP
#pragma omp parallel num_threads(nth)
#endif
  {
    trop_qos_enter(opt);
    prep_t *q = prep_new(N, T);
    ws_t *ws = ws_new(N, T);
    int *Wb = (int *) malloc(NT * sizeof(int));
    double *theta = (double *) malloc((size_t) T * sizeof(double));
    double *omega = (double *) malloc((size_t) N * sizeof(double));
    double *avg = (double *) malloc((size_t) T * sizeof(double));
    if (!q || !ws || !Wb || !theta || !omega || !avg) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
      err = 1;
    }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
    for (int m = 0; m < n_sets; m++) {
      if (err) continue;
      const int *set = sets + (size_t) m * (size_t) set_size;
      memset(Wb, 0, NT * sizeof(int));
      for (int u = 0; u < set_size; u++)
        for (int s = T - k; s < T; s++) Wb[(size_t) set[u] + (size_t) N * s] = 1;
      prep_fill(q, pc->Y, Wb, pc->obs, NULL, 0);
      double a = NAN;
      int conv = 0;
      if (global_estimate_ws(ws, q, Wb, pc->obs, NULL, set, set_size, k, lam_time, lam_unit,
                             lam_nn, opt, theta, omega, avg, NULL, &a, &conv) != TROP_OK) a = NAN;
      ate[m] = a;
    }
    prep_free(q); ws_free(ws); free(Wb); free(theta); free(omega); free(avg);
    trop_qos_leave(opt);
  }
  if (err) { free(ate); return TROP_ERR_ALLOC; }
  double s2 = 0.0;
  int v = 0;
  for (int m = 0; m < n_sets; m++) if (isfinite(ate[m])) { s2 += ate[m] * ate[m]; v++; }
  *rmse = (v > 0) ? sqrt(s2 / (double) v) : INFINITY;
  if (n_valid) *n_valid = v;
  free(ate);
  return TROP_OK;
}

int trop_global_search(const trop_panel *pc, int k, const int *sets, int n_sets,
                       int set_size, const double *g_time, int n_time,
                       const double *g_unit, int n_unit, const double *g_nn, int n_nn,
                       const double *init, int max_cycles, const trop_opts *opt,
                       trop_interrupt_fn chk, double *best, double *best_score,
                       int *n_evals, double *trace, int max_trace, int *fixed_point_out) {
  if (n_time < 1 || n_unit < 1 || n_nn < 1) return TROP_ERR_INPUT;
  cache_t cache = { NULL, 0, 0 };
  double lt = init[0], lu = init[1], ln = init[2];
  int rc = TROP_OK, fixed_point = 0;
  *best_score = INFINITY;
  for (int cyc = 0; cyc < max_cycles && !fixed_point; cyc++) {
    double old_t = lt, old_u = lu, old_n = ln;
    for (int which = 0; which < 3; which++) {
      const double *grid = (which == 0) ? g_unit : (which == 1) ? g_time : g_nn;
      int ng = (which == 0) ? n_unit : (which == 1) ? n_time : n_nn;
      double bv = grid[0], bs = INFINITY;
      for (int g = 0; g < ng; g++) {
        if (chk && chk()) { rc = TROP_ERR_INTERRUPT; goto done; }
        double a = lt, b = lu, c = ln, sc = INFINITY;
        if (which == 0) b = grid[g]; else if (which == 1) a = grid[g]; else c = grid[g];
        if (!cache_get(&cache, a, b, c, &sc)) {
          int nv = 0;
          rc = trop_global_placebo(pc, k, sets, n_sets, set_size, a, b, c, opt, &sc, &nv);
          if (rc != TROP_OK) goto done;
          if ((rc = cache_put(&cache, a, b, c, sc)) != TROP_OK) goto done;
        }
        if (sc < bs) { bs = sc; bv = grid[g]; }
      }
      if (which == 0) lu = bv; else if (which == 1) lt = bv; else { ln = bv; *best_score = bs; }
    }
    fixed_point = same(lt, old_t) && same(lu, old_u) && same(ln, old_n);
  }
  best[0] = lt; best[1] = lu; best[2] = ln;
  cache_trace(&cache, trace, max_trace, n_evals);
  if (fixed_point_out) *fixed_point_out = fixed_point;
done:
  free(cache.v);
  return rc;
}

/* ------------------------------------------------------------------------ */
/* Bootstrap (paper Algorithm 3)                                             */
/* ------------------------------------------------------------------------ */

int trop_bootstrap(const trop_panel *p, int method, int k, double lam_time,
                   double lam_unit, double lam_nn, const int *rows, int B,
                   const trop_opts *opt, trop_interrupt_fn chk, double *draws) {
  int N = p->N, T = p->T;
  int nth = nthreads(opt);
  int chunk = 32 * nth;
  int err = 0;
  for (int start = 0; start < B && !err; start += chunk) {
    int stop = start + chunk < B ? start + chunk : B;
#ifdef _OPENMP
#pragma omp parallel num_threads(nth)
#endif
    {
      trop_qos_enter(opt);
      prep_t *q = prep_new(N, T);
      ws_t *ws = ws_new(N, T);
      double *theta = (double *) malloc((size_t) T * sizeof(double));
      double *omega = (double *) malloc((size_t) N * sizeof(double));
      double *avg = (double *) malloc((size_t) T * sizeof(double));
      double *tau = (double *) malloc((size_t) N * (size_t) T * sizeof(double));
      int *Wb = (int *) malloc((size_t) N * (size_t) T * sizeof(int));
      int *ob = (int *) malloc((size_t) N * (size_t) T * sizeof(int));
      int *tr = (int *) malloc((size_t) N * sizeof(int));
      if (!q || !ws || !theta || !omega || !avg || !tau || !Wb || !ob || !tr) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
        err = 1;
      }
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 1)
#endif
      for (int b = start; b < stop; b++) {
        if (err) continue;
        const int *rb = rows + (size_t) b * (size_t) N;
        int nt = 0;
        for (int i = 0; i < N; i++) {
          int any = 0;
          for (int t = 0; t < T; t++) {
            size_t s = (size_t) rb[i] + (size_t) N * t, d = (size_t) i + (size_t) N * t;
            Wb[d] = p->W[s];
            ob[d] = p->obs[s];
            if (p->W[s] == 1) any = 1;
          }
          if (any) tr[nt++] = i;
        }
        prep_fill(q, p->Y, p->W, p->obs, rb, method == 0 && lam_unit != 0.0);
        double a = NAN;
        int rc;
        if (method == 0) {
          int nu = 0, ns = 0, nc = 0;
          rc = local_estimate_prep(q, Wb, ob, NULL, lam_time, lam_unit, lam_nn, opt, 0,
                                   tau, &a, &nu, &ns, &nc);
        } else {
          int conv = 0;
          rc = (nt > 0) ? global_estimate_ws(ws, q, Wb, ob, NULL, tr, nt, k, lam_time,
                                             lam_unit, lam_nn, opt, theta, omega, avg,
                                             NULL, &a, &conv)
                        : TROP_ERR_INPUT;
        }
        draws[b] = (rc == TROP_OK) ? a : NAN;
      }
      prep_free(q); ws_free(ws); free(theta); free(omega); free(avg); free(tau);
      free(Wb); free(ob); free(tr);
      trop_qos_leave(opt);
    }
    if (!err && chk && chk()) return TROP_ERR_INTERRUPT;
  }
  return err ? TROP_ERR_ALLOC : TROP_OK;
}
