/*
 * fasttrop: native core for the Triply Robust Panel (TROP) estimator of
 * Athey, Imbens, Qu and Viviano (arXiv:2508.21536; J. Applied Econometrics 2026).
 *
 * This header and trop_core.c use only standard C, OpenMP, and BLAS/LAPACK.
 * They make no R API calls, so the same source builds into the R package and
 * into the standalone benchmark driver in inst/bench.
 *
 * Matrices are column-major with units in rows: X[i + N * t].
 */
#ifndef FASTTROP_CORE_H
#define FASTTROP_CORE_H

typedef struct {
  int max_iter;   /* outer iterations (block) or FISTA steps / max_inner (prox) */
  int max_inner;  /* proximal (FISTA) iterations per outer iteration */
  double tol;     /* max-abs change tolerance for convergence */
  int threads;    /* OpenMP threads; values < 1 mean 1 */
  int warm;       /* 1: chain warm starts within a unit; 0: cold start every fit */
  int solver;     /* 0: profiled accelerated proximal gradient; 1: reference block scheme */
  int throttle;   /* 0: opt kernel threads out of Windows power throttling; 1: allow it */
} trop_opts;

typedef struct {
  int N, T;
  const double *Y; /* outcomes; entries with obs == 0 are ignored */
  const int *W;    /* treatment indicator, 0/1 */
  const int *obs;  /* 1 if the outcome is observed */
} trop_panel;

typedef int (*trop_interrupt_fn)(void);

enum { TROP_OK = 0, TROP_ERR_ALLOC = 1, TROP_ERR_INTERRUPT = 2, TROP_ERR_INPUT = 3 };

/* Fit min sum w (Y - a_i - b_t - L)^2 + lambda_nn ||L||_* for given weights.
 * lambda_nn = +Inf disables the low-rank term (L = 0). Returns TROP_OK or an
 * error code; *converged reports whether the tolerance was met. */
int trop_fit_weighted(int N, int T, const double *Y, const double *w,
                      double lambda_nn, const trop_opts *opt,
                      double *alpha, double *beta, double *L,
                      int *iters, int *converged);

/* Local estimator (paper Eq. 12, Algorithm 2) at fixed tuning parameters.
 * tau (length N*T) receives cell effects for treated observed cells and NaN
 * elsewhere (also NaN for non-estimable cells). */
int trop_local_estimate(const trop_panel *p, double lam_time, double lam_unit,
                        double lam_nn, const trop_opts *opt, double *tau,
                        double *att, int *n_used, int *n_skipped, int *n_nonconv);

/* Local leave-one-out criterion Q(lambda) (Eq. 4-5) over control cells.
 * cells holds n_cells linear indices (i + N*t) of control observations;
 * n_cells = 0 means every observed control cell. */
int trop_local_loocv(const trop_panel *p, double lam_time, double lam_unit,
                     double lam_nn, const int *cells, int n_cells,
                     const trop_opts *opt, double *Q, int *n_valid, int *n_nonconv);

/* Footnote-2 search for the local estimator: univariate initial searches,
 * then coordinate cycling. best = (time, unit, nn). trace (optional) receives
 * up to max_trace rows of (time, unit, nn, Q). */
int trop_local_search(const trop_panel *p, const double *g_time, int n_time,
                      const double *g_unit, int n_unit, const double *g_nn, int n_nn,
                      const int *cells, int n_cells, int max_cycles,
                      const trop_opts *opt, trop_interrupt_fn chk,
                      double *best, double *best_score, int *n_evals,
                      double *trace, int max_trace);

/* Global block estimator (official reference implementation): one weighted
 * fit with time distance to the treated-block centre and unit distance to the
 * average treated trajectory over pre-periods. treated lists treated row
 * indices; k is the number of final treated periods. */
int trop_global_estimate(const trop_panel *p, const int *treated, int n_treated,
                         int k, double lam_time, double lam_unit, double lam_nn,
                         const trop_opts *opt, double *tau, double *att, int *converged);

/* Placebo cross-validation for the global estimator on a control-only panel.
 * sets holds n_sets placebo sets of set_size row indices (row-major). Returns
 * the RMSE of placebo ATEs. */
int trop_global_placebo(const trop_panel *pc, int k, const int *sets, int n_sets,
                        int set_size, double lam_time, double lam_unit, double lam_nn,
                        const trop_opts *opt, double *rmse, int *n_valid);

/* Coordinate-cycling search for the global estimator (official cycle):
 * updates unit, then time, then nn until the triplet stops changing.
 * *fixed_point is 1 if the cycle converged within max_cycles. */
int trop_global_search(const trop_panel *pc, int k, const int *sets, int n_sets,
                       int set_size, const double *g_time, int n_time,
                       const double *g_unit, int n_unit, const double *g_nn, int n_nn,
                       const double *init, int max_cycles, const trop_opts *opt,
                       trop_interrupt_fn chk, double *best, double *best_score,
                       int *n_evals, double *trace, int max_trace, int *fixed_point);

/* Unit bootstrap at fixed tuning parameters (paper Algorithm 3). rows is a
 * B x N row-major matrix of resampled unit indices (0-based). method 0 =
 * local, 1 = global (k = treated periods). draws receives B estimates (NaN on
 * failure). */
int trop_bootstrap(const trop_panel *p, int method, int k, double lam_time,
                   double lam_unit, double lam_nn, const int *rows, int B,
                   const trop_opts *opt, trop_interrupt_fn chk, double *draws);

/* Wall-clock seconds from a high-resolution monotonic timer (trop_os.c). */
double trop_wtime(void);

/* Bracket native work on the calling thread: on Windows, unless
 * opt->throttle is set, the thread is opted out of execution-speed power
 * throttling until the matching leave call (calls nest). No-ops elsewhere.
 * Defined in trop_os.c. */
void trop_qos_enter(const trop_opts *opt);
void trop_qos_leave(const trop_opts *opt);

#endif
