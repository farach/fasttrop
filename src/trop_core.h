/*
 * fasttrop: C core for the Triply Robust Panel (TROP) estimator of Athey,
 * Imbens, Qu, and Viviano (arXiv:2508.21536, Journal of Applied Econometrics,
 * 2026).
 *
 * This header and trop_core.c use only standard C, OpenMP, BLAS, and LAPACK.
 * They do not call the R API. Thus, the same source code compiles into the R
 * package and into the C program of the benchmark in inst/bench.
 *
 * The matrices are column-major, with units in rows: X[i + N * t].
 */
#ifndef FASTTROP_CORE_H
#define FASTTROP_CORE_H

typedef struct {
  int max_iter;   /* maximum outer iterations (prox: max_iter * max_inner steps) */
  int max_inner;  /* proximal (FISTA) steps in each outer iteration */
  double tol;     /* tolerance for the largest absolute change */
  int threads;    /* OpenMP threads, with 1 for values less than 1 */
  int warm;       /* 1: start from the previous fit of the unit, 0: start from zero */
  int solver;     /* 0: accelerated proximal gradient, 1: block coordinate method */
  int throttle;   /* 0: no Windows power throttling for C threads, 1: Windows policy */
} trop_opts;

typedef struct {
  int N, T;
  const double *Y; /* outcomes, not used where obs == 0 */
  const int *W;    /* treatment indicator, 0 or 1 */
  const int *obs;  /* 1 if the cell has an outcome */
} trop_panel;

typedef int (*trop_interrupt_fn)(void);

enum { TROP_OK = 0, TROP_ERR_ALLOC = 1, TROP_ERR_INTERRUPT = 2, TROP_ERR_INPUT = 3 };

/* Find the minimum of sum w (Y - a_i - b_t - L)^2 + lambda_nn ||L||_* for the
 * given weights. With lambda_nn = +Inf, there is no low-rank term (L = 0). The
 * function returns TROP_OK or an error code. *converged is 1 if the fit got to
 * the tolerance. */
int trop_fit_weighted(int N, int T, const double *Y, const double *w,
                      double lambda_nn, const trop_opts *opt,
                      double *alpha, double *beta, double *L,
                      int *iters, int *converged);

/* Local estimator (Eq. 12 and Algorithm 2 of the paper) with fixed tuning
 * parameters. tau (length N*T) gets the effect for each treated cell with an
 * outcome. The other cells get NaN. A treated cell also gets NaN if the
 * untreated cells do not identify its unit effect and period effect. */
int trop_local_estimate(const trop_panel *p, double lam_time, double lam_unit,
                        double lam_nn, const trop_opts *opt, double *tau,
                        double *att, int *n_used, int *n_skipped, int *n_nonconv);

/* Local leave-one-out criterion Q(lambda) of Eq. 4 and Eq. 5 for untreated
 * cells. The array cells contains n_cells linear indices (i + N*t) of
 * untreated cells with an outcome. If n_cells is 0, the function uses all of
 * these cells. */
int trop_local_loocv(const trop_panel *p, double lam_time, double lam_unit,
                     double lam_nn, const int *cells, int n_cells,
                     const trop_opts *opt, double *Q, int *n_valid, int *n_nonconv);

/* Search of footnote 2 of the paper for the local estimator: first a search on
 * each axis, then coordinate cycles. best = (time, unit, nn). If trace is not
 * NULL, it gets a maximum of max_trace rows of (time, unit, nn, Q). */
int trop_local_search(const trop_panel *p, const double *g_time, int n_time,
                      const double *g_unit, int n_unit, const double *g_nn, int n_nn,
                      const int *cells, int n_cells, int max_cycles,
                      const trop_opts *opt, trop_interrupt_fn chk,
                      double *best, double *best_score, int *n_evals,
                      double *trace, int max_trace);

/* Global block estimator of the Python package from the authors. It
 * calculates one weighted fit. The time distance is the distance to the center
 * of the treated block. The unit distance is the distance to the average
 * treated unit in the periods before the treatment. The array treated contains
 * the row indices of the treated units, and k is the number of treated periods
 * at the end. */
int trop_global_estimate(const trop_panel *p, const int *treated, int n_treated,
                         int k, double lam_time, double lam_unit, double lam_nn,
                         const trop_opts *opt, double *tau, double *att, int *converged);

/* Placebo cross-validation for the global estimator on a panel with only
 * untreated units. The array sets contains n_sets placebo sets of set_size row
 * indices, in row-major order. The function gives the RMSE of the placebo
 * estimates. */
int trop_global_placebo(const trop_panel *pc, int k, const int *sets, int n_sets,
                        int set_size, double lam_time, double lam_unit, double lam_nn,
                        const trop_opts *opt, double *rmse, int *n_valid);

/* Coordinate search for the global estimator, as in the Python package from
 * the authors. Each cycle changes unit, then time, then nn. The search stops
 * when a cycle does not change the three values. *fixed_point is 1 if this
 * occurs in max_cycles cycles. */
int trop_global_search(const trop_panel *pc, int k, const int *sets, int n_sets,
                       int set_size, const double *g_time, int n_time,
                       const double *g_unit, int n_unit, const double *g_nn, int n_nn,
                       const double *init, int max_cycles, const trop_opts *opt,
                       trop_interrupt_fn chk, double *best, double *best_score,
                       int *n_evals, double *trace, int max_trace, int *fixed_point);

/* Unit bootstrap with fixed tuning parameters (Algorithm 3 of the paper). The
 * array rows is a B x N matrix in row-major order with resampled unit indices
 * that start at 0. Method 0 is local and method 1 is global, with k treated
 * periods. The array draws gets B estimates, with NaN for a failure. */
int trop_bootstrap(const trop_panel *p, int method, int k, double lam_time,
                   double lam_unit, double lam_nn, const int *rows, int B,
                   const trop_opts *opt, trop_interrupt_fn chk, double *draws);

/* Wall-clock time in seconds from a high-resolution monotonic timer, in
 * trop_os.c. */
double trop_wtime(void);

/* Put these calls before and after the C work on a thread. On Windows, if
 * opt->throttle is 0, the thread has no execution-speed power throttling until
 * the paired call to trop_qos_leave(). The calls can nest. On other
 * operating systems, the functions have no effect. The functions are in
 * trop_os.c. */
void trop_qos_enter(const trop_opts *opt);
void trop_qos_leave(const trop_opts *opt);

#endif
