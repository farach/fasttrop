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
  /* Fault test of the SVD, for the package tests only. The value 0 gives the
   * usual calculation. Bit 1 makes each dgesdd call fail, and bit 2 also makes
   * each dgesvd call fail. With bit 4, the faults occur only in bootstrap draws
   * with an odd index. */
  int fault;
} trop_opts;

/* Counts of the fits of one calculation. */
typedef struct {
  int used;     /* fits with a finite result */
  int skipped;  /* cells that the untreated cells do not identify */
  int nonconv;  /* fits with a finite result that stopped at the iteration limit */
  int failed;   /* fits with a failed SVD or a result that is not finite */
} trop_counts;

/* Columns of a trace of a tuning search: time, unit, nn, score, valid fits,
 * fits that did not converge, and fits that failed. */
#define TROP_TRACE_COLS 7

typedef struct {
  int N, T;
  const double *Y; /* outcomes, not used where obs == 0 */
  const int *W;    /* treatment indicator, 0 or 1 */
  const int *obs;  /* 1 if the cell has an outcome */
  /* Optional covariates (Eq. 13 of the paper). X has N*T*p values, with
   * covariate k at X[i + N*t + N*T*k]. If p is 0, X is not used. */
  const double *X;
  int p;
  /* Bootstrap only. 0: fit X*phi with each weighted fit (joint). 1: before
   * each draw, remove X*phi from one unweighted two-way fit on the untreated
   * cells (residualize), then fit without covariates. */
  int resid;
} trop_panel;

typedef int (*trop_interrupt_fn)(void);

enum { TROP_OK = 0, TROP_ERR_ALLOC = 1, TROP_ERR_INTERRUPT = 2, TROP_ERR_INPUT = 3,
       TROP_ERR_NUMERIC = 4 };

/* Find the minimum of sum w (Y - a_i - b_t - L)^2 + lambda_nn ||L||_* for the
 * given weights. With lambda_nn = +Inf, there is no low-rank term (L = 0). The
 * function returns TROP_OK or an error code. The code TROP_ERR_NUMERIC shows
 * that the SVD failed. *converged is 1 if the fit got to the tolerance. */
int trop_fit_weighted(int N, int T, const double *Y, const double *w,
                      double lambda_nn, const trop_opts *opt,
                      double *alpha, double *beta, double *L,
                      int *iters, int *converged);

/* Local estimator (Eq. 12 and Algorithm 2 of the paper) with fixed tuning
 * parameters. tau (length N*T) gets the effect for each treated cell with an
 * outcome. The other cells get NaN. A treated cell also gets NaN if the
 * untreated cells do not identify its unit effect and period effect. A cell
 * with a failed fit also gets NaN. The function puts the counts of the fits
 * in cnt. */
int trop_local_estimate(const trop_panel *p, double lam_time, double lam_unit,
                        double lam_nn, const trop_opts *opt, double *tau,
                        double *att, trop_counts *cnt, double *phi);

/* phi (length p, can be NULL) gets the covariate coefficients. For the local
 * estimator, it is the mean of the coefficients of the treated-cell fits. */

/* Unweighted two-way fit of Y on unit effects, period effects, and X over the
 * untreated cells with an outcome. phi (length p) gets the coefficients. A
 * column that the fixed effects absorb gets 0. */
int trop_residualize(const trop_panel *p, const trop_opts *opt, double *phi);

/* Local leave-one-out criterion Q(lambda) of Eq. 4 and Eq. 5 for untreated
 * cells. The array cells contains n_cells linear indices (i + N*t) of
 * untreated cells with an outcome. If n_cells is 0, the function uses all of
 * these cells. A cell with a failed fit is not in Q. cnt->used is the number
 * of cells in Q. */
int trop_local_loocv(const trop_panel *p, double lam_time, double lam_unit,
                     double lam_nn, const int *cells, int n_cells,
                     const trop_opts *opt, double *Q, trop_counts *cnt);

/* Search of footnote 2 of the paper for the local estimator: first a search on
 * each axis, then coordinate cycles. best = (time, unit, nn). If trace is not
 * NULL, it gets a maximum of max_trace rows with TROP_TRACE_COLS columns. */
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
 * at the end.
 *
 * The value cnt->used is the number of treated cells with an effect. A treated
 * cell is skipped if the cells with a positive weight do not identify its unit
 * effect and its period effect. The values cnt->nonconv and cnt->failed are 0
 * or 1, for the one fit. */
int trop_global_estimate(const trop_panel *p, const int *treated, int n_treated,
                         int k, double lam_time, double lam_unit, double lam_nn,
                         const trop_opts *opt, double *tau, double *att, trop_counts *cnt,
                         double *phi);

/* Placebo cross-validation for the global estimator on a panel with only
 * untreated units. The array sets contains n_sets placebo sets of set_size row
 * indices, in row-major order. The function gives the RMSE of the placebo
 * estimates of the sets. The count cnt->used is the number of placebo cells
 * with an effect in the sets with an estimate. The count cnt->skipped is the
 * number of placebo cells that the untreated cells do not identify. */
int trop_global_placebo(const trop_panel *pc, int k, const int *sets, int n_sets,
                        int set_size, double lam_time, double lam_unit, double lam_nn,
                        const trop_opts *opt, double *rmse, trop_counts *cnt);

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
 * periods. The array draws gets B estimates, with NaN for a failure.
 *
 * If diag is not NULL, it gets a B x 5 column-major matrix with the counts of
 * each draw: used, skipped, nonconv, failed, and dropped. A dropped cell is a
 * skipped cell of the draw that has eligible = 1. The array eligible has N x T
 * values, indexed by the original unit. It can be NULL.
 *
 * Group means (local method only): gid has n_layers layers of N x T cells in
 * column-major order, indexed by the original unit. Each value is a group
 * from 0 to n_groups - 1, or -1 for no group. A cell can be in one group of
 * each layer. The array gdraws gets the weighted mean effect of each group in
 * a B x n_groups column-major matrix. A group with no cell gets NaN. Use
 * gid = NULL for no groups.
 *
 * The weight of a cell is uw[i] for its original unit i, or 1 if uw is NULL.
 * If per_unit is 1, the function divides this weight by the number of cells
 * of the resampled row in the group. Then each row has the total weight uw[i]
 * in each group that it is in. */
int trop_bootstrap(const trop_panel *p, int method, int k, double lam_time,
                   double lam_unit, double lam_nn, const int *rows, int B,
                   const trop_opts *opt, trop_interrupt_fn chk, double *draws,
                   int *diag, const int *eligible, const int *gid, int n_layers,
                   int n_groups, const double *uw, int per_unit, double *gdraws);

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
