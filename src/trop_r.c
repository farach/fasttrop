/*
 * R interface for the fasttrop native core. Each entry point checks inputs,
 * passes raw pointers to trop_core.c, and reports the native kernel time so
 * callers can separate kernel cost from R-side overhead.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <math.h>
#include <string.h>
#include "trop_core.h"

/* ---- helpers ---------------------------------------------------------- */

static void chk_fn(void *data) { (void) data; R_CheckUserInterrupt(); }
static int r_interrupt(void) { return !R_ToplevelExec(chk_fn, NULL); }

static void fail_rc(int rc) {
  if (rc == TROP_ERR_ALLOC) Rf_error("fasttrop: memory allocation failed in the native core.");
  if (rc == TROP_ERR_INTERRUPT) Rf_error("fasttrop: interrupted by the user.");
  if (rc == TROP_ERR_INPUT) Rf_error("fasttrop: invalid input passed to the native core.");
  if (rc != TROP_OK) Rf_error("fasttrop: native core error code %d.", rc);
}

static double list_num(SEXP lst, const char *name, double def) {
  SEXP nms = Rf_getAttrib(lst, R_NamesSymbol);
  for (R_xlen_t k = 0; k < XLENGTH(lst); k++) {
    if (strcmp(CHAR(STRING_ELT(nms, k)), name) == 0) return Rf_asReal(VECTOR_ELT(lst, k));
  }
  return def;
}

static trop_opts read_opts(SEXP control) {
  if (!Rf_isNewList(control)) Rf_error("`control` must be a list.");
  trop_opts o;
  o.max_iter = (int) list_num(control, "max_iter", 100);
  o.max_inner = (int) list_num(control, "max_inner", 10);
  o.tol = list_num(control, "tol", 1e-6);
  o.threads = (int) list_num(control, "threads", 1);
  o.warm = (int) list_num(control, "warm", 1);
  o.solver = (int) list_num(control, "solver", 0);
  o.throttle = (int) list_num(control, "throttle", 0);
  if (o.max_iter < 1 || o.max_inner < 1 || !(o.tol > 0)) Rf_error("Invalid control settings.");
  return o;
}

static void check_matrix(SEXP x, int type, int N, int T, const char *nm) {
  if (TYPEOF(x) != type) Rf_error("`%s` has the wrong storage type.", nm);
  if (XLENGTH(x) != (R_xlen_t) N * (R_xlen_t) T) Rf_error("`%s` must have N x T entries.", nm);
}

static trop_panel read_panel(SEXP Y, SEXP W, SEXP obs) {
  if (!Rf_isMatrix(Y) || TYPEOF(Y) != REALSXP) Rf_error("`Y` must be a double matrix.");
  SEXP dim = Rf_getAttrib(Y, R_DimSymbol);
  trop_panel p;
  p.N = INTEGER(dim)[0];
  p.T = INTEGER(dim)[1];
  check_matrix(W, INTSXP, p.N, p.T, "W");
  check_matrix(obs, INTSXP, p.N, p.T, "obs");
  p.Y = REAL(Y);
  p.W = INTEGER(W);
  p.obs = INTEGER(obs);
  return p;
}

static SEXP named_list(int n, const char **names) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, n));
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, n));
  for (int k = 0; k < n; k++) SET_STRING_ELT(nm, k, Rf_mkChar(names[k]));
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(2);
  return out;
}

static SEXP lambda3(const double *v) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 3));
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, 3));
  REAL(out)[0] = v[0]; REAL(out)[1] = v[1]; REAL(out)[2] = v[2];
  SET_STRING_ELT(nm, 0, Rf_mkChar("time"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("unit"));
  SET_STRING_ELT(nm, 2, Rf_mkChar("nn"));
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(2);
  return out;
}

static SEXP trace_matrix(const double *tr, int n, int max_trace) {
  int m = n < max_trace ? n : max_trace;
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, m, 4));
  for (int j = 0; j < 4; j++)
    for (int k = 0; k < m; k++) REAL(out)[k + (R_xlen_t) m * j] = tr[k + (R_xlen_t) max_trace * j];
  UNPROTECT(1);
  return out;
}

static void read_lambda(SEXP lambda, double *lt, double *lu, double *ln) {
  if (TYPEOF(lambda) != REALSXP || XLENGTH(lambda) != 3) Rf_error("`lambda` must be a double vector of length 3.");
  *lt = REAL(lambda)[0]; *lu = REAL(lambda)[1]; *ln = REAL(lambda)[2];
  if (!(*lt >= 0) || !(*lu >= 0) || !(*ln >= 0) || isinf(*lt) || isinf(*lu))
    Rf_error("`lambda` must have finite non-negative time and unit values and a non-negative nn value.");
}

/* ---- entry points ------------------------------------------------------ */

SEXP C_fit_weighted(SEXP Y, SEXP w, SEXP lambda_nn, SEXP control) {
  if (!Rf_isMatrix(Y) || TYPEOF(Y) != REALSXP) Rf_error("`Y` must be a double matrix.");
  SEXP dim = Rf_getAttrib(Y, R_DimSymbol);
  int N = INTEGER(dim)[0], T = INTEGER(dim)[1];
  check_matrix(w, REALSXP, N, T, "w");
  trop_opts o = read_opts(control);
  SEXP alpha = PROTECT(Rf_allocVector(REALSXP, N));
  SEXP beta = PROTECT(Rf_allocVector(REALSXP, T));
  SEXP L = PROTECT(Rf_allocMatrix(REALSXP, N, T));
  int iters = 0, conv = 0;
  trop_qos_enter(&o);
  int rc = trop_fit_weighted(N, T, REAL(Y), REAL(w), Rf_asReal(lambda_nn), &o,
                             REAL(alpha), REAL(beta), REAL(L), &iters, &conv);
  trop_qos_leave(&o);
  if (rc != TROP_OK) { UNPROTECT(3); fail_rc(rc); }
  const char *nm[] = {"alpha", "beta", "L", "iterations", "converged"};
  SEXP out = PROTECT(named_list(5, nm));
  SET_VECTOR_ELT(out, 0, alpha);
  SET_VECTOR_ELT(out, 1, beta);
  SET_VECTOR_ELT(out, 2, L);
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(iters));
  SET_VECTOR_ELT(out, 4, Rf_ScalarLogical(conv));
  UNPROTECT(4);
  return out;
}

SEXP C_local_estimate(SEXP Y, SEXP W, SEXP obs, SEXP lambda, SEXP control) {
  trop_panel p = read_panel(Y, W, obs);
  trop_opts o = read_opts(control);
  double lt, lu, ln;
  read_lambda(lambda, &lt, &lu, &ln);
  SEXP tau = PROTECT(Rf_allocMatrix(REALSXP, p.N, p.T));
  double att = NAN, t0 = trop_wtime();
  int used = 0, skipped = 0, nonconv = 0;
  trop_qos_enter(&o);
  int rc = trop_local_estimate(&p, lt, lu, ln, &o, REAL(tau), &att, &used, &skipped, &nonconv);
  trop_qos_leave(&o);
  double el = trop_wtime() - t0;
  if (rc != TROP_OK) { UNPROTECT(1); fail_rc(rc); }
  const char *nm[] = {"att", "tau", "n_used", "n_skipped", "n_nonconverged", "kernel_seconds"};
  SEXP out = PROTECT(named_list(6, nm));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(att));
  SET_VECTOR_ELT(out, 1, tau);
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(used));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(skipped));
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(nonconv));
  SET_VECTOR_ELT(out, 5, Rf_ScalarReal(el));
  UNPROTECT(2);
  return out;
}

SEXP C_local_loocv(SEXP Y, SEXP W, SEXP obs, SEXP lambda, SEXP cells, SEXP control) {
  trop_panel p = read_panel(Y, W, obs);
  trop_opts o = read_opts(control);
  double lt, lu, ln;
  read_lambda(lambda, &lt, &lu, &ln);
  int n_cells = 0;
  const int *cp = NULL;
  if (!Rf_isNull(cells)) {
    if (TYPEOF(cells) != INTSXP) Rf_error("`cells` must be an integer vector.");
    n_cells = (int) XLENGTH(cells);
    cp = INTEGER(cells);
  }
  double Q = NAN, t0 = trop_wtime();
  int nv = 0, nc = 0;
  trop_qos_enter(&o);
  int rc = trop_local_loocv(&p, lt, lu, ln, cp, n_cells, &o, &Q, &nv, &nc);
  trop_qos_leave(&o);
  double el = trop_wtime() - t0;
  fail_rc(rc);
  const char *nm[] = {"Q", "n_valid", "n_nonconverged", "kernel_seconds"};
  SEXP out = PROTECT(named_list(4, nm));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(Q));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(nv));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(nc));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(el));
  UNPROTECT(1);
  return out;
}

SEXP C_local_search(SEXP Y, SEXP W, SEXP obs, SEXP g_time, SEXP g_unit, SEXP g_nn,
                    SEXP cells, SEXP max_cycles, SEXP control) {
  trop_panel p = read_panel(Y, W, obs);
  trop_opts o = read_opts(control);
  if (TYPEOF(g_time) != REALSXP || TYPEOF(g_unit) != REALSXP || TYPEOF(g_nn) != REALSXP)
    Rf_error("Grids must be double vectors.");
  int n_cells = 0;
  const int *cp = NULL;
  if (!Rf_isNull(cells)) {
    if (TYPEOF(cells) != INTSXP) Rf_error("`cells` must be an integer vector.");
    n_cells = (int) XLENGTH(cells);
    cp = INTEGER(cells);
  }
  int nt = (int) XLENGTH(g_time), nu = (int) XLENGTH(g_unit), nn = (int) XLENGTH(g_nn);
  int max_trace = 3 * (nt + nu + nn) * (Rf_asInteger(max_cycles) + 2) + 8;
  double *trace = (double *) R_alloc((size_t) max_trace * 4, sizeof(double));
  double best[3] = {NAN, NAN, NAN}, score = NAN, t0 = trop_wtime();
  int n_evals = 0;
  trop_qos_enter(&o);
  int rc = trop_local_search(&p, REAL(g_time), nt, REAL(g_unit), nu, REAL(g_nn), nn, cp, n_cells,
                             Rf_asInteger(max_cycles), &o, r_interrupt, best, &score, &n_evals,
                             trace, max_trace);
  trop_qos_leave(&o);
  double el = trop_wtime() - t0;
  fail_rc(rc);
  const char *nm[] = {"lambda", "score", "n_evals", "trace", "kernel_seconds"};
  SEXP out = PROTECT(named_list(5, nm));
  SET_VECTOR_ELT(out, 0, lambda3(best));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(score));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(n_evals));
  SET_VECTOR_ELT(out, 3, trace_matrix(trace, n_evals, max_trace));
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal(el));
  UNPROTECT(1);
  return out;
}

SEXP C_global_estimate(SEXP Y, SEXP W, SEXP obs, SEXP treated, SEXP k, SEXP lambda, SEXP control) {
  trop_panel p = read_panel(Y, W, obs);
  trop_opts o = read_opts(control);
  double lt, lu, ln;
  read_lambda(lambda, &lt, &lu, &ln);
  if (TYPEOF(treated) != INTSXP || XLENGTH(treated) < 1) Rf_error("`treated` must be a non-empty integer vector.");
  SEXP tau = PROTECT(Rf_allocMatrix(REALSXP, p.N, p.T));
  double att = NAN, t0 = trop_wtime();
  int conv = 0;
  trop_qos_enter(&o);
  int rc = trop_global_estimate(&p, INTEGER(treated), (int) XLENGTH(treated), Rf_asInteger(k),
                                lt, lu, ln, &o, REAL(tau), &att, &conv);
  trop_qos_leave(&o);
  double el = trop_wtime() - t0;
  if (rc != TROP_OK) { UNPROTECT(1); fail_rc(rc); }
  const char *nm[] = {"att", "tau", "converged", "kernel_seconds"};
  SEXP out = PROTECT(named_list(4, nm));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(att));
  SET_VECTOR_ELT(out, 1, tau);
  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(conv));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(el));
  UNPROTECT(2);
  return out;
}

static trop_panel control_panel(SEXP Yc, SEXP obs_c, SEXP Wzero) {
  return read_panel(Yc, Wzero, obs_c);
}

SEXP C_global_placebo(SEXP Yc, SEXP obs_c, SEXP Wzero, SEXP k, SEXP sets, SEXP lambda, SEXP control) {
  trop_panel pc = control_panel(Yc, obs_c, Wzero);
  trop_opts o = read_opts(control);
  double lt, lu, ln;
  read_lambda(lambda, &lt, &lu, &ln);
  if (!Rf_isMatrix(sets) || TYPEOF(sets) != INTSXP) Rf_error("`sets` must be an integer matrix.");
  SEXP sd = Rf_getAttrib(sets, R_DimSymbol);
  double rmse = NAN, t0 = trop_wtime();
  int nv = 0;
  trop_qos_enter(&o);
  int rc = trop_global_placebo(&pc, Rf_asInteger(k), INTEGER(sets), INTEGER(sd)[0], INTEGER(sd)[1],
                               lt, lu, ln, &o, &rmse, &nv);
  trop_qos_leave(&o);
  double el = trop_wtime() - t0;
  fail_rc(rc);
  const char *nm[] = {"rmse", "n_valid", "kernel_seconds"};
  SEXP out = PROTECT(named_list(3, nm));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(rmse));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(nv));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(el));
  UNPROTECT(1);
  return out;
}

SEXP C_global_search(SEXP Yc, SEXP obs_c, SEXP Wzero, SEXP k, SEXP sets, SEXP g_time,
                     SEXP g_unit, SEXP g_nn, SEXP init, SEXP max_cycles, SEXP control) {
  trop_panel pc = control_panel(Yc, obs_c, Wzero);
  trop_opts o = read_opts(control);
  if (!Rf_isMatrix(sets) || TYPEOF(sets) != INTSXP) Rf_error("`sets` must be an integer matrix.");
  if (TYPEOF(init) != REALSXP || XLENGTH(init) != 3) Rf_error("`init` must be a double vector of length 3.");
  SEXP sd = Rf_getAttrib(sets, R_DimSymbol);
  int nt = (int) XLENGTH(g_time), nu = (int) XLENGTH(g_unit), nn = (int) XLENGTH(g_nn);
  int mc = Rf_asInteger(max_cycles);
  int max_trace = (nt + nu + nn) * (mc + 1) + 8;
  double *trace = (double *) R_alloc((size_t) max_trace * 4, sizeof(double));
  double best[3] = {NAN, NAN, NAN}, score = NAN, t0 = trop_wtime();
  int n_evals = 0, fixed_point = 0;
  trop_qos_enter(&o);
  int rc = trop_global_search(&pc, Rf_asInteger(k), INTEGER(sets), INTEGER(sd)[0], INTEGER(sd)[1],
                              REAL(g_time), nt, REAL(g_unit), nu, REAL(g_nn), nn, REAL(init), mc, &o,
                              r_interrupt, best, &score, &n_evals, trace, max_trace, &fixed_point);
  trop_qos_leave(&o);
  double el = trop_wtime() - t0;
  fail_rc(rc);
  const char *nm[] = {"lambda", "score", "n_evals", "trace", "fixed_point", "kernel_seconds"};
  SEXP out = PROTECT(named_list(6, nm));
  SET_VECTOR_ELT(out, 0, lambda3(best));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(score));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(n_evals));
  SET_VECTOR_ELT(out, 3, trace_matrix(trace, n_evals, max_trace));
  SET_VECTOR_ELT(out, 4, Rf_ScalarLogical(fixed_point));
  SET_VECTOR_ELT(out, 5, Rf_ScalarReal(el));
  UNPROTECT(1);
  return out;
}

SEXP C_bootstrap(SEXP Y, SEXP W, SEXP obs, SEXP method, SEXP k, SEXP lambda, SEXP rows, SEXP control) {
  trop_panel p = read_panel(Y, W, obs);
  trop_opts o = read_opts(control);
  double lt, lu, ln;
  read_lambda(lambda, &lt, &lu, &ln);
  if (!Rf_isMatrix(rows) || TYPEOF(rows) != INTSXP) Rf_error("`rows` must be an integer matrix.");
  SEXP rd = Rf_getAttrib(rows, R_DimSymbol);
  int B = INTEGER(rd)[0];
  if (INTEGER(rd)[1] != p.N) Rf_error("`rows` must have N columns.");
  /* rows arrive as an R (column-major) B x N matrix; the core wants row-major. */
  int *rm = (int *) R_alloc((size_t) B * (size_t) p.N, sizeof(int));
  const int *ri = INTEGER(rows);
  for (int b = 0; b < B; b++)
    for (int i = 0; i < p.N; i++) {
      int v = ri[b + (R_xlen_t) B * i];
      if (v < 0 || v >= p.N) Rf_error("`rows` holds an out-of-range unit index.");
      rm[(size_t) b * (size_t) p.N + (size_t) i] = v;
    }
  SEXP draws = PROTECT(Rf_allocVector(REALSXP, B));
  double t0 = trop_wtime();
  trop_qos_enter(&o);
  int rc = trop_bootstrap(&p, Rf_asInteger(method), Rf_asInteger(k), lt, lu, ln, rm, B, &o,
                          r_interrupt, REAL(draws));
  trop_qos_leave(&o);
  double el = trop_wtime() - t0;
  if (rc != TROP_OK) { UNPROTECT(1); fail_rc(rc); }
  const char *nm[] = {"draws", "kernel_seconds"};
  SEXP out = PROTECT(named_list(2, nm));
  SET_VECTOR_ELT(out, 0, draws);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(el));
  UNPROTECT(2);
  return out;
}

SEXP C_openmp_threads(void) {
#ifdef _OPENMP
  return Rf_ScalarLogical(1);
#else
  return Rf_ScalarLogical(0);
#endif
}

static const R_CallMethodDef call_methods[] = {
  {"C_fit_weighted", (DL_FUNC) &C_fit_weighted, 4},
  {"C_local_estimate", (DL_FUNC) &C_local_estimate, 5},
  {"C_local_loocv", (DL_FUNC) &C_local_loocv, 6},
  {"C_local_search", (DL_FUNC) &C_local_search, 9},
  {"C_global_estimate", (DL_FUNC) &C_global_estimate, 7},
  {"C_global_placebo", (DL_FUNC) &C_global_placebo, 7},
  {"C_global_search", (DL_FUNC) &C_global_search, 11},
  {"C_bootstrap", (DL_FUNC) &C_bootstrap, 8},
  {"C_openmp_threads", (DL_FUNC) &C_openmp_threads, 0},
  {NULL, NULL, 0}
};

void R_init_fasttrop(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
