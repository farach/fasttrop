/*
 * Standalone C driver for the fasttrop parity benchmark.
 *
 * It compiles the package core (src/trop_core.c) into this translation unit
 * and runs the same sequence of native calls as trop_run() in R: tuning
 * search, point estimate, and bootstrap. Build it with R's compiler and
 * flags, together with src/trop_os.c, and link R's own BLAS/LAPACK (see
 * run_benchmark.R) so the only difference from the R path is the R entry
 * layer.
 *
 * Usage: trop_driver <plan.bin> <reps> <result.bin>
 */
#include "trop_core.c"
#include <stdio.h>

static FILE *fin;

static int rd_i(void) {
  int v;
  if (fread(&v, sizeof(int), 1, fin) != 1) { fprintf(stderr, "read error\n"); exit(2); }
  return v;
}
static double rd_d(void) {
  double v;
  if (fread(&v, sizeof(double), 1, fin) != 1) { fprintf(stderr, "read error\n"); exit(2); }
  return v;
}
static int *rd_iv(size_t n) {
  int *v = (int *) malloc((n ? n : 1) * sizeof(int));
  if (n && fread(v, sizeof(int), n, fin) != n) { fprintf(stderr, "read error\n"); exit(2); }
  return v;
}
static double *rd_dv(size_t n) {
  double *v = (double *) malloc((n ? n : 1) * sizeof(double));
  if (n && fread(v, sizeof(double), n, fin) != n) { fprintf(stderr, "read error\n"); exit(2); }
  return v;
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: trop_driver plan.bin reps result.bin\n"); return 1; }
  int reps = atoi(argv[2]);
  fin = fopen(argv[1], "rb");
  if (!fin) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
  if (rd_i() != 0x54524F50) { fprintf(stderr, "bad magic\n"); return 1; }
  int N = rd_i(), T = rd_i(), method = rd_i(), k = rd_i();
  trop_opts o;
  o.threads = rd_i(); o.max_iter = rd_i(); o.max_inner = rd_i(); o.solver = rd_i();
  o.warm = rd_i();
  int max_cycles = rd_i();
  o.throttle = rd_i();
  o.tol = rd_d();
  int nt = rd_i(), nu = rd_i(), nn = rd_i();
  double *gt = rd_dv((size_t) nt), *gu = rd_dv((size_t) nu), *gn = rd_dv((size_t) nn);
  int has_lambda = rd_i();
  double lam[3];
  lam[0] = rd_d(); lam[1] = rd_d(); lam[2] = rd_d();
  int n_cells = rd_i();
  int *cells = rd_iv((size_t) n_cells);
  int Nc = rd_i();
  double *Yc = rd_dv((size_t) Nc * (size_t) T);
  int *oc = rd_iv((size_t) Nc * (size_t) T);
  int n_sets = rd_i(), set_size = rd_i();
  int *sets = rd_iv((size_t) n_sets * (size_t) set_size);
  double init[3];
  init[0] = rd_d(); init[1] = rd_d(); init[2] = rd_d();
  int n_treated = rd_i();
  int *treated = rd_iv((size_t) n_treated);
  int B = rd_i();
  int *rows = rd_iv((size_t) B * (size_t) N);
  double *Y = rd_dv((size_t) N * (size_t) T);
  int *W = rd_iv((size_t) N * (size_t) T);
  int *obs = rd_iv((size_t) N * (size_t) T);
  fclose(fin);

  trop_panel p = { N, T, Y, W, obs };
  int *Wz = (int *) calloc((size_t) (Nc > 0 ? Nc : 1) * (size_t) T, sizeof(int));
  trop_panel pc = { Nc, T, Yc, Wz, oc };
  double *tau = (double *) malloc((size_t) N * (size_t) T * sizeof(double));
  double *draws = (double *) malloc((size_t) (B > 0 ? B : 1) * sizeof(double));
  double best[3] = { lam[0], lam[1], lam[2] }, score = NAN, att = NAN;
  double t_best = INFINITY, s_best = 0, e_best = 0, b_best = 0;

  /* Same thread power-throttling policy as the R entry points. */
  trop_qos_enter(&o);
  for (int r = 0; r < reps; r++) {
    double t0 = trop_wtime();
    int rc = TROP_OK, n_evals = 0;
    if (!has_lambda) {
      if (method == 0) {
        rc = trop_local_search(&p, gt, nt, gu, nu, gn, nn, n_cells ? cells : NULL, n_cells,
                               max_cycles, &o, NULL, best, &score, &n_evals, NULL, 0);
      } else {
        int fp = 0;
        rc = trop_global_search(&pc, k, sets, n_sets, set_size, gt, nt, gu, nu, gn, nn, init,
                                max_cycles, &o, NULL, best, &score, &n_evals, NULL, 0, &fp);
      }
    }
    if (rc != TROP_OK) { fprintf(stderr, "search failed %d\n", rc); return 3; }
    double t1 = trop_wtime();
    if (method == 0) {
      int nu_ = 0, ns = 0, nc = 0;
      rc = trop_local_estimate(&p, best[0], best[1], best[2], &o, tau, &att, &nu_, &ns, &nc);
    } else {
      int conv = 0;
      rc = trop_global_estimate(&p, treated, n_treated, k, best[0], best[1], best[2], &o, tau,
                                &att, &conv);
    }
    if (rc != TROP_OK) { fprintf(stderr, "estimate failed %d\n", rc); return 4; }
    double t2 = trop_wtime();
    if (B > 0) {
      rc = trop_bootstrap(&p, method, k, best[0], best[1], best[2], rows, B, &o, NULL, draws);
      if (rc != TROP_OK) { fprintf(stderr, "bootstrap failed %d\n", rc); return 5; }
    }
    double t3 = trop_wtime();
    if (t3 - t0 < t_best) { t_best = t3 - t0; s_best = t1 - t0; e_best = t2 - t1; b_best = t3 - t2; }
  }
  trop_qos_leave(&o);

  FILE *fo = fopen(argv[3], "wb");
  double head[8] = { best[0], best[1], best[2], att, t_best, s_best, e_best, b_best };
  fwrite(head, sizeof(double), 8, fo);
  fwrite(draws, sizeof(double), (size_t) B, fo);
  fclose(fo);
  printf("lambda %.10g %.10g %.10g att %.15g total %.6f search %.6f estimate %.6f boot %.6f\n",
         best[0], best[1], best[2], att, t_best, s_best, e_best, b_best);
  return 0;
}
