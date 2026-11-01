"""Reference results for the covariate oracle tests of fasttrop.

The paper (Section 6.2, Eq. 13) adds covariates to the working model:

    tau_hat = argmin_tau min_{alpha, beta, phi, L}
              sum_{j,s} theta_s omega_j (Y_js - alpha_j - beta_s - X_js phi - L_js - tau W_js)^2
              + lambda_nn ||L||_*

phi has no penalty. The weights (Eq. 3) do not change and come from the raw
outcomes Y, not from X. The packages trop (zhaonanq/TROP-Estimator, MIT) and
diff-diff (igerber/diff-diff, MIT) do not have covariates. Thus, this script
calculates these references at the commits in oracle_meta.json:

* Reference A (global method for a block design only): a copy of
  TROP_TWFE_average of trop with only one change, the term X phi in the
  objective (trop_twfe_average_cov). With p = 0, the script makes sure that
  this copy gives the ATT of the original function.
* Reference B (local and global method): a direct CVXPY/CLARABEL solve of
  Eq. 13 for given weights. CLARABEL stops at a relative gap of about 1e-9,
  which leaves errors of about 1e-5 in the fitted values (the objective is
  flat near the optimum). Thus, for a finite lambda_nn, a proximal gradient
  refinement on L (polish()) starts from the CLARABEL solution. The script
  makes sure that this gives the objective and the fitted values of diff-diff
  on the weights of oracle_solver.csv (p = 0). The weights come from an
  implementation of the conventions of diff-diff. The script makes sure that
  these weights are equal to the weights of diff-diff on the fixture panels
  D1-D4.
* A numpy weighted least squares fit with dummy variables for lambda_nn = Inf,
  as a check of Reference B.
* A "residualize once" coefficient phi_res: OLS of y on unit dummies, time
  dummies, and X over the untreated cells (no weights).

The objective of Reference B is sum w r^2 + lambda_nn ||L||_*, as in
oracle_fixtures.py (objective()), diff-diff, and trop. L is a free variable
in all cells, also in the cells with zero weight (the treated cells). The cell
effect is tau_it = Y_it - alpha_i - beta_t - X_it phi - L_it, and the ATT is
the mean of tau_it over the treated cells.

Panels (written to tests/testthat/fixtures/, long format, columns
unit,time,y,d,x1,x2; seeds of numpy.random.default_rng in PANELS below):
* XC1  Block design, two time-varying covariates.
* XC2  Staggered adoption, two time-varying covariates.
* XC3  Block design. x2 is constant within each unit, thus collinear with the
       unit fixed effects. phi_2 is not identified (NaN in the outputs), but
       tau is identified. The script checks that tau does not change if x2 is
       removed.
DGP: y = a_i + b_t + (U V')_it + x1 phi1 + x2 phi2 + e + tau d, with
phi_true = (0.8, -0.5) and tau = 1.

Outputs (tests/testthat/fixtures/):
* oracle_cov.csv        ATT (and phi for the global method) per panel,
                        configuration, method, and reference (A or B).
* oracle_cov_cells.csv  The cell effects of Reference B, with phi and the
                        objective value of each fit (for the local method,
                        each cell has its own fit).
* oracle_cov_resid.csv  phi_res per panel.
* oracle_cov_meta.json  Commits, versions, solver, tolerance, seeds, checks.

Usage, in the package folder, with clones of the two repositories at the
commits of oracle_meta.json:
    python data-raw/oracle_covariates.py <path-to-TROP-Estimator> <path-to-diff-diff>
CVXPY has no wheels for Windows on ARM64. There, use an x64 Python.
"""
import json
import math
import platform
import subprocess
import sys
import time

import numpy as np
import pandas as pd

TROP_SRC, DIFFDIFF_SRC = sys.argv[1], sys.argv[2]
sys.path.insert(0, TROP_SRC + "/src")
sys.path.insert(0, DIFFDIFF_SRC)

import clarabel  # noqa: E402
import cvxpy as cp  # noqa: E402
import scipy  # noqa: E402
import scipy.linalg  # noqa: E402
import trop  # noqa: E402
from diff_diff import TROP  # noqa: E402

FIX = "tests/testthat/fixtures"
TOL = 1e-9
SOLVER = "CLARABEL"
CLARABEL_OPTS = dict(tol_gap_abs=TOL, tol_gap_rel=TOL, tol_feas=TOL, tol_ktratio=TOL,
                     max_iter=1000)
PHI_TRUE = (0.8, -0.5)
TAU_TRUE = 1.0
N_UNITS, N_PERIODS = 12, 10
# Panel name: (seed, design, treated units (1-based), first treated period (1-based) per unit,
#              x2 constant within unit)
PANELS = {
    "XC1": (20261101, "block", [10, 11, 12], [8, 8, 8], False),
    "XC2": (20261102, "staggered", [7, 8, 9, 10, 11, 12], [6, 6, 7, 7, 9, 9], False),
    "XC3": (20261103, "block", [9, 10, 11, 12], [7, 7, 7, 7], True),
}


def commit(path):
    return subprocess.run(["git", "-C", path, "rev-parse", "HEAD"], capture_output=True,
                          text=True, check=True).stdout.strip()


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------
def make_panel(name):
    seed, _, treated, first, x2_const = PANELS[name]
    rng = np.random.default_rng(seed)
    N, T = N_UNITS, N_PERIODS
    a = rng.normal(0.0, 1.0, N)
    b = np.cumsum(rng.normal(0.0, 0.5, T))
    low = 0.5 * rng.normal(size=(N, 2)) @ rng.normal(size=(T, 2)).T
    d = np.zeros((N, T), dtype=int)
    for u, f in zip(treated, first):
        d[u - 1, f - 1:] = 1
    ever = d.any(axis=1).astype(float)
    # x1 depends on the unit effect and on treatment, thus a fit without x1 is biased.
    x1 = 0.5 * a[:, None] + rng.normal(size=(N, T)) + 0.5 * ever[:, None] + 0.5 * d
    if x2_const:
        x2 = np.repeat(rng.normal(size=(N, 1)), T, axis=1)
    else:
        x2 = rng.normal(size=(N, T)) + 0.3 * np.arange(T)[None, :] / T - 0.4 * d
    e = rng.normal(0.0, 0.3, size=(N, T))
    y = a[:, None] + b[None, :] + low + PHI_TRUE[0] * x1 + PHI_TRUE[1] * x2 + e + TAU_TRUE * d
    rows = [dict(unit=i + 1, time=t + 1, y=y[i, t], d=d[i, t], x1=x1[i, t], x2=x2[i, t])
            for t in range(T) for i in range(N)]
    df = pd.DataFrame(rows)
    df.to_csv(f"{FIX}/{name}.csv", index=False, float_format="%.15g")
    return df


def read_panel(name, covariates=("x1", "x2")):
    df = pd.read_csv(f"{FIX}/{name}.csv")
    units = sorted(df["unit"].unique())
    times = sorted(df["time"].unique())
    # The layout of diff-diff: periods in rows and units in columns.
    piv = lambda v: df.pivot(index="time", columns="unit", values=v).reindex(index=times, columns=units).values  # noqa: E731
    Y = piv("y")
    D = np.nan_to_num(piv("d")).astype(int)
    X = np.stack([piv(c) for c in covariates], axis=2) if covariates else np.zeros(Y.shape + (0,))
    return df, units, times, Y, D, X


def block_info(D):
    T, N = D.shape
    treated_units = [i for i in range(N) if D[:, i].sum() > 0]
    first = [int(np.argmax(D[:, i] == 1)) for i in treated_units]
    block = len(set(first)) == 1 and all(D[first[0]:, i].all() for i in treated_units)
    return treated_units, (T - first[0] if block else None)


# ---------------------------------------------------------------------------
# Reference A: copy of TROP_TWFE_average of trop with the term X phi.
# ---------------------------------------------------------------------------
# Copied from zhaonanq/TROP-Estimator, src/trop/estimator.py (TROP_TWFE_average),
# commit in oracle_meta.json. MIT License, Copyright (c) 2026 zhaonanq.
# Changes: the argument X (N x T x p) and the free, unpenalized variable phi in
# the residual; solver_opts is passed to prob.solve; the function returns
# (tau_hat, phi). The docstring and the input checks are shorter. All else is
# the same as in the original.
def trop_twfe_average_cov(Y, W, X, treated_units, lambda_unit, lambda_time, lambda_nn,
                          treated_periods=10, solver=None, verbose=False, solver_opts=None):
    Y = np.asarray(Y, dtype=float)
    W = np.asarray(W, dtype=float)
    X = np.asarray(X, dtype=float)
    N, T = Y.shape
    p = X.shape[2] if X.ndim == 3 else 0
    if treated_periods <= 0 or treated_periods >= T:
        raise ValueError("treated_periods must be in 1..T-1.")
    treated_units_arr = np.asarray(treated_units, dtype=int)

    center = T - treated_periods / 2.0
    dist_time = np.abs(np.arange(T, dtype=float) - center)

    average_treated = np.mean(Y[treated_units_arr, :], axis=0)
    mask = np.ones((N, T), dtype=float)
    mask[:, -treated_periods:] = 0.0
    A = np.sum(((average_treated - Y) ** 2) * mask, axis=1)
    B = np.sum(mask, axis=1)
    dist_unit = np.sqrt(A / B)

    delta_unit = np.exp(-0.5 * lambda_unit * dist_unit)
    delta_time = np.exp(-0.5 * lambda_time * dist_time)
    delta = np.outer(delta_unit, delta_time)

    unit_effects = cp.Variable((1, N))
    time_effects = cp.Variable((1, T))
    mu = cp.Variable()

    unit_factor = cp.kron(np.ones((T, 1)), unit_effects).T
    time_factor = cp.kron(np.ones((N, 1)), time_effects)

    # Covariates (the only change of the model): X phi, phi free.
    if p > 0:
        phi = cp.Variable(p)
        cov_factor = sum(phi[k] * X[:, :, k] for k in range(p))
    else:
        phi = None
        cov_factor = 0.0

    is_low_rank = not math.isinf(float(lambda_nn))
    sqrt_weights = np.multiply(1.0 - W, delta)

    if is_low_rank:
        L = cp.Variable((N, T))
        residual = Y - mu - unit_factor - time_factor - cov_factor - L
        loss = cp.sum_squares(cp.multiply(sqrt_weights, residual)) + float(lambda_nn) * cp.norm(L, "nuc")
        default_solver = "SCS"
    else:
        residual = Y - mu - unit_factor - time_factor - cov_factor
        loss = cp.sum_squares(cp.multiply(sqrt_weights, residual))
        default_solver = "OSQP"

    prob = cp.Problem(cp.Minimize(loss))
    chosen_solver = solver or default_solver
    prob.solve(solver=chosen_solver, verbose=verbose, **(solver_opts or {}))
    if prob.status not in ("optimal", "optimal_inaccurate"):
        raise RuntimeError(f"Optimization did not converge. Solver={chosen_solver}, status={prob.status}.")

    fit = mu.value + np.tile(unit_effects.value.T, (1, T)) + np.tile(time_effects.value, (N, 1))
    if p > 0:
        fit = fit + np.einsum("ntk,k->nt", X, phi.value)
    if is_low_rank:
        fit = fit + L.value
    tau_hat = float(np.mean((Y - fit)[W > 0]))
    if not np.isfinite(tau_hat):
        raise RuntimeError("Optimization did not return a valid tau.")
    return tau_hat, (np.asarray(phi.value, dtype=float) if p > 0 else np.zeros(0))


# ---------------------------------------------------------------------------
# Weights (Eq. 3), conventions of diff-diff (T x N layout).
# ---------------------------------------------------------------------------
def local_weights(Y, D, i, t, lam_time, lam_unit):
    """theta_s omega_j for the treated cell (i, t), before the (1 - W) mask."""
    T, N = Y.shape
    theta = np.exp(-lam_time * np.abs(np.arange(T) - t))
    omega = np.ones(N)
    if lam_unit != 0:
        for j in range(N):
            if j == i:
                continue
            ok = np.ones(T, dtype=bool)
            ok[t] = False
            ok &= (D[:, i] == 0) & (D[:, j] == 0) & np.isfinite(Y[:, i]) & np.isfinite(Y[:, j])
            omega[j] = np.exp(-lam_unit * np.sqrt(np.mean((Y[ok, i] - Y[ok, j]) ** 2))) if ok.any() else 0.0
    return np.outer(theta, omega)


def global_weights(Y, D, lam_time, lam_unit, k):
    """Weights of the global method, with the (1 - W) mask.

    Time: distance of the 0-based period s to T - k/2. Unit: RMS distance of
    unit j to the mean of the treated units (nanmean) over the first T - k
    periods, over the periods where both are finite.
    """
    T, N = Y.shape
    treated = np.where(D.any(axis=0))[0]
    theta = np.exp(-lam_time * np.abs(np.arange(T, dtype=float) - (T - k / 2.0)))
    avg = np.nanmean(Y[:, treated], axis=1)
    pre = np.arange(T) < T - k
    diff = avg[:, None] - Y
    ok = np.isfinite(diff) & pre[:, None]
    cnt = ok.sum(axis=0)
    dist = np.sqrt((np.where(ok, diff, 0.0) ** 2).sum(axis=0) / np.maximum(cnt, 1))
    omega = np.exp(-lam_unit * dist)
    omega[cnt == 0] = 0.0
    return np.outer(theta, omega) * (1 - D)


# ---------------------------------------------------------------------------
# Reference B: Eq. 13 for given weights (CVXPY, CLARABEL).
# ---------------------------------------------------------------------------
def design_matrix(T, N, X):
    tt, ii = np.meshgrid(np.arange(T), np.arange(N), indexing="ij")
    return np.hstack([np.eye(N)[ii.ravel()], np.eye(T)[tt.ravel()], X.reshape(T * N, X.shape[2])])


def objective(Y, w, F, L, lam_nn):
    r = np.where(w > 0, Y - F, 0.0)
    val = float(np.sum(w * r * r))
    return val + (lam_nn * float(np.linalg.svd(L, compute_uv=False).sum()) if math.isfinite(lam_nn) else 0.0)


def solve_eq13_clarabel(Y, w, X, lam_nn):
    """min sum w (Y - alpha_j - beta_s - X phi - L)^2 + lam_nn ||L||_* (T x N layout).

    w must be zero for the treated and the missing cells. Returns the fitted
    values, phi, and L.
    """
    T, N = Y.shape
    p = X.shape[2]
    alpha, beta = cp.Variable(N), cp.Variable(T)
    fit = (np.ones((T, 1)) @ cp.reshape(alpha, (1, N), order="C")
           + cp.reshape(beta, (T, 1), order="C") @ np.ones((1, N)))
    phi = cp.Variable(p) if p > 0 else None
    for c in range(p):
        fit = fit + phi[c] * X[:, :, c]
    finite = math.isfinite(lam_nn)
    if finite:
        L = cp.Variable((T, N))
        fit = fit + L
    Ys = np.where(np.isfinite(Y), Y, 0.0)
    loss = cp.sum_squares(cp.multiply(np.sqrt(w), Ys - fit))
    if finite:
        loss = loss + lam_nn * cp.normNuc(L)
    prob = cp.Problem(cp.Minimize(loss))
    prob.solve(solver=SOLVER, **CLARABEL_OPTS)
    if prob.status != "optimal":
        raise RuntimeError(f"CLARABEL status {prob.status}")
    phi_v = np.asarray(phi.value, dtype=float) if p > 0 else np.zeros(0)
    Lv = L.value if finite else np.zeros((T, N))
    F = alpha.value[None, :] + beta.value[:, None] + np.einsum("tnk,k->tn", X, phi_v) + Lv
    return F, phi_v, Lv


def polish(Y, w, X, lam_nn, L0, max_iter=200000, step_tol=1e-14):
    """Refine the CLARABEL solution of Eq. 13 to high accuracy.

    CLARABEL stops at a relative gap of about 1e-9. The objective is flat near
    the optimum, thus the fitted values can be off by about 1e-5. This function
    profiles out (alpha, beta, phi) by the weighted projection on the dummies
    and X, and runs proximal gradient steps (FISTA with restart) on L from the
    CLARABEL solution. The step is 1/(2 max w), the threshold is
    lambda_nn/(2 max w). It stops when no element of L changes by more than
    step_tol.
    """
    T, N = Y.shape
    Z = design_matrix(T, N, X)
    s = np.sqrt(w.ravel())
    SZ = Z * s[:, None]
    Q, R, _ = scipy.linalg.qr(SZ, mode="economic", pivoting=True)
    Q = Q[:, :int(np.sum(np.abs(np.diag(R)) > 1e-10 * abs(R[0, 0])))]
    Ys = np.where(np.isfinite(Y), Y, 0.0).ravel()

    def resid(L):
        v = s * (Ys - L.ravel())
        return v - Q @ (Q.T @ v)

    def value(L):
        v = resid(L)
        return float(v @ v) + lam_nn * float(np.linalg.svd(L, compute_uv=False).sum())

    def prox_step(M):
        G = M + (s * resid(M)).reshape(T, N) / w.max()
        U, d, Vt = np.linalg.svd(G, full_matrices=False)
        return (U * np.maximum(d - lam_nn / (2.0 * w.max()), 0.0)) @ Vt

    L, L_prev, tk, f_old = L0.copy(), L0.copy(), 1.0, value(L0)
    for it in range(max_iter):
        tn = (1.0 + math.sqrt(1.0 + 4.0 * tk * tk)) / 2.0
        L_new = prox_step(L + (tk - 1.0) / tn * (L - L_prev))
        f_new = value(L_new)
        if f_new > f_old:  # adaptive restart
            tn = 1.0
            L_new = prox_step(L)
            f_new = value(L_new)
        L_prev, L, tk, f_old = L, L_new, tn, f_new
        if np.max(np.abs(L - L_prev)) < step_tol:
            break
    else:
        raise RuntimeError("polish did not converge")
    coef = np.linalg.lstsq(SZ, s * (Ys - L.ravel()), rcond=None)[0]
    F = (Z @ coef).reshape(T, N) + L
    return F, coef[N + T:], L


def solve_eq13(Y, w, X, lam_nn):
    """Reference B: CLARABEL, then polish() if lambda_nn is finite.

    Returns the fitted values, phi, the objective, and the fitted values of
    CLARABEL alone.
    """
    F0, phi0, L0 = solve_eq13_clarabel(Y, w, X, lam_nn)
    if not math.isfinite(lam_nn):
        return F0, phi0, objective(Y, w, F0, L0, lam_nn), F0
    F, phi, L = polish(Y, w, X, lam_nn, L0)
    val, val0 = objective(Y, w, F, L, lam_nn), objective(Y, w, F0, L0, lam_nn)
    assert val <= val0 * (1 + 1e-12), (val, val0)
    return F, phi, val, F0


def solve_wls(Y, w, X):
    """Weighted least squares with dummies (lambda_nn = Inf): fitted values and phi."""
    T, N = Y.shape
    p = X.shape[2]
    design = design_matrix(T, N, X)
    use = w.ravel() > 0
    sw = np.sqrt(w.ravel()[use])
    coef = np.linalg.lstsq(design[use] * sw[:, None], Y.ravel()[use] * sw, rcond=None)[0]
    return (design @ coef).reshape(T, N), coef[N + T:]


def identified_phi(D_mask, X):
    """For each covariate: is phi_k identified in a two-way fit on the cells D_mask?"""
    T, N = D_mask.shape
    p = X.shape[2]
    design = design_matrix(T, N, X)[D_mask.ravel()]
    full = np.linalg.matrix_rank(design)
    return np.array([full - np.linalg.matrix_rank(np.delete(design, N + T + c, axis=1)) == 1
                     for c in range(p)])


def phi_cols(phi, ident):
    out = {}
    for c in range(2):
        out[f"phi_{c + 1}"] = float(phi[c]) if c < len(phi) and ident[c] else np.nan
    return out


# ---------------------------------------------------------------------------
# Checks on the fixture panels without covariates.
# ---------------------------------------------------------------------------
lambdas = pd.read_csv(f"{FIX}/lambdas.csv")
lambdas["nn"] = lambdas["nn"].replace("Inf", np.inf).astype(float)

pinned = json.load(open(f"{FIX}/oracle_meta.json"))
assert commit(TROP_SRC) == pinned["trop_commit"], "trop is not at the commit of oracle_meta.json"
assert commit(DIFFDIFF_SRC) == pinned["diff_diff_commit"], "diff-diff is not at the commit of oracle_meta.json"

checks = dict(weights_local_vs_diffdiff=0.0, weights_global_vs_diffdiff=0.0,
              A_p0_vs_trop=0.0, A_xzero_vs_trop=0.0)
dd = TROP(method="local")
for name in ["D1", "D2", "D3", "D4"]:
    df, units, times, Y, D, X0 = read_panel(name, covariates=())
    T, N = Y.shape
    treated_units, k = block_info(D)
    ctrl_idx = np.where(~np.any(D == 1, axis=0))[0]
    for _, lam in lambdas.iterrows():
        lt, lu, ln = float(lam["time"]), float(lam["unit"]), float(lam["nn"])
        for t in range(T):
            for i in range(N):
                if D[t, i] == 1:
                    ref = dd._compute_observation_weights(Y, D, i, t, lt, lu, ctrl_idx, N, T)
                    checks["weights_local_vs_diffdiff"] = max(
                        checks["weights_local_vs_diffdiff"],
                        float(np.max(np.abs(local_weights(Y, D, i, t, lt, lu) - ref))))
        if k is None:
            continue
        g = TROP(method="global")
        ref = g._compute_global_weights(Y, D, lt, lu, k, N, T)
        checks["weights_global_vs_diffdiff"] = max(
            checks["weights_global_vs_diffdiff"], float(np.max(np.abs(global_weights(Y, D, lt, lu, k) - ref))))
        # Reference A with p = 0 (same problem as the original) and with X = 0 (p = 2).
        orig = trop.TROP_TWFE_average(Y.T, D.T, treated_units, lu, lt, ln, treated_periods=k, solver=SOLVER)
        a0, _ = trop_twfe_average_cov(Y.T, D.T, np.zeros((N, T, 0)), treated_units, lu, lt, ln,
                                      treated_periods=k, solver=SOLVER)
        az, _ = trop_twfe_average_cov(Y.T, D.T, np.zeros((N, T, 2)), treated_units, lu, lt, ln,
                                      treated_periods=k, solver=SOLVER)
        checks["A_p0_vs_trop"] = max(checks["A_p0_vs_trop"], abs(a0 - orig))
        checks["A_xzero_vs_trop"] = max(checks["A_xzero_vs_trop"], abs(az - orig))
    print(name, "checks done", flush=True)
assert checks["weights_local_vs_diffdiff"] < 1e-14, checks
assert checks["weights_global_vs_diffdiff"] < 1e-14, checks
assert checks["A_p0_vs_trop"] < 1e-8, checks

# Reference B with p = 0 on the weights of oracle_solver.csv: objective and
# fitted values of diff-diff (diff-diff uses a tolerance of 1e-9).
checks.update(B_p0_vs_diffdiff_objective_rel=0.0, B_p0_vs_diffdiff_fitted=0.0,
              B_p0_clarabel_only_vs_diffdiff_fitted=0.0)
sol = pd.read_csv(f"{FIX}/oracle_solver.csv")
obj = pd.read_csv(f"{FIX}/oracle_objective.csv")
for _, o in obj.iterrows():
    df, units, times, Y, D, X0 = read_panel(o["dataset"], covariates=())
    s = sol[(sol.dataset == o["dataset"]) & (sol.config == o["config"]) & (sol.method == o["method"])]
    w = s.pivot(index="time", columns="unit", values="w").reindex(index=times, columns=units).values
    Fd = s.pivot(index="time", columns="unit", values="fitted").reindex(index=times, columns=units).values
    F, _, val, F0 = solve_eq13(Y, w, X0, float(o["lambda_nn"]))
    fin = np.isfinite(Y)
    checks["B_p0_vs_diffdiff_objective_rel"] = max(checks["B_p0_vs_diffdiff_objective_rel"],
                                                   (val - o["objective"]) / o["objective"])
    checks["B_p0_vs_diffdiff_fitted"] = max(checks["B_p0_vs_diffdiff_fitted"], float(np.max(np.abs(F - Fd)[fin])))
    checks["B_p0_clarabel_only_vs_diffdiff_fitted"] = max(checks["B_p0_clarabel_only_vs_diffdiff_fitted"],
                                                          float(np.max(np.abs(F0 - Fd)[fin])))
print("B (p = 0) vs diff-diff done", flush=True)
assert checks["B_p0_vs_diffdiff_objective_rel"] < 1e-12, checks
assert checks["B_p0_vs_diffdiff_fitted"] < 1e-6, checks

# ---------------------------------------------------------------------------
# Covariate panels.
# ---------------------------------------------------------------------------
att_rows, cell_rows, resid_rows = [], [], []
checks.update(A_vs_B_global_att=0.0, A_vs_B_global_phi=0.0, A_vs_B_global_att_finite_nn=0.0,
              A_vs_B_global_att_inf_nn=0.0, A_vs_B_clarabel_only_global_att=0.0,
              B_clarabel_only_vs_polished_tau=0.0, wls_vs_B_tau=0.0, wls_vs_B_phi=0.0,
              XC3_tau_with_vs_without_x2=0.0)
for name in PANELS:
    make_panel(name)
    df, units, times, Y, D, X = read_panel(name)
    T, N = Y.shape
    treated_units, k = block_info(D)
    assert (k is not None) == (PANELS[name][1] == "block")
    untreated = (D == 0) & np.isfinite(Y)
    ctrl_idx = np.where(~np.any(D == 1, axis=0))[0]

    # Residualize once: OLS on the untreated cells, no weights.
    ident_u = identified_phi(untreated, X)
    keep = np.where(ident_u)[0]
    _, phi_res = solve_wls(Y, untreated.astype(float), X[:, :, keep])
    phi_full = np.full(2, np.nan)
    phi_full[keep] = phi_res
    resid_rows.append(dict(panel=name, phi_1=phi_full[0], phi_2=phi_full[1]))

    for _, lam in lambdas.iterrows():
        cfg = lam["config"]
        lt, lu, ln = float(lam["time"]), float(lam["unit"]), float(lam["nn"])

        # Reference B, local method: one fit per treated cell.
        taus = []
        for t in range(T):
            for i in range(N):
                if D[t, i] != 1:
                    continue
                w = local_weights(Y, D, i, t, lt, lu) * untreated
                F, phi, val, F0 = solve_eq13(Y, w, X, ln)
                tau = Y[t, i] - F[t, i]
                taus.append(tau)
                checks["B_clarabel_only_vs_polished_tau"] = max(checks["B_clarabel_only_vs_polished_tau"],
                                                                abs(F0[t, i] - F[t, i]))
                ident = identified_phi(w > 0, X)
                cell_rows.append(dict(panel=name, config=cfg, method="local", unit=units[i],
                                      time=times[t], tau=tau, **phi_cols(phi, ident), objective=val))
                if not math.isfinite(ln):
                    Fw, phiw = solve_wls(Y, w, X)
                    checks["wls_vs_B_tau"] = max(checks["wls_vs_B_tau"], abs(Y[t, i] - Fw[t, i] - tau))
                    checks["wls_vs_B_phi"] = max(checks["wls_vs_B_phi"], float(np.max(np.abs(phiw - phi)[ident], initial=0)))
                if name == "XC3":
                    F1 = solve_eq13(Y, w, X[:, :, :1], ln)[0]
                    checks["XC3_tau_with_vs_without_x2"] = max(checks["XC3_tau_with_vs_without_x2"],
                                                               abs(Y[t, i] - F1[t, i] - tau))
        att_rows.append(dict(panel=name, config=cfg, method="local", reference="B",
                             att=float(np.mean(taus)), phi_1=np.nan, phi_2=np.nan))

        if k is None:
            continue
        # Reference B, global method.
        w = global_weights(Y, D, lt, lu, k) * np.isfinite(Y)
        F, phi, val, F0 = solve_eq13(Y, w, X, ln)
        ident = identified_phi(w > 0, X)
        tau_m = Y - F
        att_b = float(np.mean(tau_m[D == 1]))
        att_b0 = float(np.mean((Y - F0)[D == 1]))
        checks["B_clarabel_only_vs_polished_tau"] = max(checks["B_clarabel_only_vs_polished_tau"],
                                                        float(np.max(np.abs(F0 - F)[D == 1])))
        for t in range(T):
            for i in range(N):
                if D[t, i] == 1:
                    cell_rows.append(dict(panel=name, config=cfg, method="global", unit=units[i],
                                          time=times[t], tau=tau_m[t, i], **phi_cols(phi, ident),
                                          objective=val))
        att_rows.append(dict(panel=name, config=cfg, method="global", reference="B", att=att_b,
                             **phi_cols(phi, ident)))
        if not math.isfinite(ln):
            Fw, phiw = solve_wls(Y, w, X)
            checks["wls_vs_B_tau"] = max(checks["wls_vs_B_tau"], float(np.max(np.abs((Y - Fw - tau_m)[D == 1]))))
            checks["wls_vs_B_phi"] = max(checks["wls_vs_B_phi"], float(np.max(np.abs(phiw - phi)[ident], initial=0)))
        if name == "XC3":
            F1 = solve_eq13(Y, w, X[:, :, :1], ln)[0]
            checks["XC3_tau_with_vs_without_x2"] = max(checks["XC3_tau_with_vs_without_x2"],
                                                       float(np.max(np.abs((F1 - F)[D == 1]))))

        # Reference A (N x T layout of trop).
        att_a, phi_a = trop_twfe_average_cov(Y.T, D.T, X.transpose(1, 0, 2), treated_units, lu, lt, ln,
                                             treated_periods=k, solver=SOLVER, solver_opts=CLARABEL_OPTS)
        att_rows.append(dict(panel=name, config=cfg, method="global", reference="A", att=att_a,
                             **phi_cols(phi_a, ident)))
        d_att = abs(att_a - att_b)
        checks["A_vs_B_global_att"] = max(checks["A_vs_B_global_att"], d_att)
        key = "A_vs_B_global_att_finite_nn" if math.isfinite(ln) else "A_vs_B_global_att_inf_nn"
        checks[key] = max(checks[key], d_att)
        checks["A_vs_B_clarabel_only_global_att"] = max(checks["A_vs_B_clarabel_only_global_att"],
                                                        abs(att_a - att_b0))
        checks["A_vs_B_global_phi"] = max(checks["A_vs_B_global_phi"],
                                          float(np.max(np.abs(phi_a - phi)[ident], initial=0)))
        print(name, cfg, "done", flush=True)
    print(name, "done", flush=True)

fmt = "%.15g"
pd.DataFrame(att_rows).to_csv(f"{FIX}/oracle_cov.csv", index=False, float_format=fmt)
pd.DataFrame(cell_rows).to_csv(f"{FIX}/oracle_cov_cells.csv", index=False, float_format=fmt)
pd.DataFrame(resid_rows).to_csv(f"{FIX}/oracle_cov_resid.csv", index=False, float_format=fmt)
meta = {
    "trop_commit": commit(TROP_SRC),
    "diff_diff_commit": commit(DIFFDIFF_SRC),
    "solver": SOLVER,
    "solver_options": CLARABEL_OPTS,
    "tolerance": TOL,
    "python": platform.python_version(),
    "cvxpy": cp.__version__,
    "clarabel": clarabel.__version__,
    "numpy": np.__version__,
    "pandas": pd.__version__,
    "scipy": scipy.__version__,
    "panels": {k: dict(seed=v[0], design=v[1], treated_units=v[2], first_treated=v[3],
                       x2_constant_within_unit=v[4]) for k, v in PANELS.items()},
    "n_units": N_UNITS,
    "n_periods": N_PERIODS,
    "phi_true": PHI_TRUE,
    "tau_true": TAU_TRUE,
    "checks_max_abs_diff": checks,
    "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
}
with open(f"{FIX}/oracle_cov_meta.json", "w") as fh:
    json.dump(meta, fh, indent=2)
print(json.dumps(meta, indent=2))
