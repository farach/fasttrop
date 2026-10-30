"""Reference results for the oracle tests of fasttrop (condition C1 of the council).

The script calculates results with fixed tuning parameters on the fixture
panels. It uses two independent implementations at the commits in
oracle_meta.json:

* diff-diff (igerber/diff-diff, MIT): the local method (Eq. 12 and
  Algorithm 2 of the paper) and the global method for a block design.
* trop (zhaonanq/TROP-Estimator, MIT, the package of the authors of the
  paper): the global method for a block design (TROP_TWFE_average), with
  CVXPY and CLARABEL.

The script writes these files to tests/testthat/fixtures/:
* oracle.csv            The ATT for each panel, configuration, method, and
                        reference.
* oracle_cells.csv      The cell effects tau_it of diff-diff (local and
                        global).
* oracle_solver.csv     The weights of one weighted fit and the fitted values
                        of diff-diff, alpha_i + beta_t + L_it (+ mu for the
                        global method). The file has these values only for
                        the configurations with a finite nuclear norm
                        penalty. Thus, the tests can compare the solvers with
                        the same weights.
* oracle_objective.csv  The objective value of diff-diff for these fits.

The cell effects of diff-diff come from a copy of the loop in
TROPLocalMixin._fit_with_fixed_lambda, with the internal functions of
diff-diff. The script makes sure that their mean is equal to the ATT of that
method.

Usage, in the package folder, with clones of the two repositories:
    python data-raw/oracle_fixtures.py <path-to-TROP-Estimator> <path-to-diff-diff>
"""
import json
import math
import subprocess
import sys
import time

import numpy as np
import pandas as pd

TROP_SRC, DIFFDIFF_SRC = sys.argv[1], sys.argv[2]
sys.path.insert(0, TROP_SRC + "/src")
sys.path.insert(0, DIFFDIFF_SRC)

import trop  # noqa: E402
from diff_diff import TROP  # noqa: E402
from diff_diff.trop_local import _treated_cell_is_estimable  # noqa: E402

FIX = "tests/testthat/fixtures"
TOL = 1e-9
MAX_ITER = 200000


def commit(path):
    return subprocess.run(["git", "-C", path, "rev-parse", "HEAD"], capture_output=True,
                          text=True, check=True).stdout.strip()


def read_panel(name):
    df = pd.read_csv(f"{FIX}/{name}.csv")
    units = sorted(df["unit"].unique())
    times = sorted(df["time"].unique())
    # The layout of diff-diff: periods in rows and units in columns.
    Y = df.pivot(index="time", columns="unit", values="y").reindex(index=times, columns=units).values
    D = df.pivot(index="time", columns="unit", values="d").reindex(index=times, columns=units)
    D = D.fillna(0).astype(int).values
    return df, units, times, Y, D


def nuclear(L):
    return float(np.linalg.svd(L, compute_uv=False).sum())


def objective(Y, w, F, L, lam_nn):
    r = np.where(w > 0, Y - F, 0.0)
    val = float(np.sum(w * r * r))
    return val + (lam_nn * nuclear(L) if math.isfinite(lam_nn) else 0.0)


lambdas = pd.read_csv(f"{FIX}/lambdas.csv")
lambdas["nn"] = lambdas["nn"].replace("Inf", np.inf).astype(float)

att_rows, cell_rows, solver_rows, obj_rows = [], [], [], []
for name in ["D1", "D2", "D3", "D4"]:
    df, units, times, Y, D = read_panel(name)
    T, N = Y.shape
    treated_units = [i for i in range(N) if D[:, i].sum() > 0]
    first = [int(np.argmax(D[:, i] == 1)) for i in treated_units]
    block = len(set(first)) == 1 and all(D[first[0]:, i].all() for i in treated_units)
    k = T - first[0] if block else None
    control_mask = D == 0
    control_unit_idx = np.where(~np.any(D == 1, axis=0))[0]
    for _, lam in lambdas.iterrows():
        cfg = lam["config"]
        lt, lu, ln = float(lam["time"]), float(lam["unit"]), float(lam["nn"])

        # Local method: the loop of diff-diff for each cell (_fit_with_fixed_lambda). Keep each tau.
        est = TROP(method="local", max_iter=MAX_ITER, tol=TOL)
        t0 = time.time()
        taus, first_fit = [], None
        for t in range(T):
            for i in range(N):
                if D[t, i] != 1 or not np.isfinite(Y[t, i]):
                    continue
                wm = est._compute_observation_weights(Y, D, i, t, lt, lu, control_unit_idx, N, T)
                if not _treated_cell_is_estimable(control_mask, Y, wm, i, t):
                    continue
                alpha, beta, L = est._estimate_model(Y, control_mask, wm, ln, N, T)
                tau = Y[t, i] - alpha[i] - beta[t] - L[t, i]
                taus.append(tau)
                cell_rows.append(dict(dataset=name, config=cfg, method="local", reference="diff-diff",
                                      unit=units[i], time=times[t], tau=tau))
                if first_fit is None:
                    w_eff = wm * (control_mask & np.isfinite(Y))
                    F = alpha[np.newaxis, :] + beta[:, np.newaxis] + L
                    first_fit = (i, t, w_eff, F, L)
        att_local = float(np.mean(taus))
        seconds = time.time() - t0
        if cfg == "twfe":
            ref = est._fit_with_fixed_lambda(df, "y", "d", "unit", "time", (lt, lu, ln))
            assert abs(ref - att_local) < 1e-12, (name, ref, att_local)
        att_rows.append(dict(dataset=name, config=cfg, method="local", reference="diff-diff",
                             att=att_local, seconds=seconds))
        if math.isfinite(ln):
            i0, t0_, w_eff, F, L = first_fit
            for t in range(T):
                for i in range(N):
                    solver_rows.append(dict(dataset=name, config=cfg, method="local", unit=units[i],
                                            time=times[t], w=w_eff[t, i], fitted=F[t, i]))
            obj_rows.append(dict(dataset=name, config=cfg, method="local", target_unit=units[i0],
                                 target_time=times[t0_], lambda_nn=ln,
                                 objective=objective(Y, w_eff, F, L, ln)))

        if block:
            # Global method: the package of the authors (only the ATT).
            t0 = time.time()
            att_trop = trop.TROP_TWFE_average(Y.T, D.T, treated_units, lu, lt, ln,
                                              treated_periods=k, solver="CLARABEL")
            att_rows.append(dict(dataset=name, config=cfg, method="global", reference="trop",
                                 att=att_trop, seconds=time.time() - t0))
            # Global method: diff-diff (cell effects, weights, and fitted values).
            g = TROP(method="global", max_iter=MAX_ITER, tol=TOL)
            t0 = time.time()
            delta = g._compute_global_weights(Y, D, lt, lu, k, N, T)
            mu, alpha, beta, L = g._solve_global_model(Y, delta, ln)
            att_g, effects, _ = g._extract_posthoc_tau(Y, D, mu, alpha, beta, L)
            att_rows.append(dict(dataset=name, config=cfg, method="global", reference="diff-diff",
                                 att=att_g, seconds=time.time() - t0))
            for (i, t), tau in effects.items():
                cell_rows.append(dict(dataset=name, config=cfg, method="global", reference="diff-diff",
                                      unit=units[i], time=times[t], tau=tau))
            if math.isfinite(ln):
                w_eff = np.where(np.isfinite(Y), delta, 0.0)
                F = mu + alpha[np.newaxis, :] + beta[:, np.newaxis] + L
                for t in range(T):
                    for i in range(N):
                        solver_rows.append(dict(dataset=name, config=cfg, method="global",
                                                unit=units[i], time=times[t], w=w_eff[t, i],
                                                fitted=F[t, i]))
                obj_rows.append(dict(dataset=name, config=cfg, method="global", target_unit=np.nan,
                                     target_time=np.nan, lambda_nn=ln,
                                     objective=objective(Y, w_eff, F, L, ln)))
        print(name, cfg, "done", flush=True)

fmt = "%.15g"
pd.DataFrame(att_rows).to_csv(f"{FIX}/oracle.csv", index=False, float_format=fmt)
pd.DataFrame(cell_rows).to_csv(f"{FIX}/oracle_cells.csv", index=False, float_format=fmt)
pd.DataFrame(solver_rows).to_csv(f"{FIX}/oracle_solver.csv", index=False, float_format=fmt)
pd.DataFrame(obj_rows).to_csv(f"{FIX}/oracle_objective.csv", index=False, float_format=fmt)
meta = {
    "trop_commit": commit(TROP_SRC),
    "diff_diff_commit": commit(DIFFDIFF_SRC),
    "diff_diff_tol": TOL,
    "diff_diff_max_iter": MAX_ITER,
    "global_solver_trop": "CLARABEL",
    "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
}
with open(f"{FIX}/oracle_meta.json", "w") as fh:
    json.dump(meta, fh, indent=2)
print(json.dumps(meta, indent=2))
