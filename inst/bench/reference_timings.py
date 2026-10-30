"""Matched timings for the Python TROP references (council condition C4, D-A-O1).

Times the same fixed-lambda computations that inst/bench/reference_timings.R
times for fasttrop and cfcompare, on the same panels, with each package's
default solver settings:

* diff-diff (pure Python, no Rust extension): local TROP ATT, global ATT, and
  one leave-one-out evaluation Q(lambda) over every untreated cell.
* trop (the paper authors' package): global ATT via CVXPY.

Usage (from the package root, after reference_timings.R has written the panels):
    python inst/bench/reference_timings.py <path-to-TROP-Estimator> <path-to-diff-diff>
"""
import sys
import time

import numpy as np
import pandas as pd

TROP_SRC, DIFFDIFF_SRC = sys.argv[1], sys.argv[2]
sys.path.insert(0, TROP_SRC + "/src")
sys.path.insert(0, DIFFDIFF_SRC)

import trop  # noqa: E402
from diff_diff import TROP  # noqa: E402

BENCH = "inst/bench"
LAM = (0.3, 0.8, 0.5)  # (time, unit, nn)
REPS = 3


def best_of(fn, reps=REPS):
    best, val = np.inf, None
    for _ in range(reps):
        t0 = time.perf_counter()
        val = fn()
        best = min(best, time.perf_counter() - t0)
    return best, val


rows = []
for panel in ["panel_20x15", "panel_40x20"]:
    df = pd.read_csv(f"{BENCH}/{panel}.csv")
    units = sorted(df["unit"].unique())
    times = sorted(df["time"].unique())
    Y = df.pivot(index="time", columns="unit", values="y").reindex(index=times, columns=units).values
    D = df.pivot(index="time", columns="unit", values="d").reindex(index=times, columns=units)
    D = D.fillna(0).astype(int).values
    T, N = Y.shape
    treated = [i for i in range(N) if D[:, i].sum() > 0]
    k = int(D[:, treated[0]].sum())
    lt, lu, ln = LAM

    est = TROP(method="local")
    s, v = best_of(lambda: est._fit_with_fixed_lambda(df, "y", "d", "unit", "time", LAM))
    rows.append(dict(panel=panel, task="local ATT, fixed lambda", implementation="diff-diff (Python)",
                     seconds=s, value=v))

    g = TROP(method="global")
    s, v = best_of(lambda: g._fit_global_with_fixed_lambda(df, "y", "d", "unit", "time", LAM, k))
    rows.append(dict(panel=panel, task="global ATT, fixed lambda", implementation="diff-diff (Python)",
                     seconds=s, value=v))

    s, v = best_of(lambda: trop.TROP_TWFE_average(Y.T, D.T, treated, lu, lt, ln, treated_periods=k))
    rows.append(dict(panel=panel, task="global ATT, fixed lambda",
                     implementation="trop (authors, Python/CVXPY)", seconds=s, value=v))

    if panel == "panel_20x15":
        control_mask = D == 0
        control_unit_idx = np.where(~np.any(D == 1, axis=0))[0]
        s, v = best_of(lambda: est._loocv_score_obs_specific(Y, D, control_mask, control_unit_idx,
                                                             lt, lu, ln, N, T), reps=1)
        rows.append(dict(panel=panel, task="LOOCV Q(lambda), all untreated cells",
                         implementation="diff-diff (Python)", seconds=s, value=v))
    print(panel, "done", flush=True)

pd.DataFrame(rows).to_csv(f"{BENCH}/reference-timings-python.csv", index=False, float_format="%.10g")
print(pd.DataFrame(rows).to_string())
