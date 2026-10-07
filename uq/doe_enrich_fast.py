"""Fast-family enrichment DOE: +N rows in the direct-descent region.

sg-fast is the one substantial gate miss (0.906 R2 / 5.7%); the fast
family is smooth with a proven 1.0000 ceiling given t_final -- under-
resolved at n=898, not misspecified (build_surrogate.py v2 diagnostics).
This script samples from the same joint design+uncertainty machinery,
keeps candidates the current classifiers call fast and captured
(P(fast|x) > 0.98, P(capture|x) > 0.99), stages them in
uq/doe_fast1000.csv (append + flush + resume, doe.py conventions), and --
only for a full run (N >= 100) -- appends to doe_results.csv (backup
uq/doe_results_1800.csv). Fast rows are cheap (~10-15 s each).

Follow with: python -m uq.doe_descriptors  (resumes at row 1800)
             python -m uq.build_surrogate
Usage: python -m uq.doe_enrich_fast [N] [PFAST_MIN]
"""

import csv
import ctypes
import math
import os
import shutil
import sys
import time

import joblib
import numpy as np
from scipy.stats import qmc

from capsule_opt.optimization.objectives import evaluate_shape
from capsule_opt.uncertainty.variables import sample
from uq.doe import IN_COLS, OUT_COLS

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
PFAST = float(sys.argv[2]) if len(sys.argv) > 2 else 0.98
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, 'data/doe_results.csv')
STAGE = os.path.join(HERE, 'data/doe_fast1000.csv')
BACKUP = os.path.join(HERE, 'data/doe_results_1800.csv')
T_MAX = 14400.0
N_CAND = 8000


def run_point(row):
    """uq.doe.run_point + t_max=14400 passthrough (post-censoring-fix)."""
    rn, rs, rt, g, V, kcd, kcl = row
    try:
        r = evaluate_shape(rn, rs, rt, aero_scales=(kcd, kcl),
                           entry_conditions={'gamma0_deg': g, 'Vo': V,
                                             't_max': T_MAX})
        det = r['constraints']['details']
        vals = [r['objectives']['Qs'], r['objectives']['sg'],
                det['q_stag_max'], det['q_shldr_max'], det['n_max'],
                float(r['feasible'])] + \
               [float(r['constraints'][k]) for k in
                ('g_q_stag', 'g_q_shldr', 'g_n_max', 'g_pitch_stable',
                 'g_trim_on_nose')]
        return [v if isinstance(v, float) else float(v) for v in vals]
    except Exception as e:  # record, keep going; build drops NaN rows
        print(f"  point FAILED: {type(e).__name__}: {e}", flush=True)
        return [math.nan] * len(OUT_COLS)


def build_candidates():
    d = qmc.LatinHypercube(d=3, seed=21).random(N_CAND)
    z = sample(N_CAND, seed=22)
    X = np.column_stack([d, z])
    mix = joblib.load(os.path.join(HERE, 'artifacts/sg_mixture.joblib'))
    cap = joblib.load(os.path.join(HERE, 'artifacts/capture_clf.joblib'))
    Xs = mix['scaler'].transform(X)
    p_fast = 1.0 - mix['family_clf'].predict_proba(Xs)[:, 1]
    p_cap = cap['clf'].predict_proba(Xs)[:, 1]
    keep = np.where((p_fast > PFAST) & (p_cap > 0.99))[0]
    print(f"candidates: {N_CAND} -> {len(keep)} with P(fast)>{PFAST}, "
          f"P(cap)>0.99", flush=True)
    return X[keep[:N]]


def main():
    if os.name == 'nt':  # keep awake while running
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000001)
    with open(RESULTS, newline='') as f:
        n_rows = sum(1 for _ in csv.reader(f)) - 1
    if n_rows > 1800:
        print(f"results already enriched: {n_rows} rows; nothing to do.",
              flush=True)
        return
    cand = build_candidates()
    n_eff = min(N, len(cand))
    done = 0
    if os.path.exists(STAGE):
        with open(STAGE, newline='') as f:
            done = sum(1 for _ in csv.reader(f)) - 1
    print(f"ENRICH-FAST: N={n_eff} (P(fast)>{PFAST}), staged={done}, "
          f"~12 s/row ~ {n_eff * 12 / 3600:.1f} h", flush=True)
    t0 = time.time()
    newf = not os.path.exists(STAGE)
    with open(STAGE, 'a', newline='') as f:
        w = csv.writer(f)
        if newf:
            w.writerow(IN_COLS + OUT_COLS)
        for i in range(done, n_eff):
            w.writerow(list(cand[i]) + run_point(cand[i]))
            f.flush()
            if (i + 1) % 10 == 0 or i + 1 == n_eff:
                el = time.time() - t0
                eta = el / (i + 1 - done) * (n_eff - i - 1)
                print(f"  {i+1}/{n_eff}  elapsed {el/60:.1f} min  "
                      f"ETA {eta/60:.1f} min", flush=True)
    with open(STAGE, newline='') as f:
        staged = sum(1 for _ in csv.reader(f)) - 1
    if N < 100:  # smoke run: stage only, never touch the master csv
        print(f"SMOKE DONE: staged {staged}/{n_eff}; master csv untouched.",
              flush=True)
        return
    if staged < n_eff:
        print(f"INCOMPLETE: {staged}/{n_eff} staged; re-run to resume.",
              flush=True)
        return
    if not os.path.exists(BACKUP):
        shutil.copy2(RESULTS, BACKUP)
    with open(STAGE, newline='') as f:
        rows = list(csv.reader(f))
    with open(RESULTS, 'a', newline='') as f:
        csv.writer(f).writerows(rows[1:])
    with open(RESULTS, newline='') as f:
        n2 = sum(1 for _ in csv.reader(f)) - 1
    print(f"APPENDED: doe_results.csv now {n2} rows "
          f"(backup {os.path.basename(BACKUP)}).", flush=True)


if __name__ == '__main__':
    main()
