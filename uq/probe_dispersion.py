"""Probe-0 for Step 6: true uncertainty-driven dispersion at fixed geometry.

Option A (moments) makes std[Qs] a fitness objective, but the surrogate's
own gate error (~2.3% of range for Qs) may be the same size as the true
sigma_z. This probe measures the real thing: 5 fixed designs x 60
uncertainty samples (z ~ uncertainty model, t_max=14400), then reports
per design the sample mean/std of each output against the final-gate
surrogate RMSE. Verdict rule: sigma/rmse > ~2 -> std[Qs] objective is
resolvable as-is; ~1 -> use GP-variance-inflated std in the optimizer.

Designs (fixed, span the space incl. the loft/slow corner):
  D0 (0.50, 0.50, 0.50)  center, known feasible (PLAN probe point)
  D1 (0.35, 0.50, 0.50)  low rn  (heating-critical)
  D2 (0.65, 0.50, 0.50)  high rn (blunt)
  D3 (0.50, 0.40, 0.65)  geometry spread
  D4 (0.45, 0.45, 0.25)  small r_theta -> loft/slow-prone (slow sims ~50 s)

Usage (from capsule_opt/): python -m uq.probe_dispersion [N_SAMP] [N_DESIGN]
  N_SAMP   samples per design   (default 60; smoke = 2)
  N_DESIGN first N designs only (default 5;   smoke = 1)
Resume-safe: staged rows keyed on (design_id, isamp) in
uq/data/probe_dispersion.csv; smoke rows are a prefix of the full z-LHS
(sample(60, seed=101) sliced), so they are reused, not wasted.
"""

import csv
import ctypes
import math
import os
import sys
import time

import numpy as np

from capsule_opt.optimization.objectives import evaluate_shape
from capsule_opt.uncertainty.variables import sample
from uq.doe import IN_COLS, OUT_COLS

N_SAMP = int(sys.argv[1]) if len(sys.argv) > 1 else 60
N_DESIGN = int(sys.argv[2]) if len(sys.argv) > 2 else 5
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'data/probe_dispersion.csv')
RESULTS = os.path.join(HERE, 'data/doe_results.csv')
T_MAX = 14400.0
SLOW_CUT = 2250.0          # build v2 diagnostics: log-t gap 2241-2297 s
N_FULL, SEED = 60, 101     # full z-LHS, sliced -> smoke rows stay valid

DESIGNS = [(0.50, 0.50, 0.50), (0.35, 0.50, 0.50), (0.65, 0.50, 0.50),
           (0.50, 0.40, 0.65), (0.45, 0.45, 0.25)]

# final gate (n=2800, CV_RESTARTS=2, 2026-10-08): RMSE as % of DOE range
GATE_RMSE_PCT = {'Qs': 2.332, 'q_stag_max': 2.294, 'q_shldr_max': 3.009,
                 'n_max': 2.225}
COLS = ['design_id', 'isamp'] + IN_COLS + OUT_COLS + ['t_final', 'captured']


def run_point(design, z):
    rn, rs, rt = design
    g, V, kcd, kcl = z
    try:
        r = evaluate_shape(rn, rs, rt, aero_scales=(kcd, kcl),
                           entry_conditions={'gamma0_deg': g, 'Vo': V,
                                             't_max': T_MAX})
        det = r['constraints']['details']
        tr = r['traj']
        vals = [r['objectives']['Qs'], r['objectives']['sg'],
                det['q_stag_max'], det['q_shldr_max'], det['n_max'],
                float(r['feasible'])] + \
               [float(r['constraints'][k]) for k in
                ('g_q_stag', 'g_q_shldr', 'g_n_max', 'g_pitch_stable',
                 'g_trim_on_nose')]
        vals += [float(tr['t'][-1]),
                 float('termination event' in tr['status'])]
        return [v if isinstance(v, float) else float(v) for v in vals]
    except Exception as e:
        print(f"  point FAILED: {type(e).__name__}: {e}", flush=True)
        return [math.nan] * (len(OUT_COLS) + 2)


def summarize():
    import pandas as pd
    df = pd.read_csv(OUT)
    try:
        rng = pd.read_csv(RESULTS)
        ranges = {k: float(rng[k].max() - rng[k].min())
                  for k in GATE_RMSE_PCT}
    except Exception:
        ranges = None
    print("\n===== Probe-0 summary (captured rows only) =====", flush=True)
    for did, d in enumerate(DESIGNS[:N_DESIGN]):
        sub = df[(df['design_id'] == did) & (df['captured'] == 1)]
        n_all = int((df['design_id'] == did).sum())
        n_slow = int((sub['t_final'] > SLOW_CUT).sum())
        p_feas = df[df['design_id'] == did]['feasible'].mean()
        print(f"\nD{did} rn/rs/rt = {d}   n={n_all} "
              f"(captured {len(sub)}, slow {n_slow}, "
              f"P(feasible) {p_feas:.2f})", flush=True)
        print("  output        mean        sigma     CV %    gate RMSE   "
              "sigma/RMSE", flush=True)
        for k, pct in GATE_RMSE_PCT.items():
            v = sub[k].dropna()
            if len(v) < 2:
                continue
            mu, sd = float(v.mean()), float(v.std(ddof=1))
            rmse = pct / 100.0 * ranges[k] if ranges else math.nan
            print(f"  {k:<12} {mu:11.4g} {sd:10.4g} {100*sd/mu:6.1f}   "
                  f"{rmse:9.4g}   {sd/rmse:8.2f}" if ranges else
                  f"  {k:<12} {mu:11.4g} {sd:10.4g} {100*sd/mu:6.1f}",
                  flush=True)
    print("\nVerdict rule: sigma/RMSE > ~2 -> plain std[Qs] objective OK; "
          "~1 -> inflate with GP predictive variance.", flush=True)


def main():
    if os.name == 'nt':  # keep awake while running
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000001)
    Z = sample(N_FULL, SEED)[:N_SAMP]          # slice of the fixed LHS
    done = set()
    if os.path.exists(OUT):
        with open(OUT, newline='') as f:
            for r in csv.DictReader(f):
                done.add((int(r['design_id']), int(r['isamp'])))
    todo = [(d, i) for d in range(N_DESIGN) for i in range(N_SAMP)
            if (d, i) not in done]
    print(f"PROBE-0: {N_DESIGN} designs x {N_SAMP} samples, {len(todo)} to "
          f"run (~15-50 s each), t_max={T_MAX:.0f} s", flush=True)
    t0 = time.time()
    newf = not os.path.exists(OUT)
    with open(OUT, 'a', newline='') as f:
        w = csv.writer(f)
        if newf:
            w.writerow(COLS)
        for n, (d, i) in enumerate(todo, 1):
            rn, rs, rt = DESIGNS[d]
            w.writerow([d, i, rn, rs, rt] + list(Z[i]) +
                       run_point(DESIGNS[d], Z[i]))
            f.flush()
            if n % 10 == 0 or n == len(todo):
                el = time.time() - t0
                eta = el / n * (len(todo) - n)
                print(f"  {n}/{len(todo)}  elapsed {el/60:.1f} min  "
                      f"ETA {eta/60:.1f} min", flush=True)
    summarize()


if __name__ == '__main__':
    main()
