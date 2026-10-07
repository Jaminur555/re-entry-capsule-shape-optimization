"""DOE runner: LHS over 7-dim (design + uncertainty) space -> doe_results.csv.

Usage: python -m uq.doe [N] [outfile]
Appends + flushes after every point; restart skips completed rows.
"""

import csv
import ctypes
import math
import os
import sys
import time

import numpy as np
from scipy.stats import qmc

from capsule_opt.optimization.objectives import evaluate_shape
from capsule_opt.uncertainty.variables import sample

N = int(sys.argv[1]) if len(sys.argv) > 1 else 500
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   sys.argv[2] if len(sys.argv) > 2 else 'data/doe_results.csv')

IN_COLS = ['rn', 'rs', 'r_theta', 'gamma0_deg', 'V0', 'k_CD', 'k_CL']
OUT_COLS = ['Qs', 'sg', 'q_stag_max', 'q_shldr_max', 'n_max', 'feasible',
            'g_q_stag', 'g_q_shldr', 'g_n_max', 'g_pitch_stable', 'g_trim_on_nose']


def build_inputs():
    d = qmc.LatinHypercube(d=3, seed=1).random(N)
    z = sample(N, seed=2)
    return np.column_stack([d, z])


def run_point(row):
    rn, rs, rt, g, V, kcd, kcl = row
    try:
        r = evaluate_shape(rn, rs, rt, aero_scales=(kcd, kcl),
                           entry_conditions={'gamma0_deg': g, 'Vo': V})
        det = r['constraints']['details']
        vals = [r['objectives']['Qs'], r['objectives']['sg'],
                det['q_stag_max'], det['q_shldr_max'], det['n_max'],
                float(r['feasible'])] + \
               [float(r['constraints'][k]) for k in
                ('g_q_stag', 'g_q_shldr', 'g_n_max', 'g_pitch_stable', 'g_trim_on_nose')]
        return [v if isinstance(v, float) else float(v) for v in vals]
    except Exception as e:  # record, keep going; count at end
        print(f"  point FAILED: {type(e).__name__}: {e}", flush=True)
        return [math.nan] * len(OUT_COLS)


def main():
    if os.name == 'nt':  # no idle sleep while the DOE runs
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000001)
    X = build_inputs()
    done = 0
    if os.path.exists(OUT):
        with open(OUT, newline='') as f:
            done = sum(1 for _ in csv.reader(f) if _) - 1
    new_file = done == 0
    print(f"DOE: N={N}, resuming at {done}, out={OUT}", flush=True)
    t0 = time.time()
    with open(OUT, 'a', newline='') as f:
        w = csv.writer(f)
        if new_file:
            w.writerow(IN_COLS + OUT_COLS)
        for i in range(done, N):
            w.writerow(list(X[i]) + run_point(X[i]))
            f.flush()
            if (i + 1) % 10 == 0 or i + 1 == N:
                el = time.time() - t0
                eta = el / (i + 1 - done) * (N - i - 1) if i + 1 > done else 0
                print(f"  {i+1}/{N}  elapsed {el/60:.1f} min  ETA {eta/60:.1f} min", flush=True)
    with open(OUT, newline='') as f:
        rows = list(csv.reader(f))[1:]
    bad = sum(1 for r in rows if r[7] in ('', 'nan'))  # Qs column
    print(f"DONE. {len(rows)} rows, {bad} failed/NaN ({100*bad/max(len(rows),1):.1f}%)", flush=True)


if __name__ == '__main__':
    main()
