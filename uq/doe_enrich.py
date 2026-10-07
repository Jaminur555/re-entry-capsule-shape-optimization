"""Targeted enrichment DOE: loft-corner box where OOF residuals concentrate.

diag_residuals.py (09-22) showed worst GP errors at shallow gamma0 x small
r_theta — the loft/skip corner. Adds N points (default 300) with r_theta
restricted to [0, RT_HI], gamma0 truncated to [G_LO, G_HI] of its N(-2, 0.5)
law; V0/k_CD/k_CL from the frozen sampler, unchanged. Appends to
doe_results.csv with doe.py's flush/resume conventions.
Usage: python -m uq.doe_enrich [N]
"""

import csv
import ctypes
import os
import sys
import time

import numpy as np
from scipy.stats import norm, qmc

from capsule_opt.uncertainty.variables import sample
from uq.doe import OUT_COLS, run_point

N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'data/doe_results.csv')
IN_COLS = ['rn', 'rs', 'r_theta', 'gamma0_deg', 'V0', 'k_CD', 'k_CL']

MU, SD = -2.0, 0.5        # gamma0 law (uncertainty/variables.py)
G_LO, G_HI = -2.4, -0.25  # shallow-loft enrichment window [deg]
RT_HI = 0.45              # r_theta restricted to [0, RT_HI]


def build_inputs():
    d = qmc.LatinHypercube(d=3, seed=11).random(N)
    d[:, 2] *= RT_HI
    z = sample(N, seed=12)  # (gamma0, V0, k_CD, k_CL); gamma0 overwritten below
    lo, hi = norm.cdf([(G_LO - MU) / SD, (G_HI - MU) / SD])
    u = qmc.LatinHypercube(d=1, seed=13).random(N).ravel()
    g = MU + SD * norm.ppf(lo + u * (hi - lo))
    return np.column_stack([d, g, z[:, 1:]])


def main():
    if os.name == 'nt':  # no idle sleep while the DOE runs
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000001)
    X = build_inputs()
    with open(OUT, newline='') as f:
        done = sum(1 for _ in csv.reader(f)) - 1
    print(f"ENRICH: N={N} (r_theta<= {RT_HI}, gamma0 in [{G_LO},{G_HI}]), "
          f"existing rows={done}, out={OUT}", flush=True)
    t0 = time.time()
    with open(OUT, 'a', newline='') as f:
        w = csv.writer(f)
        for i in range(N):
            w.writerow(list(X[i]) + run_point(X[i]))
            f.flush()
            if (i + 1) % 10 == 0 or i + 1 == N:
                el = time.time() - t0
                eta = el / (i + 1) * (N - i - 1)
                print(f"  {i+1}/{N}  elapsed {el/60:.1f} min  ETA {eta/60:.1f} min",
                      flush=True)
    with open(OUT, newline='') as f:
        rows = list(csv.reader(f))[1:]
    print(f"DONE. {len(rows)} total rows now.", flush=True)


if __name__ == '__main__':
    main()
