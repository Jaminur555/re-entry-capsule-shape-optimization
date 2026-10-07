"""Descriptor run: re-run every DOE point recording what the DOE dropped.

For each row of doe_results.csv (identical inputs), integrate with a relaxed
max_step (labels do not need fine sampling; families are separated by
t=1900s vs 4100s, far above coarse error) and record:
  captured (exact solver status), t_final, h_final, V_final,
  gamma_min, h_max_after_60km (skip indicator), sigma_engaged (bank law
  ever modulated past 5 deg), n_max (cross-check vs CSV).
Writes uq/doe_descriptors.csv (append + flush + resume, doe.py conventions).

Usage: python -m uq.doe_descriptors   (~3-4 h at dt_max=15; t_max=14400
per uq/diag_censor_tail.py: censored runs capture by 13.0 ks)
"""

import csv
import ctypes
import os
import time

import numpy as np
import pandas as pd

from capsule_opt import config
from capsule_opt.aerodynamics.aero_database import build_aero_database
from capsule_opt.properties import get_capsule_properties
from capsule_opt.trajectory import propagate_trajectory

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'data/doe_results.csv')
OUT = os.path.join(HERE, 'data/doe_descriptors.csv')
DT_MAX = 15.0   # relaxed step cap: labels only, families separated by >>1% in t
COLS = ['idx', 'captured', 't_final', 'h_final', 'V_final', 'gamma_min_deg',
        'h_max_after60', 'sigma_engaged']


def describe(idx, r):
    """evaluate_shape's plumbing (objectives.py) with relaxed dt_max."""
    shape_props = get_capsule_properties(r.rn, r.rs, r.r_theta)
    interp_CD, interp_CL, _, alpha_trim = build_aero_database(
        r.rn, r.rs, r.r_theta)

    def _scaled(interp, k):
        return lambda x: k * interp(x)

    interp_CD = _scaled(interp_CD, r.k_CD)
    interp_CL = _scaled(interp_CL, r.k_CL)
    m = config.CAPSULE_DENSITY * shape_props['volume']

    traj = propagate_trajectory(
        r.rn, r.rs, r.r_theta, m,
        gamma0_deg=r.gamma0_deg, Vo=r.V0,
        shape_props=shape_props, interp_CD=interp_CD, interp_CL=interp_CL,
        alpha_trim_arr=alpha_trim, dt_max=DT_MAX, t_max=14400.0)
    h = traj['h']
    # first time below 60 km, then max altitude afterwards (skip indicator)
    low = np.where(h < 60000.0)[0]
    h_post = h[low[0]:].max() if len(low) else h.max()
    return [idx,
            float('termination event' in traj['status']),
            float(traj['t'][-1]), float(h[-1]), float(traj['V'][-1]),
            float(traj['gamma_deg'].min()), float(h_post),
            float((traj['bank_deg'] > 5.0).any())]


def main():
    if os.name == 'nt':  # keep awake while running
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000001)
    df = pd.read_csv(SRC)
    done = 0
    if os.path.exists(OUT):
        done = sum(1 for _ in open(OUT)) - 1
    print(f"DESCRIPTORS: {len(df)} rows, resuming at {done}, out={OUT}",
          flush=True)
    t0 = time.time()
    with open(OUT, 'a', newline='') as f:
        w = csv.writer(f)
        if done == 0:
            w.writerow(COLS)
        for i in range(done, len(df)):
            w.writerow(describe(i, df.iloc[i]))
            f.flush()
            if (i + 1) % 25 == 0 or i + 1 == len(df):
                el = time.time() - t0
                eta = el / (i + 1 - done) * (len(df) - i - 1)
                print(f"  {i+1}/{len(df)}  elapsed {el/60:.1f} min  "
                      f"ETA {eta/60:.1f} min", flush=True)
    print(f"DONE. {len(df)} descriptor rows.", flush=True)


if __name__ == '__main__':
    main()
