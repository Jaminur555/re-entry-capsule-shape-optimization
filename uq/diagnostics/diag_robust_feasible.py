"""Robust-feasibility scan (Step-6 pre-flight): 300 LHS designs through the
moment engine at n_z=500 -- how big is the feasible region under each
constraint form, and which constraints bind?

Prints feasibility fractions under: det point constraints; E+2*sd_plain;
E+2*sd_infl; P_slow <= 0.05; geometry; and the joint robust set -- plus
per-constraint active fractions and g ranges.

Usage: python -m uq.diagnostics.diag_robust_feasible [N=300] [N_Z=500]
"""

import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import qmc

from capsule_opt import config
from uq.engine import MomentEngine

N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
N_Z = int(sys.argv[2]) if len(sys.argv) > 2 else 500
HERE = os.path.dirname(os.path.abspath(__file__))
LIMITS = {'q_stag_max': config.C_Q_STAG_MAX,
          'q_shldr_max': config.C_Q_SHLDR_MAX,
          'n_max': config.C_N_MAX}


def gform(m, o, kind):
    L = LIMITS[o]
    sd = m[f'sd_{o}_{kind}']
    return (m[f'E_{o}'] + 2.0 * sd - L) / L


def main():
    X = qmc.LatinHypercube(d=3, seed=33).random(N)
    eng = MomentEngine(n_z=N_Z)
    t0 = time.time()
    m = eng.moments(X)
    print(f"{N} designs x n_z={N_Z}: {time.time() - t0:.0f} s "
          f"({(time.time() - t0) / N:.2f} s/cand)", flush=True)

    print("\n--- feasibility by constraint form (fraction g <= 0) ---",
          flush=True)
    rows = []
    for o in LIMITS:
        rows.append((f'det {o}', float((m[f'E_{o}'] <= LIMITS[o]).mean())))
        rows.append((f'plain {o}', float((gform(m, o, "plain") <= 0).mean())))
        rows.append((f'infl   {o}', float((gform(m, o, "infl") <= 0).mean())))
    p_slow_ok = m['P_slow'] <= 0.05
    rows.append(('P_slow<=0.05', float(p_slow_ok.mean())))
    geom = (m['g_pitch_stable'] <= 0) & (m['g_trim_on_nose'] <= 0)
    rows.append(('geometry', float(geom.mean())))
    for k, v in rows:
        print(f"  {k:<18} {v:.1%}", flush=True)

    joint_plain = geom & p_slow_ok & \
        np.all([gform(m, o, 'plain') <= 0 for o in LIMITS], axis=0)
    joint_infl = geom & p_slow_ok & \
        np.all([gform(m, o, 'infl') <= 0 for o in LIMITS], axis=0)
    joint_nopslow_infl = geom & np.all(
        [gform(m, o, 'infl') <= 0 for o in LIMITS], axis=0)
    print(f"\nJOINT robust (plain):   {joint_plain.mean():.1%}", flush=True)
    print(f"JOINT robust (infl):    {joint_infl.mean():.1%}", flush=True)
    print(f"JOINT infl w/o P_slow:  {joint_nopslow_infl.mean():.1%}",
          flush=True)

    print("\n--- binding info on the joint-infl feasible set ---", flush=True)
    for o in LIMITS:
        g = gform(m, o, 'infl')[joint_infl]
        if len(g):
            print(f"  {o:<12} min {g.min():+.3f}  p50 {np.median(g):+.3f}  "
                  f"max {g.max():+.3f}  active(>-0.01) "
                  f"{(g > -0.01).mean():.0%}", flush=True)
    ps = m['P_slow'][joint_infl]
    if len(ps):
        print(f"  P_slow       min {ps.min():.3f}  p50 {np.median(ps):.3f}  "
              f"max {ps.max():.3f}", flush=True)
    print(f"\nE_Qs on feasible: {m['E_Qs'][joint_infl].min() / 1e6:.1f} - "
          f"{m['E_Qs'][joint_infl].max() / 1e6:.1f} MJ/m2" if joint_infl.any()
          else "  (empty set)", flush=True)
    pd.DataFrame({**{f'x{i + 1}': X[:, i] for i in range(3)}, **m}).to_csv(
        os.path.join(HERE, '..', 'data', 'diag_robust_feasible.csv'),
        index=False)


if __name__ == '__main__':
    main()
