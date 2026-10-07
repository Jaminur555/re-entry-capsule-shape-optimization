"""Diagnostic: how long do the censored (non-capture) runs need?

Re-runs ONLY the captured=False rows of diag_capture_sample.csv with a large
t_max (default 43200 s = 12 h) via the new entry_conditions['t_max'] hook.
Reports whether each eventually captures and at what time, so we can pick the
final t_max for the full DOE re-label.  No files modified.
Usage: python -m uq.diag_censor_tail
"""

import os
import time

import pandas as pd

from capsule_opt.optimization.objectives import evaluate_shape

HERE = os.path.dirname(os.path.abspath(__file__))
T_MAX = 43200.0


def run_row(r):
    r_ = evaluate_shape(
        r.rn, r.rs, r.r_theta, aero_scales=(r.k_CD, r.k_CL),
        entry_conditions={'gamma0_deg': r.gamma0, 'Vo': r.V0,
                          't_max': T_MAX})
    traj = r_['traj']
    return {'captured': 'termination event' in traj['status'],
            't_final': traj['t'][-1],
            'h_final_km': traj['h'][-1] / 1000.0,
            'V_final': traj['V'][-1]}


def main():
    doe = pd.read_csv(os.path.join(HERE, 'doe_results.csv'))
    smp = pd.read_csv(os.path.join(HERE, 'diag_capture_sample.csv'))
    cen = smp[~smp.captured].copy()
    cen = cen.join(doe[['rn', 'rs', 'r_theta', 'k_CD', 'k_CL', 'V0']],
                   on='idx')
    print(f"re-running {len(cen)} censored rows with t_max={T_MAX:.0f} s\n",
          flush=True)

    rows = []
    for i, (_, r) in enumerate(cen.iterrows(), 1):
        t0 = time.time()
        out = run_row(r)
        out.update(idx=r.idx, sg_km=r.sg_km)
        rows.append(out)
        print(f"  {i:2d}/{len(cen)}  sg={r.sg_km:8.0f} km  "
              f"captured={out['captured']}  t={out['t_final']:6.0f}s  "
              f"h_end={out['h_final_km']:6.1f} km  [{time.time()-t0:.0f}s]",
              flush=True)

    d = pd.DataFrame(rows)
    cap = d[d.captured]
    print(f"\neventually captured: {len(cap)}/{len(d)}")
    if len(cap):
        print(f"capture-time quantiles (s): "
              f"{cap.t_final.quantile([.5, .9, 1.0]).round(0).to_dict()}")
    d.to_csv(os.path.join(HERE, 'diag_censor_tail.csv'), index=False)
    print("\nsaved uq/diag_censor_tail.csv")


if __name__ == '__main__':
    main()
