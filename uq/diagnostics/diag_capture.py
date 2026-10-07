"""Diagnostic: measure the t_max-censoring (non-capture) fraction of the DOE.

Root-cause hypothesis: propagate_trajectory stops at a terminal Mach-3 event
with t_max=7200 s; runs that never reach Mach 3 are censored and their 'sg'
is endpoint-arc-at-7200s, a different response branch.  The DOE CSV never
recorded the status.  This re-runs a stratified sample (5 per sg-decile)
through evaluate_shape with IDENTICAL arguments and reports the true capture
label, t_final, final state.  No files modified.
Usage: python -m uq.diag_capture
"""

import os
import time

import numpy as np
import pandas as pd

from capsule_opt.optimization.objectives import evaluate_shape

HERE = os.path.dirname(os.path.abspath(__file__))
PER_DECILE = 5


def run_row(r):
    r_ = evaluate_shape(
        r.rn, r.rs, r.r_theta, aero_scales=(r.k_CD, r.k_CL),
        entry_conditions={'gamma0_deg': r.gamma0_deg, 'Vo': r.V0})
    traj = r_['traj']
    return {
        'captured': 'termination event' in traj['status'],
        't_final': traj['t'][-1],
        'h_final': traj['h'][-1],
        'V_final': traj['V'][-1],
        'sg_re': traj['sg'],
    }


def main():
    df = pd.read_csv(os.path.join(HERE, 'doe_results.csv'))
    df['decile'] = pd.qcut(df['sg'] / 1000.0, 10, labels=False)
    sample = df.groupby('decile', group_keys=False).apply(
        lambda g: g.sample(PER_DECILE, random_state=7))
    print(f"re-running {len(sample)} of {len(df)} rows "
          f"({PER_DECILE} per sg-decile)\n", flush=True)

    rows = []
    for i, (idx, r) in enumerate(sample.iterrows(), 1):
        t0 = time.time()
        out = run_row(r)
        out.update(idx=idx, sg_km=r.sg / 1000.0, decile=r.decile,
                   gamma0=r.gamma0_deg, rt=r.r_theta)
        rows.append(out)
        print(f"  {i:2d}/{len(sample)}  sg={out['sg_km']:8.0f} km  "
              f"captured={out['captured']}  t={out['t_final']:6.0f}s  "
              f"h_end={out['h_final']/1000:6.1f} km  "
              f"[{time.time()-t0:.0f}s]", flush=True)

    d = pd.DataFrame(rows)
    assert np.allclose(d.sg_km * 1000, df.loc[d.idx, 'sg'].to_numpy(), rtol=1e-9), \
        "re-run sg mismatch -- non-determinism!"
    print(f"\ndeterminism check vs CSV: PASS (max rel diff "
          f"{np.max(np.abs(d.sg_km*1000 - df.loc[d.idx,'sg'].to_numpy())/df.loc[d.idx,'sg'].to_numpy()):.1e})")
    print(f"\ncaptured: {d.captured.sum()}/{len(d)} "
          f"({100*d.captured.mean():.0f}%)")
    print("\nby sg-decile (0 = shortest):")
    print(d.groupby('decile').captured.agg(['mean', 'count']).round(2).to_string())
    print("\ncaptured rows: t_final med %.0f s   capped rows: t_final %s" % (
        d.loc[d.captured, 't_final'].median(),
        sorted(d.loc[~d.captured, 't_final'].round(0).unique())))
    d.to_csv(os.path.join(HERE, 'diag_capture_sample.csv'), index=False)
    print("\nsaved uq/diag_capture_sample.csv")


if __name__ == '__main__':
    main()
