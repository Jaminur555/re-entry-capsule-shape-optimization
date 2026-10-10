"""Post-hoc fragility audit (Step 6): push a front CSV through the moment
engine at report precision (n_z=2000) and tabulate what the deterministic
view hides -- sigma_Qs, peak exceedance probabilities, P(slow), E_sg.

Defaults to the det front (the paper's "hidden cost of determinism" table);
works on any opt_*.csv with rn/rs/r_theta columns.

Usage: python -m uq.audit [uq/data/opt_det.csv]   (run AFTER opt runs;
~1-2 s/design at n_z=2000 -- do not run while an optimization is active.)
Writes <csv>_audit.csv and prints the summary table.
"""

import os
import sys
import time

import numpy as np
import pandas as pd

from uq.engine import MomentEngine

HERE = os.path.dirname(os.path.abspath(__file__))
N_Z = 2000


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        HERE, 'data/opt_det.csv')
    df = pd.read_csv(src)
    X = df[['rn', 'rs', 'r_theta']].to_numpy()

    t0 = time.time()
    m = MomentEngine(n_z=N_Z).moments(X)
    print(f"audit: {len(df)} designs @ n_z={N_Z} in {time.time() - t0:.0f} s",
          flush=True)
    for k, v in m.items():
        df[f'aud_{k}'] = v

    q = lambda c: m[c]
    print("\n===== det-front fragility audit (moments under z-dispersion) =====",
          flush=True)
    print(f"  sigma_Qs (infl): median {np.median(q('sd_Qs_infl'))/1e6:.1f} "
          f"MJ/m2  range {q('sd_Qs_infl').min()/1e6:.1f} - "
          f"{q('sd_Qs_infl').max()/1e6:.1f}  (CV median "
          f"{np.median(q('sd_Qs_infl')/q('E_Qs')):.1%})", flush=True)
    for o in ('q_stag_max', 'q_shldr_max', 'n_max'):
        pe = q(f'Pexc_{o}')
        print(f"  P({o} > limit): median {np.median(pe):.3f}  max {pe.max():.3f}"
              f"  frac > 1% {(pe > 0.01).mean():.1%}", flush=True)
    ps = q('P_slow')
    print(f"  P(slow): median {np.median(ps):.3f}  range {ps.min():.3f} - "
          f"{ps.max():.3f}", flush=True)
    rob_feas = np.ones(len(df), bool)
    from uq.optimize import LIMITS
    for o, L in LIMITS.items():
        rob_feas &= (q(f'E_{o}') + 2 * q(f'sd_{o}_infl')) <= L
    print(f"  robust-feasible fraction of det front (E+2sd_infl): "
          f"{rob_feas.mean():.1%}", flush=True)

    out = src.replace('.csv', '_audit.csv')
    df.to_csv(out, index=False)
    print(f"saved {out}", flush=True)


if __name__ == '__main__':
    main()
