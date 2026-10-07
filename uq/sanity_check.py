"""Sanity checks for the GP surrogates (PLAN.md Step 5).

1. Corner check: surrogate vs direct evaluate_shape at nominal z
   (gamma0=-2 deg, V0=7830 m/s, k=1) over fresh random design points;
   deviations compared against the Step-4 RMSE/range bands.
2. Envelope check: DOE sg range vs thesis Table 6.2 (2930-4695 km).
Usage: python -m uq.sanity_check
"""

import os

import joblib
import numpy as np
import pandas as pd

from capsule_opt.optimization.objectives import evaluate_shape

HERE = os.path.dirname(os.path.abspath(__file__))
BANDS = {'Qs': 0.02, 'sg': 0.02, 'q_stag_max': 0.03, 'q_shldr_max': 0.03,
         'n_max': 0.03}
N_PTS = 10


def main():
    arts = {o: joblib.load(os.path.join(HERE, 'artifacts', f'gp_{o}.joblib'))
            for o in BANDS}

    # --- 1. corner check at nominal z -------------------------------------
    d = np.random.RandomState(7).uniform(size=(N_PTS, 3))
    z_nom = [-2.0, 7830.0, 1.0, 1.0]
    dev = {o: [] for o in BANDS}
    for i, (rn, rs, rt) in enumerate(d):
        r = evaluate_shape(rn, rs, rt)
        det = r['constraints']['details']
        direct = {'Qs': r['objectives']['Qs'], 'sg': r['objectives']['sg'],
                  'q_stag_max': det['q_stag_max'],
                  'q_shldr_max': det['q_shldr_max'], 'n_max': det['n_max']}
        x = np.array([[rn, rs, rt] + z_nom])
        for o, art in arts.items():
            xs = art['scaler'].transform(x)
            p = art['gp'].predict(xs)[0]
            if art['log']:
                p = np.exp(p)
            dev[o].append(abs(p - direct[o]) / np.diff(art['y_range'])[0])
        print(f"  point {i + 1}/{N_PTS}", flush=True)

    print(f"\ncorner check (|surrogate-direct| / range, {N_PTS} pts):")
    print(f"{'output':<12}{'median':>9}{'max':>9}{'band':>7}  verdict")
    for o, band in BANDS.items():
        a = np.array(dev[o])
        ok = a.mean() < band  # RMSE-level criterion
        print(f"{o:<12}{np.median(a):>9.3%}{a.max():>9.3%}{band:>7.0%}"
              f"  {'PASS' if ok else 'FAIL'}")

    # --- 2. envelope check vs thesis Table 6.2 -----------------------------
    df = pd.read_csv(os.path.join(HERE, 'data/doe_results.csv'))
    feas = df[df['feasible'] > 0]
    t62 = (2930.0, 4695.0)
    lo, hi = feas['sg'].min() / 1e3, feas['sg'].max() / 1e3
    med = feas['sg'].median() / 1e3
    print(f"\nenvelope: feasible-DOE sg [{lo:.0f}, {hi:.0f}] km "
          f"(median {med:.0f}) vs thesis Table 6.2 {t62} km")
    print("contains thesis range" if lo <= t62[0] and hi >= t62[1]
          else "DOES NOT contain thesis range")


if __name__ == '__main__':
    main()
