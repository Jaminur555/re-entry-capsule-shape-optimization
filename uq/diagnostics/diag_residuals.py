"""One-off diagnosis (surrogate redesign): where do OOF residuals live in X?

5-fold OOF per output (same kernel/CV settings as build_surrogate, restarts=0),
then compares the worst-20% |residual| points vs the rest, input by input.
Localized errors -> targeted enrichment / regime split.
Diffuse errors -> trend (PCK-style) or kernel redesign.
Usage: python -m uq.diag_residuals
"""
import os

import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
IN = ['rn', 'rs', 'r_theta', 'gamma0_deg', 'V0', 'k_CD', 'k_CL']


def oof_predict(X, y):
    oof = np.empty_like(y)
    for tr, te in KFold(5, shuffle=True, random_state=0).split(X):
        k = (ConstantKernel(1.0, (1e-3, 1e3))
             * Matern(np.ones(X.shape[1]), (1e-2, 1e2), nu=2.5)
             + WhiteKernel(1e-2, (1e-6, 1e1)))
        gp = GaussianProcessRegressor(kernel=k, normalize_y=True,
                                      random_state=0).fit(X[tr], y[tr])
        oof[te] = gp.predict(X[te])
    return oof


def main():
    df = pd.read_csv(os.path.join(HERE, 'doe_results.csv')).dropna()
    sc = StandardScaler().fit(df[IN].to_numpy())
    X = sc.transform(df[IN].to_numpy())
    for out in ['sg', 'Qs', 'n_max']:
        y = np.log(df[out].to_numpy())
        r = np.abs(oof_predict(X, y) - y)
        cut = np.quantile(r, 0.8)
        bad, ok = r >= cut, r < cut
        print(f"\n=== {out}: |log-resid| worst20% median {np.median(r[bad]):.4f} "
              f"(~{np.expm1(np.median(r[bad])):.1%} rel) vs rest {np.median(r[ok]):.4f}")
        print(f"{'input':<12}{'worst20% med':>14}{'rest med':>12}{'sep':>7}")
        for c in IN:
            bw, ro = np.median(df.loc[bad, c]), np.median(df.loc[ok, c])
            spread = np.percentile(df[c], 90) - np.percentile(df[c], 10)
            print(f"{c:<12}{bw:>14.4g}{ro:>12.4g}{(bw-ro)/spread:>7.2f}")


if __name__ == '__main__':
    main()
