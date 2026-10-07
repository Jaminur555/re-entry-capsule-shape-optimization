"""Decisive test: capture-labelled, branch-conditional surrogate viability.

1. Re-run the 8 ambiguous rows (n_max in the [0.3, 1.5) g moat) for their true
   solver capture status; all other rows labelled by n_max > 0.5 g (verified
   exact on the 50-pt diag_capture sample: capped <= 0.12 g, captured >= 1.8).
2. Capture classifier: logistic deg-2, 5-fold OOF AUC / acc / Brier.
3. Regressors on CAPTURED rows only: 5 outputs x 5-fold OOF, raw + log stats
   vs the v1 gate bands.
Usage: python -m uq.diag_gate3
"""

import os
import time

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

from capsule_opt.optimization.objectives import evaluate_shape
from uq.build_surrogate import FOLDS, IN_COLS, SEED, fit_one

HERE = os.path.dirname(os.path.abspath(__file__))
G_LO, G_HI = 0.3, 1.5   # the empirical moat between the populations [g]
BANDS = {'Qs': 0.02, 'sg': 0.02, 'q_stag_max': 0.03, 'q_shldr_max': 0.03,
         'n_max': 0.03}


def true_status(r):
    r_ = evaluate_shape(r.rn, r.rs, r.r_theta, aero_scales=(r.k_CD, r.k_CL),
                        entry_conditions={'gamma0_deg': r.gamma0_deg,
                                          'Vo': r.V0})
    return 'termination event' in r_['traj']['status']


def main():
    df = pd.read_csv(os.path.join(HERE, 'doe_results.csv'))
    moat = df.n_max.between(G_LO, G_HI)
    print(f"n={len(df)}  moat rows to relabel exactly: {moat.sum()}")
    label = (df.n_max > 0.5).to_numpy()
    for i in np.where(moat)[0]:
        t0 = time.time()
        lab = true_status(df.iloc[i])
        label[i] = lab
        print(f"  row {i}: n_max={df.n_max.iloc[i]:.2f} g -> captured={lab} "
              f"[{time.time()-t0:.0f}s]", flush=True)
    print(f"captured: {label.sum()} ({100*label.mean():.1f}%)  "
          f"capped: {(~label).sum()}\n")

    X = StandardScaler().fit_transform(df[IN_COLS].to_numpy())
    kf = KFold(FOLDS, shuffle=True, random_state=SEED)

    clf = Pipeline([('poly', PolynomialFeatures(2, include_bias=False)),
                    ('lr', LogisticRegression(max_iter=2000))])
    p = cross_val_predict(clf, X, label, cv=kf, method='predict_proba')[:, 1]
    print("capture classifier (logistic deg-2, 5-fold OOF):")
    print(f"   AUC={roc_auc_score(label, p):.4f}  "
          f"acc={accuracy_score(label, p > .5):.4f}  "
          f"brier={brier_score_loss(label, p):.4f}\n")

    sub = df[label]
    Xs = StandardScaler().fit_transform(sub[IN_COLS].to_numpy())
    print(f"regressors on captured rows (n={len(sub)}), 5-fold OOF:")
    print(f"{'output':<12}{'R2_raw':>8}{'RMSE/rng':>10}{'logR2':>8}"
          f"{'logRMSE%':>10}   gate (R2>0.95, band)")
    for out, band in BANDS.items():
        t0 = time.time()
        y = sub[out].to_numpy()
        oof = np.full(len(sub), np.nan)
        for tr, te in KFold(FOLDS, shuffle=True, random_state=SEED).split(Xs):
            oof[te] = fit_one(Xs[tr], np.log(y[tr]), 0).predict(Xs[te])
        pred = np.exp(oof)
        r2 = 1 - np.sum((pred - y) ** 2) / np.sum((y - y.mean()) ** 2)
        frac = np.sqrt(np.mean((pred - y) ** 2)) / (y.max() - y.min())
        lr2 = 1 - np.sum((np.log(pred) - np.log(y)) ** 2) / \
            np.sum((np.log(y) - np.log(y).mean()) ** 2)
        lrm = 100 * np.sqrt(np.mean((np.log(pred) - np.log(y)) ** 2))
        ok = 'PASS' if (r2 > 0.95 and frac < band) else 'fail'
        print(f"{out:<12}{r2:>8.4f}{frac:>10.2%}{lr2:>8.4f}{lrm:>10.2f}"
              f"   {ok} ({band:.0%})  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == '__main__':
    main()
