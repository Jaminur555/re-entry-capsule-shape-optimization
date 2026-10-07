"""Diagnostic: three views of the surrogate gate to decide the redesign.

View A  log-space stats   -- is the failure an artifact of raw-space scoring?
View B  corridor-restricted -- do outputs pass once sg>CAP rows are excluded?
View C  tail classifier   -- is P(sg > CAP) predictable (skip front as prob.)?
Plus a censored-sg fit (target clipped at CAP) -- learnable smooth surface?

No artifacts saved; prints tables only.  ~10 min at n=1800.
Usage: python -m uq.diag_gate2 [CAP_km]
"""

import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, brier_score_loss,
                             roc_auc_score)
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import (PolynomialFeatures, StandardScaler)

from uq.build_surrogate import FOLDS, IN_COLS, SEED, fit_one

CAP_KM = float(sys.argv[1]) if len(sys.argv) > 1 else 8000.0
CAP = CAP_KM * 1000.0
HERE = os.path.dirname(os.path.abspath(__file__))


def reg_stats(y, pred, space):
    if space == 'log':
        y, pred = np.log(y), np.log(pred)
    rmse = float(np.sqrt(np.mean((pred - y) ** 2)))
    r2 = 1 - float(np.sum((pred - y) ** 2)) / float(np.sum((y - y.mean()) ** 2))
    return rmse, r2


def main():
    df = pd.read_csv(os.path.join(HERE, 'doe_results.csv')).dropna()
    X = StandardScaler().fit_transform(df[IN_COLS].to_numpy())
    tail = (df['sg'] > CAP).to_numpy()
    print(f"n={len(df)}  CAP={CAP_KM:.0f} km  tail rows={tail.sum()} "
          f"({100 * tail.mean():.1f}%)\n")

    kf = KFold(FOLDS, shuffle=True, random_state=SEED)

    # ---- View C: tail classifier (fast, do first) --------------------
    clf = Pipeline([
        ('poly', PolynomialFeatures(2, include_bias=False)),
        ('lr', LogisticRegression(max_iter=2000, C=1.0))])
    proba = cross_val_predict(clf, X, tail, cv=kf, method='predict_proba')[:, 1]
    print("C. tail classifier P(sg > CAP), logistic deg-2, 5-fold OOF:")
    print(f"   AUC={roc_auc_score(tail, proba):.4f}  "
          f"acc={accuracy_score(tail, proba > .5):.3f}  "
          f"brier={brier_score_loss(tail, proba):.4f}\n")

    # ---- Views A/B + censored fit -------------------------------------
    hdr = (f"{'output':<12}{'logR2':>8}{'logRMSE%':>10}"
           f"{'corR2':>8}{'corRMSE%':>10}   corridor = raw-space on sg<=CAP")
    print("A/B. regressors, 5-fold OOF:")
    print(hdr)
    for out in ['Qs', 'sg', 'q_stag_max', 'q_shldr_max', 'n_max']:
        t0 = time.time()
        y = df[out].to_numpy()
        oof = np.full(len(df), np.nan)
        for tr, te in kf.split(X):
            oof[te] = fit_one(X[tr], np.log(y[tr]), 0).predict(X[te])
        pred = np.exp(oof)
        r2_log = reg_stats(y, pred, 'log')[1]
        rel_log = 100 * float(np.sqrt(np.mean((np.log(pred) - np.log(y)) ** 2)))
        keep = ~tail
        rmse_c, r2_c = reg_stats(y[keep], pred[keep], 'raw')
        rng = df.loc[keep, out].max() - df.loc[keep, out].min()
        print(f"{out:<12}{r2_log:>8.4f}{rel_log:>10.2f}"
              f"{r2_c:>8.4f}{100 * rmse_c / rng:>10.2f}"
              f"   [{time.time() - t0:.0f}s]", flush=True)

    # censored sg: learn min(sg, CAP); scored vs censored target
    t0 = time.time()
    yc = np.minimum(df['sg'].to_numpy(), CAP)
    oof = np.full(len(df), np.nan)
    for tr, te in kf.split(X):
        oof[te] = fit_one(X[tr], np.log(yc[tr]), 0).predict(X[te])
    pred = np.minimum(np.exp(oof), CAP)
    rmse, r2 = reg_stats(yc, pred, 'raw')
    r2_log = reg_stats(yc, pred, 'log')[1]
    print(f"\ncensored sg (min(sg,CAP)) vs censored target: "
          f"R2={r2:.4f}  RMSE/range={100 * rmse / (CAP - yc.min()):.2f}%  "
          f"logR2={r2_log:.4f}  [{time.time() - t0:.0f}s]")


if __name__ == '__main__':
    main()
