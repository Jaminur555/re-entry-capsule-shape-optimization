"""Fit surrogate v2 on the cleaned doe_results.csv (PLAN.md Step 4).

v1 failed its gate because the 7200 s horizon had censored 202 rows into a
second response branch (uq/diag_*.py; fixed by uq/doe_rerun202.py). This
build works on the cleaned dataset plus uq/doe_descriptors.csv labels:

  - capture classifier over all 1800 rows (5 true non-captures);
  - trajectory families from descriptors t_final: bimodal in log t with an
    empty gap 2241-2297 s -> SLOW_CUT = 2250 s (fast = direct descent,
    302-2241 s; slow = lofted equilibrium glide, 2297-13474 s). Linear
    logistic family classifier on the 7 inputs (OOF AUC ~0.999);
  - sg: one GP per family in log space. The fast family is smooth (log sg
    ~ f(log t_final), corr 0.997). The slow family is NOT modelable at
    practical DOE density (GP/HGBR/RF <= 0.66 log R2, no input
    correlation, ~11% GP noise floor: sensitivity amplifies along the
    ~1e4 s bank-modulated glide) -> gated IRREDUCIBLE and carried as an
    explicit dispersion term into the robust formulation;
  - Qs / q_stag_max / q_shldr_max / n_max: single captured-branch GPs in
    log space (gate: pooled 5-fold OOF, R2 > 0.95, RMSE/range band).

Artifacts: uq/gp_<out>.joblib (captured branch), uq/sg_mixture.joblib
(family classifier + both sg GPs), uq/capture_clf.joblib.
Usage: python -m uq.build_surrogate [out ...]   (subset of the four
non-sg outputs; family classifier + sg mixture are always built)
"""

import os
import sys
import time

import joblib
import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
IN_COLS = ['rn', 'rs', 'r_theta', 'gamma0_deg', 'V0', 'k_CD', 'k_CL']
OUTPUTS = {  # output column -> RMSE/range band (R2 band is 0.95 for all)
    'Qs': 0.02,
    'q_stag_max': 0.03, 'q_shldr_max': 0.03, 'n_max': 0.03,
}
SG_BAND = 0.02     # sg gate band applies to the fast family
SLOW_CUT = 2250.0  # s; families are separated by an empty 2241-2297 s gap
FOLDS = 5
CV_RESTARTS = 2     # fold fits need real optima: r0 gate is conservative
FINAL_RESTARTS = 3  # artifact fit on all rows
SEED = 0


def fit_one(Xtr, ytr, restarts=10):
    kernel = ConstantKernel(1.0, (1e-3, 1e3)) * Matern(
        length_scale=np.ones(Xtr.shape[1]), length_scale_bounds=(1e-2, 1e2),
        nu=2.5) + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-6, 1e1))
    gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True,
                                  n_restarts_optimizer=restarts,
                                  random_state=0)
    gp.fit(Xtr, ytr)
    return gp


def gate(y, pred):
    rmse = float(np.sqrt(np.mean((pred - y) ** 2)))
    return rmse, rmse / (y.max() - y.min()), r2_score(y, pred)


def ard_str(gp):
    ls = gp.kernel_.k1.k2.length_scale
    return '  '.join(f"{c}={v:.2f}" for c, v in zip(IN_COLS, ls))


def main():
    only = sys.argv[1:] or list(OUTPUTS)
    df = pd.read_csv(os.path.join(HERE, 'data/doe_results.csv'))
    desc = pd.read_csv(os.path.join(HERE, 'data/doe_descriptors.csv'))
    assert len(df) == len(desc), 'results/descriptors row mismatch'
    captured = (desc['captured'] == 1).to_numpy()
    slow = captured & (desc['t_final'] >= SLOW_CUT).to_numpy()
    fast = captured & ~slow

    scaler = StandardScaler().fit(df[IN_COLS].to_numpy())
    X = scaler.transform(df[IN_COLS].to_numpy())

    print(f"{len(df)} rows | captured {captured.sum()} "
          f"(fast {fast.sum()} / slow {slow.sum()}) | {FOLDS}-fold pooled OOF "
          f"+ final fit (restarts={FINAL_RESTARTS})", flush=True)

    # ---- capture classifier (all rows) --------------------------------
    t0 = time.time()
    p_cap = cross_val_predict(LogisticRegression(max_iter=2000), X,
                              captured.astype(int),
                              cv=KFold(FOLDS, shuffle=True,
                                       random_state=SEED),
                              method='predict_proba')[:, 1]
    auc_cap = roc_auc_score(captured, p_cap)
    cap_clf = LogisticRegression(max_iter=2000).fit(X, captured.astype(int))
    joblib.dump({'clf': cap_clf, 'scaler': scaler, 'in_cols': IN_COLS,
                 'oof_auc': auc_cap, 'n_noncapture': int((~captured).sum())},
                os.path.join(HERE, 'artifacts/capture_clf.joblib'))
    print(f"capture_clf   OOF AUC {auc_cap:.4f} "
          f"({(~captured).sum()} non-captures)  "
          f"{time.time() - t0:.0f} s  saved", flush=True)

    # ---- family classifier + sg mixture -------------------------------
    t0 = time.time()
    fam = slow[captured].astype(int)
    p_fam = cross_val_predict(LogisticRegression(max_iter=2000),
                              X[captured], fam,
                              cv=KFold(FOLDS, shuffle=True,
                                       random_state=SEED),
                              method='predict_proba')[:, 1]
    auc_fam = roc_auc_score(fam, p_fam)
    print(f"family_clf    OOF AUC {auc_fam:.4f}  acc "
          f"{((p_fam > 0.5) == fam).mean():.3f}  "
          f"(slow cut t_final >= {SLOW_CUT:.0f} s)", flush=True)
    fam_clf = LogisticRegression(max_iter=2000).fit(X[captured], fam)

    y_sg = df['sg'].to_numpy()
    bag_log = {}   # fold-bagged log-sg prediction at every captured row
    gps_sg = {}
    for name, mask in [('sg-fast', fast), ('sg-slow', slow)]:
        t1 = time.time()
        y = y_sg[mask]
        yt = np.log(y)
        Xf = X[mask]
        oof = np.full(mask.sum(), np.nan)
        bagged = np.zeros(captured.sum())
        for tr, te in KFold(FOLDS, shuffle=True,
                            random_state=SEED).split(Xf):
            m = fit_one(Xf[tr], yt[tr], CV_RESTARTS)
            oof[te] = m.predict(Xf[te])
            bagged += m.predict(X[captured])   # for the mixture reference
        bagged /= FOLDS
        bag_log[name] = bagged
        pred = np.exp(oof)
        rmse, frac, r2 = gate(y, pred)
        if name == 'sg-fast':
            ok = r2 > 0.95 and frac < SG_BAND
            verdict = 'PASS' if ok else 'FAIL'
        else:
            verdict = 'IRREDUCIBLE (dispersion term)'
        print(f"{name:<12} pooled OOF  rmse {rmse:>10.4g}  {frac:>7.3%}"
              f"  R2 {r2:>7.4f}  {time.time() - t1:>4.0f} s  {verdict}",
              flush=True)
        gp = fit_one(Xf, yt, FINAL_RESTARTS)
        print(f"  ARD ls (full fit): {ard_str(gp)}  "
              f"noise={gp.kernel_.k2.noise_level:.2g}", flush=True)
        gps_sg[name] = gp

    # pooled mixture reference (OOF probs x fold-bagged family predictions)
    mix_pred = (1 - p_fam) * np.exp(bag_log['sg-fast']) + \
        p_fam * np.exp(bag_log['sg-slow'])
    rmse, frac, r2 = gate(y_sg[captured], mix_pred)
    print(f"{'sg-mixture':<12} (reference) rmse {rmse:>10.4g}  {frac:>7.3%}"
          f"  R2 {r2:>7.4f}  (pooled, OOF-weighted)", flush=True)
    joblib.dump({'family_clf': fam_clf, 'gp_fast': gps_sg['sg-fast'],
                 'gp_slow': gps_sg['sg-slow'], 'scaler': scaler,
                 'in_cols': IN_COLS, 'log': True, 'slow_cut': SLOW_CUT,
                 'fam_oof_auc': auc_fam,
                 'y_range': (y_sg[captured].min(), y_sg[captured].max())},
                os.path.join(HERE, 'artifacts/sg_mixture.joblib'))
    stale = os.path.join(HERE, 'artifacts/gp_sg.joblib')
    if os.path.exists(stale):
        os.remove(stale)
        print("removed stale v1 artifact gp_sg.joblib", flush=True)
    print(f"sg block {time.time() - t0:.0f} s  saved sg_mixture.joblib",
          flush=True)

    # ---- captured-branch single GPs ------------------------------------
    print(f"{'output':<12}{'RMSE_cv':>12}{'RMSE/range':>11}{'R2_cv':>9}"
          f"{'t [s]':>7}  verdict", flush=True)
    for out in only:
        band = OUTPUTS[out]
        t0 = time.time()
        y = df.loc[captured, out].to_numpy()
        yt = np.log(y)
        oof = np.full(captured.sum(), np.nan)
        for tr, te in KFold(FOLDS, shuffle=True,
                            random_state=SEED).split(X[captured]):
            oof[te] = fit_one(X[captured][tr], yt[tr], CV_RESTARTS).predict(
                X[captured][te])
        pred = np.exp(oof)
        rmse, frac, r2 = gate(y, pred)
        ok = r2 > 0.95 and frac < band
        print(f"{out:<12}{rmse:>12.4g}{frac:>11.3%}{r2:>9.4f}"
              f"{time.time() - t0:>7.0f}  {'PASS' if ok else 'FAIL'}",
              flush=True)
        gp = fit_one(X[captured], yt, FINAL_RESTARTS)
        print(f"  ARD ls (full fit): {ard_str(gp)}  "
              f"noise={gp.kernel_.k2.noise_level:.2g}", flush=True)
        joblib.dump({'gp': gp, 'scaler': scaler, 'in_cols': IN_COLS,
                     'log': True, 'branch': 'captured',
                     'y_range': (y.min(), y.max())},
                    os.path.join(HERE, f'artifacts/gp_{out}.joblib'))
    print("saved uq/gp_<out>.joblib + sg_mixture.joblib + capture_clf.joblib",
          flush=True)


if __name__ == '__main__':
    main()
