"""Step-6 moment engine: surrogate statistics over a fixed z-LHS (CRN).

Loads the v2 artifacts (build_surrogate.py) and, for a batch of designs
(rn, rs, r_theta), predicts over ONE fixed uncertainty LHS reused for every
candidate (common random numbers -> smooth NSGA-II fitness):

  E_z[f], std_z[f] (plain + GP-variance-inflated via the law of total
  variance), P(peak > limit), P(slow family), P(capture).

The two z-free geometry constraints (g_pitch_stable, g_trim_on_nose) are
surrogated from the DOE columns (exact evaluation would need an aero
database build per candidate): a trim-existence classifier carries pitch
stability (in the DOE, pitch-stable <=> a trim exists) and a within-trim
HGBR carries the nose margin. eta_V stays exact (pure geometry).

Self-test (conventions guard): python -m uq.engine  -- point-predicts a few
DOE rows and compares against the recorded CSV values.
"""

import os

import joblib
import numpy as np
from scipy.stats import norm as _norm

from capsule_opt import config
from capsule_opt.optimization.metrics import volumetric_efficiency
from capsule_opt.uncertainty.variables import sample as _sample_z

HERE = os.path.dirname(os.path.abspath(__file__))
ART = os.path.join(HERE, 'artifacts')
IN7 = ['rn', 'rs', 'r_theta', 'gamma0_deg', 'V0', 'k_CD', 'k_CL']
IN3 = ['rn', 'rs', 'r_theta']
GP_OUTPUTS = ['Qs', 'q_stag_max', 'q_shldr_max', 'n_max']
LIMITS = {'q_stag_max': config.C_Q_STAG_MAX,
          'q_shldr_max': config.C_Q_SHLDR_MAX,
          'n_max': config.C_N_MAX}
Z0 = np.array([-2.0, 7830.0, 1.0, 1.0])   # nominal (gamma0, V0, k_CD, k_CL)
N_Z, Z_SEED = 2000, 42                    # PLAN.md Step 6: fixed CRN
CHUNK = 10000                             # GP predict rows per block


# ---------------------------------------------------------------- geometry
def build_geom_surrogate():
    """Geometry constraints from the DOE columns (z-free in (rn, rs, r_theta)).

    DOE facts: pitch-stability <=> a trim exists (all 2327 trimmed rows are
    stable, g ~ -0.003; the 473 no-trim rows carry the sentinel 1.0/90.0) ->
    classification problem. The nose margin is continuous WITHIN the trimmed
    region only -> regressor on trimmed rows. Gates: pooled 5-fold OOF
    accuracy / sign-accuracy. Cached to artifacts/geom_margins.joblib.
    """
    import pandas as pd
    from sklearn.ensemble import (GradientBoostingClassifier,
                                  HistGradientBoostingRegressor)
    from sklearn.model_selection import KFold, cross_val_predict

    df = pd.read_csv(os.path.join(HERE, 'data/doe_results.csv')).dropna(
        subset=IN3 + ['g_pitch_stable', 'g_trim_on_nose'])
    X = df[IN3].to_numpy()
    has_trim = (df['g_pitch_stable'] <= 0).to_numpy().astype(int)

    p = cross_val_predict(GradientBoostingClassifier(random_state=0), X,
                          has_trim, cv=KFold(5, shuffle=True,
                                             random_state=0),
                          method='predict_proba')[:, 1]
    acc_trim = float(((p > 0.5) == has_trim).mean())
    print(f"geom trim-clf      OOF acc {acc_trim:.4f} "
          f"(base {max(has_trim.mean(), 1 - has_trim.mean()):.4f})",
          flush=True)

    tr = has_trim == 1
    y = df.loc[tr, 'g_trim_on_nose'].to_numpy()
    oof = cross_val_predict(HistGradientBoostingRegressor(), X[tr], y,
                            cv=KFold(5, shuffle=True, random_state=0))
    acc_nose = float(((oof <= 0) == (y <= 0)).mean())
    print(f"geom nose-within   OOF sign-acc {acc_nose:.4f} "
          f"(base {max((y <= 0).mean(), (y > 0).mean()):.4f})", flush=True)

    out = {'in_cols': IN3,
           'trim_clf': GradientBoostingClassifier(random_state=0).fit(
               X, has_trim),
           'nose': HistGradientBoostingRegressor().fit(X[tr], y),
           'acc_trim': acc_trim, 'acc_nose': acc_nose}
    joblib.dump(out, os.path.join(ART, 'geom_margins.joblib'))
    return out


def load_geom():
    p = os.path.join(ART, 'geom_margins.joblib')
    return joblib.load(p) if os.path.exists(p) else build_geom_surrogate()


# ------------------------------------------------------------------ engine
def _log_moments(mu, sd):
    """Per-sample lognormal (mu, sd) -> unit-space mean & variance."""
    m = np.exp(mu + 0.5 * sd ** 2)
    v = (np.exp(sd ** 2) - 1.0) * np.exp(2.0 * mu + sd ** 2)
    return m, v


class MomentEngine:
    """All surrogate predictions for Step 6. One instance per optimization."""

    def __init__(self, n_z=N_Z, seed=Z_SEED):
        self.n_z = n_z
        self.Z = _sample_z(n_z, seed)                    # (n_z, 4)
        self.gp = {o: joblib.load(os.path.join(ART, f'gp_{o}.joblib'))
                   for o in GP_OUTPUTS}
        self.mix = joblib.load(os.path.join(ART, 'sg_mixture.joblib'))
        self.cap = joblib.load(os.path.join(ART, 'capture_clf.joblib'))
        self.geom = load_geom()

    # -- low level ------------------------------------------------------
    def _geom(self, Xd):
        """Effective geometry margins at (n,3) designs.

        g_pitch = 0.5 - P(trim)  (<=0 feasible; pitch-stability <=> trim
        exists in the DOE). g_nose = within-trim HGBR margin, sentinel 90
        where the classifier blocks the no-trim region.
        """
        p_trim = self.geom['trim_clf'].predict_proba(Xd)[:, 1]
        g_nose = self.geom['nose'].predict(Xd)
        return 0.5 - p_trim, np.where(p_trim > 0.5, g_nose, 90.0)

    def _expand(self, Xd):
        """(n,3) designs + fixed z-LHS -> (n*n_z, 7) raw input rows."""
        n = len(Xd)
        X7 = np.empty((n * self.n_z, 7))
        X7[:, :3] = np.repeat(Xd, self.n_z, axis=0)
        X7[:, 3:] = np.tile(self.Z, (n, 1))
        return X7

    def _predict_log(self, art, X7, chunk=CHUNK):
        """Scaled GP predict with std, chunked; returns (mu, sd) log-space."""
        Xs = art['scaler'].transform(X7)
        mu = np.empty(len(Xs))
        sd = np.empty(len(Xs))
        for a in range(0, len(Xs), chunk):
            b = a + chunk
            m, s = art['gp'].predict(Xs[a:b], return_std=True)
            mu[a:b], sd[a:b] = m, s
        return mu, sd

    # -- deterministic baseline ------------------------------------------
    def det_predict(self, Xd):
        """Point predictions at nominal z0 (thesis formulation on surrogate)."""
        X7 = np.empty((len(Xd), 7))
        X7[:, :3] = Xd
        X7[:, 3:] = Z0
        out = {'eta_V': np.array([volumetric_efficiency(*x) for x in Xd])}
        out['g_pitch_stable'], out['g_trim_on_nose'] = self._geom(Xd)
        for o in GP_OUTPUTS:
            mu, _ = self._predict_log(self.gp[o], X7)
            out[o] = np.exp(mu)
        out.update(self._mixture_at(X7))
        out['P_slow'] = out.pop('p_slow')
        out['P_cap'] = out.pop('p_cap')
        return out

    def _mixture_at(self, X7):
        """sg mixture + family/capture probabilities at arbitrary input rows."""
        Xs = self.mix['scaler'].transform(X7)
        p_slow = self.mix['family_clf'].predict_proba(Xs)[:, 1]
        p_cap = self.cap['clf'].predict_proba(Xs)[:, 1]
        mu_f, sd_f = self._predict_log(
            {'gp': self.mix['gp_fast'], 'scaler': self.mix['scaler']}, X7)
        mu_s, sd_s = self._predict_log(
            {'gp': self.mix['gp_slow'], 'scaler': self.mix['scaler']}, X7)
        m_f, _ = _log_moments(mu_f, sd_f)
        m_s, _ = _log_moments(mu_s, sd_s)
        return {'sg': (1.0 - p_slow) * m_f + p_slow * m_s,
                'p_slow': p_slow, 'p_cap': p_cap}

    # -- robust moments ---------------------------------------------------
    def moments(self, Xd):
        """CRN-LHS moments for a batch of designs (n,3) -> dict of (n,) arrays.

        std_plain  = std over z of the point predictions exp(mu);
        std_infl   = sqrt( mean_z var_pred + var_z mean_pred )  (law of
        total variance: includes the GP's own predictive uncertainty).
        """
        n, nz = len(Xd), self.n_z
        X7 = self._expand(Xd)
        shp = lambda a: a.reshape(n, nz)
        out = {'eta_V': np.array([volumetric_efficiency(*x) for x in Xd])}
        out['g_pitch_stable'], out['g_trim_on_nose'] = self._geom(Xd)

        mix = self._mixture_at(X7)
        out['E_sg'] = shp(mix['sg']).mean(axis=1)
        out['sd_sg_plain'] = shp(mix['sg']).std(axis=1, ddof=0)
        out['P_slow'] = shp(mix['p_slow']).mean(axis=1)
        out['P_cap'] = shp(mix['p_cap']).mean(axis=1)

        for o in GP_OUTPUTS:
            mu, sd = self._predict_log(self.gp[o], X7)
            m, v = _log_moments(mu, sd)
            pm, pv = shp(m), shp(v)
            mean = pm.mean(axis=1)
            var_plain = shp(np.exp(mu)).var(axis=1, ddof=0)
            out[f'E_{o}'] = mean
            out[f'sd_{o}_plain'] = np.sqrt(var_plain)
            out[f'sd_{o}_infl'] = np.sqrt(pv.mean(axis=1) +
                                          pm.var(axis=1, ddof=0))
            if o in LIMITS:  # P(peak > limit), lognormal tail per sample
                lnL = np.log(LIMITS[o])
                out[f'Pexc_{o}'] = shp(_norm.sf((lnL - mu) / sd)).mean(axis=1)
        return out


# ------------------------------------------------------------- self-test
def selftest():
    """Conventions guard: point-predict DOE rows, compare to recorded values."""
    import pandas as pd
    df = pd.read_csv(os.path.join(HERE, 'data/doe_results.csv')).dropna(
        subset=GP_OUTPUTS + ['sg'])
    idx = [0, 777, 1500, 2799]
    eng = MomentEngine()
    print(f"selftest on doe_results rows {idx} (rel diff vs CSV; "
          f"expect ~ gate RMSE, i.e. a few %)", flush=True)
    for i in idx:
        row = df.iloc[i]
        Xd = np.array([[row[c] for c in IN3]])
        X7 = np.array([[row[c] for c in IN7]])
        mix = eng._mixture_at(X7)
        pred = {}
        for o in GP_OUTPUTS:
            mu, _ = eng._predict_log(eng.gp[o], X7)
            pred[o] = float(np.exp(mu[0]))
        pred['sg'] = float(mix['sg'][0])
        print(f"  row {i}: " + "  ".join(
            f"{k} {abs(pred[k]-row[k])/abs(row[k]):.3f}"
            for k in GP_OUTPUTS + ['sg']), flush=True)
    m = eng.moments(Xd := np.array([[0.5, 0.5, 0.5]]))
    print("  moments(D0) shapes ok:", all(len(v) == 1 for v in m.values()),
          f"| E_Qs {m['E_Qs'][0]:.3e}  sd_infl {m['sd_Qs_infl'][0]:.3e}"
          f"  E_sg {m['E_sg'][0]/1e3:.0f} km  P_slow {m['P_slow'][0]:.3f}",
          flush=True)


if __name__ == '__main__':
    selftest()
