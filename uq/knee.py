"""Step 6: knee-point designs + strict non-dominated fronts (paper table).

Dominance is re-checked at report precision (robust objectives from the
n_z=2000 re-evaluation in wide_frame, not the loop's n_z=500 F columns),
the strict non-dominated fronts are written to uq/data/opt_<mode>_nd.csv,
and each front's knee is picked in min-max-normalized objective space as
the point with maximum perpendicular distance from the hyperplane through
the M anchor points (Branke et al. 2004, generalized to M>2); pseudo-weight
deviation ||w - 1/M|| is the cross-check / fallback when the anchors are
affinely dependent (robust front: one design anchors both E[Qs] and sd[Qs]).

Nominal (z0) table values come from the SIMULATOR (evaluate_shape,
t_max=14400, ~15 s/design, cached in uq/data/knee_sim_z0.csv), not from the
surrogate: the det front's surrogate sg(z0) inflates via slow-branch GP
extrapolation at loft-prone shapes (knee overpredicted 3.5x; see PLAN.md).

Prints the knee table (physical units, robust lens for both designs) and
saves uq/data/knee_designs.csv.

Usage: python -m uq.knee     (~40 s fresh incl. 2 simulator runs; needs
uq.optimize + uq.audit done)
"""

import os
import time

import numpy as np
import pandas as pd

from capsule_opt.geometry.parameterization import cap_params
from uq.engine import LIMITS, MomentEngine

HERE = os.path.dirname(os.path.abspath(__file__))
SIM_CACHE = os.path.join(HERE, 'data/knee_sim_z0.csv')
OBJ = {'det': ('eta_V', 'Qs', 'sg'),                  # minimized as
       'robust': ('eta_V', 'E_Qs', 'E_sg', 'sd_Qs_infl')}  # (-, +, -, +)
SIGN = np.array([-1.0, 1.0, -1.0, 1.0])


def load_front(mode):
    src = 'data/opt_det_audit.csv' if mode == 'det' else 'data/opt_robust.csv'
    df = pd.read_csv(os.path.join(HERE, src))
    if mode == 'det':   # audit columns carry the robust view at n_z=2000
        df = df.rename(columns={c: c[4:] for c in df.columns
                                if c.startswith('aud_')})
        df = df.loc[:, ~df.columns.duplicated()]
    return df


def F_min(df, mode):
    return df[list(OBJ[mode])].to_numpy(float) * SIGN[:len(OBJ[mode])]


def strict_nd(F):
    """True where no other point weakly dominates with one strict objective."""
    keep = np.ones(len(F), bool)
    for i in range(len(F)):
        keep[i] = not ((F <= F[i]).all(axis=1) &
                       (F < F[i]).any(axis=1)).any()
    return keep


def knee(F):
    """(idx, method, dist, pw_dev) -- see module docstring."""
    fmin, fmax = F.min(0), F.max(0)
    rng = np.where(fmax - fmin > 0, fmax - fmin, 1.0)
    Z = (F - fmin) / rng
    n, m = Z.shape
    W = (fmax - F) / rng                       # pseudo-weights (sum 1)
    pw_dev = ((W / W.sum(1, keepdims=True) - 1.0 / m) ** 2).sum(1)
    anchors = np.array([int(np.argmin(Z[:, j])) for j in range(m)])
    A = np.column_stack([Z[anchors], np.ones(m)])
    if np.linalg.matrix_rank(A) == m:          # anchors affinely independent
        v = np.linalg.svd(A)[2][-1]            # null vector = (w, b)
        w, b = v[:m], v[m]
        d = (Z @ w + b) / np.linalg.norm(w)    # signed distance, utopia at b
        sgn = -np.sign(b) if b != 0 else 1.0
        i = int(np.argmax(sgn * d))
        return i, 'anchor-hyperplane', float(sgn * d[i]), float(pw_dev[i])
    i = int(np.argmin(pw_dev))
    return i, 'pseudo-weights', float('nan'), float(pw_dev[i])


def sim_z0(mode, src_idx, Xn):
    """Simulator at nominal z0 for one knee design (cached, ~15 s)."""
    cols = ['mode', 'src_idx', 'rn', 'rs', 'r_theta', 'Qs_MJ', 'sg_km',
            'q_stag_kW', 'q_shldr_kW', 'n_max_g', 't_final_s', 'captured',
            'feasible']
    hit = None
    if os.path.exists(SIM_CACHE):
        c = pd.read_csv(SIM_CACHE)
        m = ((c['mode'] == mode) & (c['src_idx'] == src_idx) &
             np.allclose(c[['rn', 'rs', 'r_theta']].to_numpy(), Xn,
                         atol=1e-9))
        if m.any():
            hit = c.loc[m].iloc[0]
    if hit is None:
        from capsule_opt.optimization.objectives import evaluate_shape
        t0 = time.time()
        res = evaluate_shape(*Xn, aero_scales=(1.0, 1.0),
                             entry_conditions={'gamma0_deg': -2.0,
                                               'Vo': 7830.0,
                                               't_max': 14400.0})
        o, det, tr = (res['objectives'], res['constraints']['details'],
                      res['traj'])
        hit = {'Qs_MJ': o['Qs'] / 1e6, 'sg_km': o['sg'] / 1e3,
               'q_stag_kW': det['q_stag_max'] / 1e3,
               'q_shldr_kW': det['q_shldr_max'] / 1e3,
               'n_max_g': det['n_max'], 't_final_s': tr['t'][-1],
               'captured': float('termination event' in tr['status']),
               'feasible': float(res['feasible'])}
        row = {'mode': mode, 'src_idx': src_idx, 'rn': Xn[0], 'rs': Xn[1],
               'r_theta': Xn[2], **hit}
        old = pd.read_csv(SIM_CACHE) if os.path.exists(SIM_CACHE) else None
        pd.concat([old, pd.DataFrame([row])])[cols].to_csv(
            SIM_CACHE, index=False)
        print(f"  simulator z0 [{mode} row {src_idx}] {time.time()-t0:.0f} s"
              f" | sg {hit['sg_km']:.0f} km | Qs {hit['Qs_MJ']:.1f} MJ/m2",
              flush=True)
    return hit


def design_row(mode, df_nd, i, eng):
    """Knee row: geometry, simulator nominal view, robust view, margins."""
    row = df_nd.loc[i]
    Xn = row[['rn', 'rs', 'r_theta']].to_numpy(float)
    Rn, Rs, tc = cap_params(*Xn)
    sim = sim_z0(mode, int(row.src_idx), Xn)
    sur = eng.det_predict(Xn[None, :])
    r = {'mode': mode, 'src_idx': int(row.src_idx),
         'rn': Xn[0], 'rs': Xn[1], 'r_theta': Xn[2],
         'Rn_m': Rn, 'Rs_m': Rs, 'theta_c_deg': tc, 'eta_V': row.eta_V,
         'Qs_z0_MJ': sim['Qs_MJ'], 'sg_z0_km': sim['sg_km'],
         'q_stag_z0_kW': sim['q_stag_kW'], 'q_shldr_z0_kW': sim['q_shldr_kW'],
         'n_max_z0_g': sim['n_max_g'], 't_final_z0_s': sim['t_final_s'],
         'captured_z0': sim['captured'],
         'sg_sur_z0_km': sur['sg'][0] / 1e3,     # diagnostic: slow-GP
         'E_Qs_MJ': row.E_Qs / 1e6, 'sd_Qs_MJ': row.sd_Qs_infl / 1e6,
         'CV_Qs': row.sd_Qs_infl / row.E_Qs, 'E_sg_km': row.E_sg / 1e3,
         'P_slow': row.P_slow}
    for o, L in LIMITS.items():
        e2s = row[f'E_{o}'] + 2 * row[f'sd_{o}_infl']
        u = 1e-3 if o.startswith('q') else 1.0
        r[f'{o}_E2S'] = e2s * u
        r[f'{o}_lim'] = L * u
        r[f'Pexc_{o}'] = row[f'Pexc_{o}']
    return r


def main():
    eng = MomentEngine()      # exact eta_V + surrogate diagnostics
    rows, nds = [], {}
    for mode in ('det', 'robust'):
        df = load_front(mode)
        F = F_min(df, mode)
        keep = strict_nd(F)
        nd = df.loc[keep].reset_index(names='src_idx')
        nd.to_csv(os.path.join(HERE, f'data/opt_{mode}_nd.csv'), index=False)
        nds[mode] = nd
        i, method, dist, pw = knee(F[keep])
        print(f"{mode}: {len(df)} -> {keep.sum()} strictly non-dominated | "
              f"knee = row {int(nd.loc[i, 'src_idx'])} "
              f"({method}, dist {dist:.3f}, pw_dev {pw:.4f})", flush=True)
        if mode == 'det':    # slow-GP extrapolation diagnostic (see docstring)
            print(f"  det front: surrogate sg(z0) above the 20,015 km "
                  f"great-circle max on {(df.sg > 20.015e6).sum()}"
                  f"/{len(df)} designs", flush=True)
        r = design_row(mode, nd, i, eng)
        r.update(knee_method=method, knee_dist=dist, pw_dev=pw,
                 n_front=int(keep.sum()))
        rows.append(r)

    # matched-robust counterpart of the det knee: nearest robust design in
    # the (E[Qs], E[sg]) plane (normalized by robust-front spreads) -- the
    # fig2 match window (+2%/-2%) misses the knee by 0.05% on E[sg], so the
    # parameter-free nearest design is used and its deltas reported.
    kdet = rows[0]
    m = nds['robust']
    qk = kdet['E_Qs_MJ'] * 1e6
    sk = kdet['E_sg_km'] * 1e3
    dist = np.hypot((m.E_Qs - qk) / (m.E_Qs.max() - m.E_Qs.min()),
                    (m.E_sg - sk) / (m.E_sg.max() - m.E_sg.min()))
    i_m = int(dist.idxmin())
    r = design_row('robust_m', m, i_m, eng)
    r['dE_Qs_pct'] = 100 * (r['E_Qs_MJ'] / kdet['E_Qs_MJ'] - 1)
    r['dE_sg_pct'] = 100 * (r['E_sg_km'] / kdet['E_sg_km'] - 1)
    r['sd_cut_pct'] = 100 * (1 - r['sd_Qs_MJ'] / kdet['sd_Qs_MJ'])
    rows.append(r)
    print(f"robust matched counterpart: row {int(r['src_idx'])} "
          f"(dE[Qs] {r['dE_Qs_pct']:+.1f}%, dE[sg] {r['dE_sg_pct']:+.1f}%, "
          f"sigma_Qs cut {r['sd_cut_pct']:.1f}%)", flush=True)

    out = os.path.join(HERE, 'data/knee_designs.csv')
    pd.DataFrame(rows).to_csv(out, index=False)
    cols = [('det knee', rows[0]), ('robust knee', rows[1])]
    if len(rows) > 2:
        cols.append(('robust matched', rows[2]))
    f2 = lambda v: f"{v:8.3f}"
    f3 = lambda v: f"{v:8.4f}"
    TAB = [
        ('Rn [m]', f2, 'Rn_m'), ('Rs [m]', f2, 'Rs_m'),
        ('theta_c [deg]', f2, 'theta_c_deg'),
        ('eta_V [-]', f3, 'eta_V'),
        ('Qs at z0, sim [MJ/m2]', f2, 'Qs_z0_MJ'),
        ('sg at z0, sim [km]', f2, 'sg_z0_km'),
        ('q_stag at z0, sim [kW/m2]', f2, 'q_stag_z0_kW'),
        ('q_shldr at z0, sim [kW/m2]', f2, 'q_shldr_z0_kW'),
        ('n_max at z0, sim [g]', f2, 'n_max_z0_g'),
        ('E[Qs] [MJ/m2]', f2, 'E_Qs_MJ'),
        ('sigma_Qs,infl [MJ/m2]', f2, 'sd_Qs_MJ'),
        ('CV(Qs) [%]', f2, 'CV_Qs'),
        ('E[sg] [km]', f2, 'E_sg_km'),
        ('q_stag E+2s [kW/m2] (lim 700)', f2, 'q_stag_max_E2S'),
        ('q_shldr E+2s [kW/m2] (lim 1000)', f2, 'q_shldr_max_E2S'),
        ('n_max E+2s [g] (lim 5)', f2, 'n_max_E2S'),
        ('P(q_stag > lim) [%]', f2, 'Pexc_q_stag_max'),
        ('P(q_shldr > lim) [%]', f2, 'Pexc_q_shldr_max'),
        ('P(n_max > lim) [%]', f3, 'Pexc_n_max'),
        ('P(slow) [%]', f2, 'P_slow'),
        ('sigma_Qs cut vs det knee [%]', f2, 'sd_cut_pct'),
    ]
    head = f"{'quantity':<34}" + ''.join(f"{c:>15}" for c, _ in cols)
    print("\n===== knee-point designs (robust lens, n_z=2000) =====",
          flush=True)
    print(head, flush=True)
    for lab, f, k in TAB:
        cells = []
        for _, rr in cols:
            v = rr.get(k)
            cells.append('      --' if v is None or (isinstance(v, float)
                                                     and np.isnan(v))
                         else f(v * 100 if k.startswith(('CV', 'Pexc',
                                                        'P_slow')) else v))
        print(f"{lab:<34}" + ''.join(f"{c:>15}" for c in cells), flush=True)
    d, r = rows[0], rows[1]
    print(f"\ndiagnostic: surrogate sg(z0) {d['sg_sur_z0_km']:.0f} / "
          f"{r['sg_sur_z0_km']:.0f} km vs simulator {d['sg_z0_km']:.0f} / "
          f"{r['sg_z0_km']:.0f} km (slow-GP extrapolation; det knee "
          f"overpredicted {d['sg_sur_z0_km'] / d['sg_z0_km']:.1f}x)",
          flush=True)
    print(f"saved {out} (+ simulator cache {SIM_CACHE})", flush=True)


if __name__ == '__main__':
    main()
