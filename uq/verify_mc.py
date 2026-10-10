"""Verification Monte Carlo (Step 6): true-simulator check of the surrogate
moment engine at the paper's headline designs.

Designs (fixed order; argv[2] runs only the first K):
  0 det_knee     uq/data/knee_designs.csv mode=det        (sg(z0)-inflation case)
  1 rob_knee     uq/data/knee_designs.csv mode=robust
  2 rob_matched  uq/data/knee_designs.csv mode=robust_m
  3 det_fragile  det front row with max audit P_exc(q_shldr)   (audit headline)
  4 rob_floor    robust front row with min sd_Qs_infl           (sigma floor)

Each design is evaluated at N z-samples through the TRUE simulator
(t_max=14400, same dt as DOE/probe), staged to uq/data/verify_mc.csv
(resume-safe, keyed on design_id x isamp). The z-LHS is drawn once at
N_MAX=500 (seed 2026 -- independent of the CRN seed 42 and probe seed 101)
and sliced to N, so partial runs stay valid prefixes of the full MC.

The summary compares MC moments against the moment engine at report
precision (MomentEngine default: n_z=2000, CRN seed 42 -- the numbers used
in the paper): E/sd[Qs], E[sg], P_slow, P_cap, peak mean+2sd and exceedance
counts (rule of three: 0 exceedances in n runs -> P <= 3/n).

Usage (from capsule_opt/):  python -m uq.verify_mc [N=300] [N_DESIGN=5]
  smoke: python -m uq.verify_mc 2 1        knee-only: python -m uq.verify_mc 300 3
Timing: ~15-50 s per sim (probe average ~26 s) -> 5 designs x 300 ~ 10-12 h.
"""

import csv
import ctypes
import math
import os
import sys
import time

import numpy as np
import pandas as pd

from capsule_opt.optimization.objectives import evaluate_shape
from capsule_opt.uncertainty.variables import sample

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'data/verify_mc.csv')
SUMMARY = os.path.join(HERE, 'data/verify_mc_summary.csv')
T_MAX = 14400.0
SLOW_CUT = 2250.0        # build v2 diagnostics: log-t gap 2241-2297 s
N_MAX, SEED_MC = 500, 2026
ZC = ['gamma0_deg', 'V0', 'k_CD', 'k_CL']
PEAKS = ['q_stag_max', 'q_shldr_max', 'n_max']
PEAK_UNIT = {'q_stag_max': 1e-3, 'q_shldr_max': 1e-3, 'n_max': 1.0}  # -> kW/m2, g
COLS = ['design_id', 'isamp', 'rn', 'rs', 'r_theta'] + ZC + \
       ['Qs', 'sg', 'q_stag_max', 'q_shldr_max', 'n_max', 'feasible',
        't_final', 'captured']


def load_designs(k):
    """First K verification designs as (name, design tuple, provenance)."""
    knee = pd.read_csv(os.path.join(HERE, 'data/knee_designs.csv'))
    by = knee.set_index('mode')
    kn = lambda m: f"knee_designs.csv mode={m} (front row {int(by.loc[m, 'src_idx'])})"
    out = [('det_knee', (by.loc['det', 'rn'], by.loc['det', 'rs'],
                         by.loc['det', 'r_theta']), kn('det')),
           ('rob_knee', (by.loc['robust', 'rn'], by.loc['robust', 'rs'],
                         by.loc['robust', 'r_theta']), kn('robust')),
           ('rob_matched', (by.loc['robust_m', 'rn'], by.loc['robust_m', 'rs'],
                            by.loc['robust_m', 'r_theta']), kn('robust_m'))]
    if k >= 4:
        d = pd.read_csv(os.path.join(HERE, 'data/opt_det_audit.csv'))
        r = d.loc[d.aud_Pexc_q_shldr_max.idxmax()]
        out.append(('det_fragile', (r.rn, r.rs, r.r_theta),
                    f"det front row {int(r.name)}, max audit Pexc(q_shldr) "
                    f"= {r.aud_Pexc_q_shldr_max:.3f}"))
    if k >= 5:
        d = pd.read_csv(os.path.join(HERE, 'data/opt_robust.csv'))
        r = d.loc[d.sd_Qs_infl.idxmin()]
        out.append(('rob_floor', (r.rn, r.rs, r.r_theta),
                    f"robust front row {int(r.name)}, min sd_Qs_infl "
                    f"= {r.sd_Qs_infl / 1e6:.1f} MJ/m2"))
    return out[:k]


def run_point(design, z):
    """One true-simulator evaluation (probe convention, g's dropped)."""
    rn, rs, rt = design
    g, V, kcd, kcl = z
    try:
        r = evaluate_shape(rn, rs, rt, aero_scales=(kcd, kcl),
                           entry_conditions={'gamma0_deg': g, 'Vo': V,
                                             't_max': T_MAX})
        det = r['constraints']['details']
        tr = r['traj']
        return [r['objectives']['Qs'], r['objectives']['sg'],
                det['q_stag_max'], det['q_shldr_max'], det['n_max'],
                float(r['feasible']), float(tr['t'][-1]),
                float('termination event' in tr['status'])]
    except Exception as e:
        print(f"  point FAILED: {type(e).__name__}: {e}", flush=True)
        return [math.nan] * 8


def summarize(k):
    """MC moments vs moment engine (n_z=2000 CRN) per design; save table."""
    from uq.engine import LIMITS, MomentEngine

    df = pd.read_csv(OUT)
    designs = load_designs(k)
    eng = MomentEngine()                       # n_z=2000, seed 42 (paper)
    m = eng.moments(np.array([list(d[1]) for d in designs]))
    rows = []
    for i, (name, _, prov) in enumerate(designs):
        d = df[df.design_id == i]
        cap = d[(d.captured == 1) & d.Qs.notna()]
        n, nc = len(d), len(cap)
        slow = (cap.t_final > SLOW_CUT)
        q = cap.Qs / 1e6
        se = q.std(ddof=1) / math.sqrt(nc) if nc > 1 else float('nan')
        sgk = cap.sg / 1e3
        print(f"\n-- {name}  [{prov}]", flush=True)
        print(f"   N={n}  captured {nc}  (P_cap sur {m['P_cap'][i]:.3f})"
              f"  slow {int(slow.sum())} = {slow.mean():.1%}"
              f"  (P_slow sur {m['P_slow'][i]:.3f})", flush=True)

        def line(lab, mc, se_, sur, unit='', scale=1.0):
            s_ = sur * scale
            dp = 100 * (mc - s_) / s_ if s_ else float('nan')
            z_ = (mc - s_) / se_ if se_ else float('nan')
            print(f"   {lab:<14}{mc:9.1f} ± {se_:5.1f} SE   surrogate "
                  f"{s_:9.1f} {unit:<7} diff {dp:+6.1f}% ({z_:+4.1f} SE)",
                  flush=True)
            rows.append({'design': name, 'quantity': lab, 'mc': mc,
                         'mc_se': se_, 'surrogate': s_, 'diff_pct': dp,
                         'diff_in_se': z_})

        line('E[Qs]', q.mean(), se, m['E_Qs'][i], 'MJ/m2', 1e-6)
        sdp, sdi = m['sd_Qs_plain'][i] / 1e6, m['sd_Qs_infl'][i] / 1e6
        print(f"   {'sd[Qs]':<14}{q.std(ddof=1):9.1f}            surrogate"
              f" plain {sdp:7.1f}  infl {sdi:7.1f}  MJ/m2", flush=True)
        rows.append({'design': name, 'quantity': 'sd[Qs]',
                     'mc': q.std(ddof=1), 'mc_se': float('nan'),
                     'surrogate': sdi, 'diff_pct': 100 * (q.std(ddof=1) - sdi)
                     / sdi, 'diff_in_se': float('nan')})
        line('E[sg]', sgk.mean(),
             sgk.std(ddof=1) / math.sqrt(nc) if nc > 1 else float('nan'),
             m['E_sg'][i], 'km', 1e-3)
        for o in PEAKS:
            u = PEAK_UNIT[o]
            v = cap[o].to_numpy() * u
            lim = LIMITS[o] * u
            e2s = (m[f'E_{o}'][i] + 2 * m[f'sd_{o}_infl'][i]) * u
            exc = int((cap[o] > LIMITS[o]).sum())
            print(f"   {o:<14}mean {v.mean():6.1f}  MC max {v.max():6.1f}"
                  f"  MC mean+2sd {v.mean() + 2 * v.std(ddof=1):6.1f}"
                  f"  (lim {lim:.0f}; sur E+2s {e2s:6.1f})", flush=True)
            ub = min(1.0, 3.0 / nc) if nc else float('nan')
            print(f"   {'':14}exceed {exc}/{nc} = {exc / nc if nc else 0:.1%}"
                  f" (0 -> P<= {ub:.1%})   sur Pexc "
                  f"{m[f'Pexc_{o}'][i]:.3f}", flush=True)
            rows.append({'design': name, 'quantity': f'Pexc_{o}',
                         'mc': exc / nc if nc else float('nan'),
                         'mc_se': float('nan'),
                         'surrogate': m[f'Pexc_{o}'][i],
                         'diff_pct': float('nan'), 'diff_in_se': float('nan')})

    pd.DataFrame(rows).to_csv(SUMMARY, index=False)
    print(f"\nsaved {SUMMARY}", flush=True)


def main():
    if os.name == 'nt':  # keep awake while running
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000001)
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    assert n <= N_MAX, f"N <= {N_MAX} (z-LHS is fixed at N_MAX and sliced)"
    designs = load_designs(k)
    Z = sample(N_MAX, SEED_MC)[:n]
    done = set()
    if os.path.exists(OUT):
        with open(OUT, newline='') as f:
            for r in csv.DictReader(f):
                done.add((int(r['design_id']), int(r['isamp'])))
    todo = [(i, j) for i in range(k) for j in range(n)
            if (i, j) not in done]
    print(f"VERIFY-MC: {k} designs x {n} samples, {len(todo)} to run "
          f"(~15-50 s each; 5x300 ~ 10-12 h), t_max={T_MAX:.0f} s, "
          f"z-seed {SEED_MC}", flush=True)
    for i, (name, _, prov) in enumerate(designs):
        print(f"  design {i} {name}: {prov}", flush=True)

    t0 = time.time()
    newf = not os.path.exists(OUT)
    with open(OUT, 'a', newline='') as f:
        w = csv.writer(f)
        if newf:
            w.writerow(COLS)
        for c, (i, j) in enumerate(todo, 1):
            name, design, _ = designs[i]
            w.writerow([i, j] + list(design) + list(Z[j]) +
                       run_point(design, Z[j]))
            f.flush()
            if c % 10 == 0 or c == len(todo):
                el = time.time() - t0
                print(f"  {c}/{len(todo)} done  {el / 60:.0f} min elapsed, "
                      f"ETA {(el / c * (len(todo) - c)) / 3600:.1f} h",
                      flush=True)
    summarize(k)


if __name__ == '__main__':
    main()
