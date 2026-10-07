"""Re-run the 202 horizon-mislabelled DOE rows at t_max=14400 s, fine dt.

uq/doe_descriptors.csv (extended horizon, coarse dt) proved 202 rows capture
between 7219 and 13474 s -- past the 7200 s horizon under which
doe_results.csv was produced. Their stored objectives are censored-branch
values: partial heat load Qs, mid-skip ground range sg, possibly stale
peaks. This script re-evaluates exactly those rows with uq.doe.run_point's
plumbing (evaluate_shape, default fine dt_max=2 s) plus t_max=14400:

  1. stage each row in uq/doe_rerun202.csv (append + flush + resume,
     doe.py conventions), adding fine-dt verification columns
     (t_final, captured_fine);
  2. once all 202 are staged and none NaN, splice the new output fields
     into doe_results.csv in place -- input fields byte-identical, the
     other 1598 rows untouched; backup first as uq/doe_results_pre202.csv.

Usage: python -m uq.doe_rerun202 [limit]
  limit : optional int -- run only the first <limit> pending rows and do
          NOT merge (smoke test). Full run ~50 s/row x 202 ~ 2.8 h
          (measured; these are the slowest trajectories in the DOE).
"""

import csv
import ctypes
import math
import os
import shutil
import statistics
import sys
import time

import pandas as pd

from capsule_opt.optimization.objectives import evaluate_shape
from uq.doe import IN_COLS, OUT_COLS

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, 'data/doe_results.csv')
DESCRIPTORS = os.path.join(HERE, 'data/doe_descriptors.csv')
STAGE = os.path.join(HERE, 'data/doe_rerun202.csv')
BACKUP = os.path.join(HERE, 'data/doe_results_pre202.csv')
T_MAX = 14400.0     # diag_censor_tail.py: slowest capture 13474 s + margin
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else None

STAGE_COLS = ['idx'] + IN_COLS + OUT_COLS + ['t_final', 'captured_fine']


def run_point(row):
    """uq.doe.run_point verbatim + t_max=14400 + fine-dt verification."""
    rn, rs, rt, g, V, kcd, kcl = row
    try:
        r = evaluate_shape(rn, rs, rt, aero_scales=(kcd, kcl),
                           entry_conditions={'gamma0_deg': g, 'Vo': V,
                                             't_max': T_MAX})
        det = r['constraints']['details']
        tr = r['traj']
        vals = [r['objectives']['Qs'], r['objectives']['sg'],
                det['q_stag_max'], det['q_shldr_max'], det['n_max'],
                float(r['feasible'])] + \
               [float(r['constraints'][k]) for k in
                ('g_q_stag', 'g_q_shldr', 'g_n_max', 'g_pitch_stable',
                 'g_trim_on_nose')]
        vals = [v if isinstance(v, float) else float(v) for v in vals]
        vals += [float(tr['t'][-1]),
                 float('termination event' in tr['status'])]
        return vals
    except Exception as e:  # record, keep going; merge refuses NaN rows
        print(f"  point FAILED: {type(e).__name__}: {e}", flush=True)
        return [math.nan] * (len(OUT_COLS) + 2)


def merge(targets, t_coarse):
    """Splice staged outputs into doe_results.csv (backup + atomic)."""
    with open(STAGE, newline='') as f:
        srows = list(csv.reader(f))
    assert srows[0] == STAGE_COLS, 'unexpected stage header'
    assert all(len(r) == len(STAGE_COLS) for r in srows[1:])
    staged = {int(r[0]): r for r in srows[1:]}
    assert set(staged) >= set(targets), 'stage file incomplete'

    with open(RESULTS, newline='') as f:
        lines = list(csv.reader(f))
    header, data = lines[0], lines[1:]
    assert header == IN_COLS + OUT_COLS, 'unexpected doe_results header'
    assert len(data) == 1800, f'expected 1800 data rows, got {len(data)}'

    nan_idx, nc_idx = [], []
    dq, ds, dtf = [], [], []
    max_in = 0.0
    for i in targets:
        s = staged[i]
        if any(v.lower() == 'nan' for v in s[8:19]):
            nan_idx.append(i)
            continue
        # pandas' C parser can be ~1 ULP off Python float(); tolerance, not ==
        ins = max(abs(float(a) - float(b)) / max(abs(float(b)), 1e-300)
                  for a, b in zip(s[1:8], data[i][:7]))
        max_in = max(max_in, ins)
        assert ins < 1e-9, f'idx {i}: stage/results input mismatch (rel {ins:.2e})'
        if float(s[-1]) == 0.0:
            nc_idx.append(i)
        oqs, osg, nqs, nsg = (float(data[i][7]), float(data[i][8]),
                              float(s[8]), float(s[9]))
        dq.append((nqs - oqs) / oqs)
        ds.append((nsg - osg) / osg)
        if float(s[-1]) == 1.0:
            dtf.append(abs(float(s[-2]) - t_coarse[i]) / t_coarse[i])
        data[i] = data[i][:7] + s[8:19]
    if nan_idx:
        print(f"MERGE ABORTED: NaN objectives in staged rows {nan_idx}. "
              f"Delete those lines from {os.path.basename(STAGE)} and "
              f"re-run.", flush=True)
        return

    if not os.path.exists(BACKUP):
        shutil.copy2(RESULTS, BACKUP)
        print(f"  backup written: {os.path.basename(BACKUP)}", flush=True)
    tmp = RESULTS + '.tmp'
    with open(tmp, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(data)
    os.replace(tmp, RESULTS)
    with open(RESULTS, newline='') as f:
        check = list(csv.reader(f))
    assert len(check) == 1801 and check[0] == header

    print(f"MERGED: {len(targets)} rows updated in doe_results.csv "
          f"(inputs byte-identical, {1800 - len(targets)} rows untouched).",
          flush=True)
    # both Qs and sg are monotone along a trajectory -> must not decrease
    print(f"  input round-trip max rel diff (stage vs results): {max_in:.2e}",
          flush=True)
    print(f"  dQs/Qs: median {statistics.median(dq)*100:+.1f}%, "
          f"min {min(dq)*100:+.1f}%, max {max(dq)*100:+.1f}%", flush=True)
    print(f"  dsg/sg: median {statistics.median(ds)*100:+.1f}%, "
          f"min {min(ds)*100:+.1f}%, max {max(ds)*100:+.1f}%", flush=True)
    if dtf:
        print(f"  t_final(fine) vs descriptors(coarse): max rel diff "
              f"{max(dtf)*100:.2f}%", flush=True)
    if nc_idx:
        print(f"  WARNING: {len(nc_idx)} re-run rows did NOT capture at fine "
              f"dt (idx {nc_idx}) -- descriptors disagreed; their values are "
              f"censored at {T_MAX:.0f} s.", flush=True)


def main():
    if os.name == 'nt':  # keep awake while running
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000001)

    desc = pd.read_csv(DESCRIPTORS)
    df = pd.read_csv(RESULTS)
    assert len(desc) == len(df), 'descriptors/results row-count mismatch'
    mis = desc[(desc['captured'] == 1) & (desc['t_final'] > 7200.0)]
    targets = [int(i) for i in mis['idx']]
    assert len(targets) == 202, \
        f'expected 202 mislabelled rows, got {len(targets)}'
    t_coarse = {int(r['idx']): float(r['t_final']) for _, r in desc.iterrows()}
    print(f"RERUN202: {len(targets)} rows (captured, t_final>7200 s per "
          f"descriptors), t_max={T_MAX:.0f} s, fine dt", flush=True)

    done = set()
    if os.path.exists(STAGE):
        with open(STAGE, newline='') as f:
            srows = list(csv.reader(f))
        done = {int(r[0]) for r in srows[1:]}
    pending = [i for i in targets if i not in done]
    print(f"  staged: {len(done)}, pending: {len(pending)}", flush=True)
    if LIMIT is not None:
        print(f"  SMOKE mode: limit={LIMIT}, will NOT merge", flush=True)

    n_todo = min(len(pending), LIMIT) if LIMIT is not None else len(pending)
    t0 = time.time()
    n_new = 0
    if n_todo:
        new_file = not os.path.exists(STAGE)
        with open(STAGE, 'a', newline='') as f:
            w = csv.writer(f)
            if new_file:
                w.writerow(STAGE_COLS)
            for i in pending:
                if n_new >= n_todo:
                    break
                x = [float(df.at[i, c]) for c in IN_COLS]
                w.writerow([i] + x + run_point(x))
                f.flush()
                n_new += 1
                if n_new % 10 == 0:
                    el = time.time() - t0
                    eta = el / n_new * (n_todo - n_new)
                    print(f"  {n_new}/{n_todo}  elapsed {el/60:.1f} min  "
                          f"ETA {eta/60:.1f} min", flush=True)
        el = time.time() - t0
        print(f"  {n_new} row(s) run in {el/60:.1f} min "
              f"({el / max(n_new, 1):.0f} s/row)", flush=True)

    if LIMIT is not None:
        print(f"SMOKE DONE: staged, NOT merged. Full run: re-launch without "
              f"the limit argument.", flush=True)
        return
    if n_new < len(pending):
        print(f"INCOMPLETE: staged {len(done) + n_new}/{len(targets)}; "
              f"re-run this script to resume. No merge.", flush=True)
        return
    merge(targets, t_coarse)


if __name__ == '__main__':
    main()
