#!/usr/bin/env python
"""Convert the PlotDigitizer extraction of Dirkx (2017) into reference CSVs.

Source: ../../Extracted Data.txt (Figs 3.5/3.6, 3.9, 7.7, 7.17, Table 7.1)
Output: validation/reference_data/*.csv (one row per digitized point,
consumed by the V0-V4 scripts). Hygiene per README section 5.

Run:  python validation/build_reference_data.py
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import numpy as np

SRC = Path(__file__).resolve().parents[2] / "Extracted Data.txt"
NEWFIG717 = SRC.parent / "digitized_data.csv"
OUT = Path(__file__).resolve().parent / "reference_data"

# section starters (lower-cased prefix -> mode), checked before block parsing
SECTIONS = ("cd vs alpha", "cl vs alpha", "cm vs alpha", "l/d vs alpha",
            "cp vs aoa", "zoby", "table 7.1", "table 7.2", "trajectory")

RE_MACH_LABEL = re.compile(r"^M\s*=\s*([0-9.]+)\s*\((.*)\)\s*$", re.I)
RE_MACH_CP = re.compile(r"^M\s*=\s*(infinte|infinite|[0-9.]+)\s*:\s*$", re.I)
RE_ARRAY_ANY = re.compile(
    r"^(?P<name>[A-Za-z][\w ./]*)\s*(?P<unit>\([^)]*\))?\s*=\s*(?P<arr>\[.*\])\s*$")
RE_ARRAY = re.compile(
    r"^(?P<name>M|alpha|L/D|bank_angle|Qcs|n|t|h|theta|cp|Rm\s*/\s*Rn|Rm/Reff)"
    r"\s*(?P<unit>\([^)]*\))?\s*=\s*(?P<arr>\[.*\])\s*$", re.I)


def split_row(line):
    """'a,\tb' -> [a, b] as floats; [] if not exactly two numbers."""
    parts = [p.strip() for p in re.split(r"[,\t]+", line.strip()) if p.strip()]
    try:
        vals = [float(p) for p in parts]
    except ValueError:
        return []
    return vals if len(vals) == 2 else []


def parse_array(text):
    return [float(v) for v in text.strip().strip("[]").split(",") if v.strip()]


def source_of(label):
    lab = label.lower()
    if "local incl" in lab:
        return "lim"
    if "wind" in lab:
        return "wind_tunnel"
    return None


def is_section(line):
    low = line.lower()
    return any(low.startswith(s) for s in SECTIONS)


def main():
    lines = SRC.read_text(encoding="utf-8").splitlines()
    n = len(lines)

    aero_rows = []    # coeff, mach, source, alpha_deg, value
    cp_rows = []      # method, mach, theta_deg, cp
    zs_rows = []      # rs_over_rm, rm_over_rn, rm_over_reff
    table71 = []      # quantity, wind_tunnel, local_inclination, pct_diff
    traj_rows = []    # profile, aero_source, indep_name, indep_value, value
    entry_ics = {}

    mode, coeff = None, None
    method = None
    zs_key, zs_x = None, None
    traj_src = None
    pend = {}         # pending arrays for cp / trajectory pairing

    i = 0
    while i < n:
        s = lines[i].strip()
        low = s.lower()

        # ---------------- section switching ----------------
        if low.startswith(("cd vs alpha", "cl vs alpha", "cm vs alpha", "l/d vs alpha")):
            coeff = ("CD" if low.startswith("cd") else
                     "CL" if low.startswith("cl") else
                     "Cm" if low.startswith("cm") else "L/D")
            mode = "coeff"
            i += 1
            continue
        if low.startswith("cp vs aoa"):
            mode, method, pend = "cp", None, {}
            i += 1
            continue
        if low.startswith("zoby"):
            mode, zs_key, zs_x = "zs", None, None
            i += 1
            continue
        if low.startswith("table 7.1"):
            mode = "t71"
            i += 1
            continue
        if low.startswith("table 7.2"):
            mode = "t72"          # Shuttle: parsed over, not emitted
            i += 1
            continue
        if low.startswith("trajectory"):
            mode, traj_src = "traj", None
            i += 1
            continue
        if not s:
            i += 1
            continue

        # ---------------- coefficient sections (Fig 7.7) ----------------
        if mode == "coeff":
            m = RE_MACH_LABEL.match(s)
            if m:
                mach, src = float(m.group(1)), source_of(m.group(2))
                i += 1
                # two export styles: legacy 'alpha, Cd' header + pair rows,
                # or PlotDigitizer arrays 'alpha/x = [...]' + 'CD/Cl/Cm/y = [...]'
                pend_alpha, pend_vals = None, None
                while i < n:
                    t = lines[i].strip()
                    if RE_MACH_LABEL.match(t) or is_section(t):
                        break
                    am = RE_ARRAY_ANY.match(t)
                    if am:
                        key = re.sub(r"\s+", "", am.group("name").lower())
                        arr = parse_array(am.group("arr"))
                        if key.startswith("alpha") or key == "x":
                            pend_alpha = arr
                        else:
                            pend_vals = arr
                        if pend_alpha and pend_vals:
                            for a, v in zip(pend_alpha, pend_vals):
                                aero_rows.append((coeff, mach, src, a, v))
                            pend_alpha = pend_vals = None
                        i += 1
                        continue
                    if not t and not (pend_alpha or pend_vals):
                        i += 1
                        continue
                    vals = split_row(t)
                    if len(vals) == 2:
                        aero_rows.append((coeff, mach, src, vals[0], vals[1]))
                    i += 1
                continue

        # ---------------- Cp vs theta (Figs 3.5/3.6) ----------------
        if mode == "cp":
            if low.startswith("iii)"):
                method = "tangent_cone"
            elif low.startswith("ii)"):
                method = "tangent_wedge"
            elif low.startswith("i)"):
                method = "modified_newtonian"
            else:
                mm = RE_MACH_CP.match(s)
                if mm:
                    g = mm.group(1).lower()
                    pend = {"mach": "inf" if g.startswith("inf") else float(g),
                            "theta": None, "cp": None}
                else:
                    am = RE_ARRAY.match(s)
                    if am and pend:
                        key = re.sub(r"\s+", "", am.group("name").lower())
                        if key == "theta":
                            pend["theta"] = parse_array(am.group("arr"))
                        elif key == "cp":
                            pend["cp"] = parse_array(am.group("arr"))
                        if pend.get("theta") and pend.get("cp") and "mach" in pend:
                            order = np.argsort(pend["theta"])
                            th = np.clip(np.array(pend["theta"])[order], None, 90.0)
                            cpv = np.array(pend["cp"])[order]
                            mach = pend["mach"]
                            for a, b in zip(th, cpv):
                                cp_rows.append((method, mach, a, b))
                            pend = {"mach": mach, "theta": None, "cp": None}
            i += 1
            continue

        # ---------------- Zoby-Sullivan chart (Fig 3.9) ----------------
        if mode == "zs":
            mk = re.match(r"^Rs/Rm\s*=\s*([0-9.]+)\s*:", s)
            mx = re.match(r"^Rm\s*/\s*Rn\s*=", s)
            my = re.match(r"^Rm/Reff\s*=", s)
            if mk:
                zs_key = float(mk.group(1))
                zs_x = None
            elif mx:
                am = RE_ARRAY.match(s)
                if am:
                    zs_x = parse_array(am.group("arr"))
            elif my and zs_key is not None:
                am = RE_ARRAY.match(s)
                if am and zs_x:
                    ys = parse_array(am.group("arr"))
                    order = np.argsort(zs_x)
                    for x, y in zip(np.array(zs_x)[order], np.array(ys)[order]):
                        zs_rows.append((zs_key, x, y))
            i += 1
            continue

        # ---------------- Table 7.1 ----------------
        if mode == "t71":
            m1 = re.match(r"^(Stagnation-point heat load),\s*([0-9.]+)\s*MJ/m2,\s*"
                          r"([0-9.]+)\s*MJ/m2,\s*\+?(-?[0-9.]+)", s)
            m2 = re.match(r"^(Ground track length),\s*([0-9.]+)\s*km,\s*"
                          r"([0-9.]+)\s*km,\s*\+?(-?[0-9.]+)", s)
            m = m1 or m2
            if m:
                table71.append((m.group(1), float(m.group(2)),
                                float(m.group(3)), float(m.group(4))))
            i += 1
            continue

        # ---------------- Trajectory (Fig 7.17) ----------------
        if mode == "traj":
            if low.startswith("initial condition"):
                clean = s.replace("−", "-")
                for key, pat in [("h_km", r"h=([0-9.]+)km"), ("lon_deg", r"τ=([0-9.]+)"),
                                 ("lat_deg", r"δ=(-?[0-9.]+)"), ("V_R_km_s", r"V_R=([0-9.]+)"),
                                 ("gamma_deg", r"γ=(-?[0-9.]+)"), ("chi_deg", r"χ=([0-9.]+)"),
                                 ("mass_kg", r"m=([0-9.]+)")]:
                    mm = re.search(pat, clean)
                    if mm:
                        entry_ics[key] = float(mm.group(1))
                i += 1
                continue
            if low.startswith("wind-tunnel"):
                traj_src = "wind_tunnel"
                i += 1
                continue
            if low.replace(" ", "_").startswith("local_inclination"):
                traj_src = "lim"
                i += 1
                continue
            am = RE_ARRAY.match(s)
            if am and traj_src:
                name = re.sub(r"\s+", "", am.group("name").lower())
                arr = parse_array(am.group("arr"))
                profiles = {"alpha": "alpha_trim_deg", "l/d": "lift_to_drag",
                            "bank_angle": "bank_angle_deg",
                            "qcs": "heat_rate_w_m2", "n": "load_factor_g"}
                if name == "m":
                    pend["mach"] = arr
                elif name == "t":
                    pend["time_s"] = arr
                elif name == "h" and "time_s" in pend:
                    for t, h in zip(pend.pop("time_s"), arr):
                        traj_rows.append(("altitude_km", traj_src, "time_s", t, h))
                elif name in profiles and "mach" in pend:
                    for mm_, v in zip(pend["mach"], arr):
                        traj_rows.append((profiles[name], traj_src, "mach", mm_, v))
            i += 1
            continue

        i += 1

    # ---------------- emit CSVs ----------------
    OUT.mkdir(exist_ok=True)

    def write(name, header, rows, fmt=lambda v: f"{v:.12g}"):
        with open(OUT / name, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(header)
            for r in rows:
                w.writerow([fmt(v) if isinstance(v, float) else v for v in r])
        print(f"  {name:<45s} {len(rows):4d} rows")

    aero_rows.sort(key=lambda r: (r[0], r[1], r[2] or "", r[3]))
    write("apollo_flight_wind_tunnel_data.csv",
          ("coeff", "mach", "source", "alpha_deg", "value"), aero_rows)

    cp_rows.sort(key=lambda r: (r[0], np.inf if r[1] == "inf" else float(r[1]), r[2]))
    write("cp_methods_dm2017_digitized.csv",
          ("method", "mach", "theta_deg", "cp"), cp_rows)

    zs_rows.sort(key=lambda r: (r[0], r[1]))
    write("zoby_sullivan_chart_digitized.csv",
          ("rs_over_rm", "rm_over_rn", "rm_over_reff"), zs_rows)

    write("dm2017_table7_1.csv",
          ("quantity", "wind_tunnel", "local_inclination", "pct_diff"), table71)

    # Fig 7.17 re-digitization (PlotDigitizer CSV export) supersedes the
    # Extracted Data.txt trajectory blocks when present
    if NEWFIG717.exists():
        profmap = {"time_height": ("altitude_km", "time_s"),
                   "mach_alpha": ("alpha_trim_deg", "mach"),
                   "mach_ld": ("lift_to_drag", "mach"),
                   "mach_bank": ("bank_angle_deg", "mach"),
                   "mach_qcs": ("heat_rate_w_m2", "mach"),
                   "mach_load_factor": ("load_factor_g", "mach")}
        traj_rows = []
        with open(NEWFIG717, newline="", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                ds = row["dataset"].strip().lower()
                if ds in profmap and row.get("x") and row.get("y"):
                    p, ind = profmap[ds]
                    traj_rows.append((p, "lim", ind,
                                      float(row["x"]), float(row["y"])))
        traj_rows.sort(key=lambda r: (r[0], r[3]))

        # stray-dropout filter for the heat-rate curve: PlotDigitizer mis-clicks
        # sit far below the smooth curve; flag points >1.8x off a rolling
        # median (window +-2 pts) in log space
        idx = [i for i, r in enumerate(traj_rows) if r[0] == "heat_rate_w_m2"]
        y = np.log10(np.array([traj_rows[i][4] for i in idx]))
        keep = np.array([abs(y[k] - np.median(y[max(0, k - 2):k + 3]))
                         <= np.log10(1.8) for k in range(len(y))])
        dropped = [(traj_rows[i][3], traj_rows[i][4])
                   for i, k in zip(idx, range(len(y))) if not keep[k]]
        traj_rows = [r for i, r in enumerate(traj_rows) if i not in set(idx)] \
            + [traj_rows[i] for i, k in zip(idx, range(len(y))) if keep[k]]
        traj_rows.sort(key=lambda r: (r[0], r[3]))
        print(f"  Fig 7.17: {len(traj_rows)} rows from {NEWFIG717.name} "
              f"(supersedes Extracted Data.txt blocks)")
        if dropped:
            print(f"  heat-rate dropouts removed ({len(dropped)}): "
                  + ", ".join(f"M={m:.2f}:{v/1e3:.0f}kW" for m, v in dropped))

    write("apollo_entry_reference_dm2017.csv",
          ("profile", "aero_source", "indep_name", "indep_value", "value"), traj_rows)

    # thesis Table 4.1 shoulder-ratio source points (manual, not digitized)
    shoulder = [
        ("marvin_sinclair", 10.5, 0.0, 0.05, 1.35),
        ("marvin_sinclair", 10.5, 0.0, 0.15, 1.28),
        ("marvin_sinclair", 10.5, 0.0, 0.25, 1.20),
        ("jones", 8.0, 0.0, 0.00, 1.35),
        ("jones", 8.0, 15.0, 0.00, 1.55),
        ("jones", 8.0, 30.0, 0.00, 2.00),
        ("wadham", 10.0, 0.0, 0.05, 1.25),
        ("wadham", 10.0, 20.0, 0.05, 1.50),
        ("wadham", 10.0, 28.0, 0.05, 1.70),
    ]
    write("shoulder_ratio_reference.csv",
          ("source", "mach", "alpha_deg", "rs_over_rm", "reference_ratio"), shoulder)

    # ---------------- summary + digitization coherence check ----------------
    print("\n  curves per (coeff, mach, source):")
    from collections import Counter
    for (c, m, s), k in sorted(Counter((r[0], r[1], r[2]) for r in aero_rows).items()):
        print(f"    {c:4s} M={m:<4g} {s:<12s} {k:3d} pts")
    print("  cp curves per (method, mach):")
    for (meth, m), k in sorted(Counter((r[0], str(r[1])) for r in cp_rows).items()):
        print(f"    {meth:<18s} M={m:<5s} {k:3d} pts")
    print("  trajectory profiles per source:",
          dict(Counter((r[1], r[0]) for r in traj_rows)))
    print("  Table 7.1 rows:", len(table71), "| entry ICs:", entry_ics)

    # L/D vs CL/CD coherence (nearest-alpha) -- digitization self-check
    by_key = {}
    for c, m, s, a, v in aero_rows:
        by_key.setdefault((c, m, s), ([], []))
        by_key[(c, m, s)][0].append(a)
        by_key[(c, m, s)][1].append(v)
    worst = 0.0
    for (c, m, s) in list(by_key):
        if c != "L/D":
            continue
        for key in (("CL", m, s), ("CD", m, s)):
            if key not in by_key:
                break
        else:
            la, lv = by_key[("CL", m, s)]
            da, dv = by_key[("CD", m, s)]
            xa, xv = by_key[("L/D", m, s)]
            la, da = np.array(la), np.array(da)
            lo = max(la.min(), da.min())          # skip endpoints where the
            hi = min(la.max(), da.max())          # alpha grids do not overlap
            for a, ld in zip(xa, xv):
                if not (lo <= a <= hi):
                    continue
                cl = np.interp(a, la, np.array(lv))
                cd = np.interp(a, da, np.array(dv))
                worst = max(worst, abs(cl / cd - ld))
    print(f"\n  digitization coherence |L/D - CL/CD| max = {worst:.4f} "
          f"({'good' if worst < 0.01 else 'CHECK'})")


if __name__ == "__main__":
    main()
