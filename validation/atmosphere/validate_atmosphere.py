"""V0 - US-76 atmosphere verification (thesis section 4.5).

Checks published layer-base anchors, ideal-gas consistency of the analytic
0-86 km model, the continuity of the 86 km join to the tabulated extension,
density monotonicity and extension-node fidelity.

Run from the repo root:  python -m validation.atmosphere.validate_atmosphere
"""

import numpy as np

from validation import utils
from capsule_opt import config
from capsule_opt.atmosphere.atmosphere import atmosphere

R = config.r_earth


def h_geometric(H_km):
    """Geometric altitude [m] corresponding to a geopotential altitude [km]."""
    H = H_km * 1000.0
    return H * R / (R - H)


def main():
    d = utils.load_numeric("US76_reference_table.csv",
                           ["altitude_km", "T_K", "p_Pa", "rho_kg_m3"])

    # --- 1. published anchors (layer bases; sea-level triple) -------------
    ok, worst = True, 0.0
    for H_km, T_ref, p_ref, rho_ref in zip(d["altitude_km"], d["T_K"],
                                           d["p_Pa"], d["rho_kg_m3"]):
        atm = atmosphere(h_geometric(H_km))
        for val, ref, rel_band, tag in ((atm["T"], T_ref, 1e-4, "T"),
                                        (atm["P"], p_ref, 1e-3, "p"),
                                        (atm["rho"], rho_ref, 3e-3, "rho")):
            if not np.isfinite(ref):
                continue
            rel = abs(val - ref) / abs(ref)
            worst = max(worst, rel)
            ok &= rel <= rel_band
    utils.report("V0-anchors", "US-76 layer-base anchors (T/p/rho)",
                 ok, f"worst rel dev {worst:.2e}")

    # --- 2. ideal-gas consistency of the analytic model (0-86 km) --------
    hs = np.linspace(0.0, 85_999.0, 200)
    rel = max(abs(atmosphere(h)["P"] /
                  (atmosphere(h)["rho"] * config.R_air * atmosphere(h)["T"]) - 1.0)
              for h in hs)
    utils.report("V0-gas-law", "p = rho R T below the homopause",
                 rel < 1e-6, f"max rel dev {rel:.2e}")

    # --- 3. continuity of the 86 km join ---------------------------------
    a, b = atmosphere(85_999.0), atmosphere(86_001.0)
    steps = {k: abs(b[k] - a[k]) / a[k] for k in ("T", "P", "rho")}
    ok = all(v < 3e-3 for v in steps.values())
    utils.report("V0-join86", "continuity at the 86 km table join",
                 ok, " ".join(f"{k} {v:.2e}" for k, v in steps.items()))

    # --- 4. density monotone over the full model -------------------------
    hs = np.linspace(0.0, 120_000.0, 400)
    rho = np.array([atmosphere(h)["rho"] for h in hs])
    utils.report("V0-monotone", "rho strictly decreasing 0-120 km",
                 bool(np.all(np.diff(rho) < 0)))

    # --- 5. extension-node fidelity --------------------------------------
    worst = 0.0
    for h, T_ref, p_ref, rho_ref in ((86_000, config.T_EXT[0], config.P_EXT[0], config.RHO_EXT[0]),
                                     (100_000, config.T_EXT[28], config.P_EXT[28], config.RHO_EXT[28]),
                                     (120_000, config.T_EXT[-1], config.P_EXT[-1], config.RHO_EXT[-1])):
        atm = atmosphere(h)
        for val, ref in ((atm["T"], T_ref), (atm["P"], p_ref), (atm["rho"], rho_ref)):
            worst = max(worst, abs(val - ref) / ref)
    utils.report("V0-ext-nodes", "86-120 km interpolation hits table nodes",
                 worst < 5e-4, f"worst rel dev {worst:.2e}")

    # --- 6. sea-level speed of sound -------------------------------------
    a0 = atmosphere(0.0)["a"]
    utils.report("V0-sound", "sea-level speed of sound",
                 abs(a0 - 340.294) < 0.1, f"a0 = {a0:.3f} m/s (ref 340.294)")


if __name__ == "__main__":
    main()
