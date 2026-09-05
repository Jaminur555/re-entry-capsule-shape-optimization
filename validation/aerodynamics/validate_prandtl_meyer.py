"""V1c - Prandtl-Meyer expansion branch: analytic properties.

No D&M curve exists for the lee side, so the check is against closed-form
gas dynamics: nu(M) anchors, Cp(theta=0)=0, monotone expansion, and the
vacuum-limit clamp Cp_vac = -2/(gamma M^2) beyond nu(M).

Run from the repo root:  python -m validation.aerodynamics.validate_prandtl_meyer
"""

import numpy as np

from validation import utils
from capsule_opt import config
from capsule_opt.aerodynamics.expansion import nu, cp_prandtl_meyer
from capsule_opt.aerodynamics.interpolation import cp_prandtl_interp

GAMMA = config.gamma


def main():
    # --- 1. nu(M) anchors (independent published values) ------------------
    d = utils.load_numeric("prandtl_meyer_reference.csv", ["mach", "nu_deg"])
    worst = 0.0
    for M, nu_ref in zip(d["mach"], d["nu_deg"]):
        if np.isfinite(nu_ref):
            worst = max(worst, abs(np.rad2deg(nu(M)) - nu_ref))
    utils.report("V1c-nu", "Prandtl-Meyer angle nu(M) vs published anchors",
                 worst < 0.15, f"max dev {worst:.3f} deg")

    # --- 2. closed-form Cp properties (theta_turn in RADIANS) -------------
    # Tolerances reflect the implementation's numerical epsilons: Cp(0)
    # carries ~1e-11 from the nu() root-find, and the clamp beyond nu(M)
    # sits ~2.5e-4 from the exact vacuum value (0.16% of |Cp_vac|).
    ok = True
    for M in (3.0, 5.0, 10.0, 25.0):
        th = np.linspace(0.0, 89.0, 179)
        cp = np.array([cp_prandtl_meyer(M, np.deg2rad(t)) for t in th])
        ok &= abs(cp[0]) < 1e-9                        # Cp(0) = 0
        ok &= bool(np.all(np.diff(cp) <= 1e-12))       # monotone expansion
        n_max_deg = np.rad2deg(nu(M))
        beyond = th > n_max_deg
        cp_vac = -2.0 / (GAMMA * M * M)
        if beyond.any():
            ok &= bool(np.all(np.abs(cp[beyond] - cp_vac) < 1e-3))  # clamp
    utils.report("V1c-closed-form",
                 "Cp(0)=0, monotone in theta, vacuum clamp beyond nu(M)",
                 bool(ok))

    # --- 3. production table vs closed form (interpolation error) ---------
    # The 2-deg theta grid cannot resolve the bend at the vacuum asymptote,
    # so points within 3 deg of nu(M) are assessed against a looser band.
    rng = np.random.default_rng(0)
    worst, worst_bend = 0.0, 0.0
    for M in rng.uniform(3.0, 25.0, 40):
        n_max = np.rad2deg(nu(M))
        for t in rng.uniform(0.0, min(89.0, n_max), 10):
            err = abs(cp_prandtl_interp(M, t)
                      - cp_prandtl_meyer(M, np.deg2rad(t)))
            if n_max - t < 3.0:
                worst_bend = max(worst_bend, err)
            else:
                worst = max(worst, err)
    utils.report("V1c-table", "cached PM table vs closed form (smooth region)",
                 worst < 1e-2, f"max|dCp| {worst:.2e} (2-deg theta x 1-Mach grid)")
    utils.report("V1c-table-bend",
                 "PM table near the vacuum asymptote (within 3 deg of nu)",
                 worst_bend < 1e-2, f"max|dCp| {worst_bend:.2e} (2-deg grid)")


if __name__ == "__main__":
    main()
