"""V1d - Mach-regime blend continuity (Eq 3.48-3.49 implementation).

The low/high hypersonic databases are cubic-blended over M in [5, 12] with
the continuous weight w = 3 s^2 - 2 s^3, s = (M - 5) / 7. Checks the weight
endpoints, monotone transition, per-panel Cp continuity across the overlap,
and finiteness of the full law across the production Mach range.

Run from the repo root:
    python -m validation.aerodynamics.validate_blend_continuity
"""

import numpy as np

from validation import utils
from capsule_opt.aerodynamics.lim import high_hyp_weight, lim_pressure_coefficient
from capsule_opt.properties import PART_AFTERBODY


def main():
    # --- 1. weight endpoints + monotone transition ------------------------
    w5, w12 = high_hyp_weight(5.0), high_hyp_weight(12.0)
    ms = np.linspace(5.0, 12.0, 141)
    w = np.array([high_hyp_weight(m) for m in ms])
    ok = (abs(w5) < 1e-12 and abs(w12 - 1.0) < 1e-12
          and bool(np.all(np.diff(w) >= -1e-12))
          and 0.45 < high_hyp_weight(8.5) < 0.55)     # smooth midpoint
    utils.report("V1d-weight", "blend weight endpoints + monotone + midpoint",
                 bool(ok), f"w(5)={w5:.1e}, w(12)={w12:.12f}, w(8.5)={high_hyp_weight(8.5):.3f}")

    # --- 2. Cp continuity across the overlap for sample inclinations ------
    ok_all, worst = True, 0.0
    for theta in (10.0, 30.0, 55.0, 70.0):
        ms = np.linspace(4.90, 12.10, 721)
        cp = np.array([lim_pressure_coefficient(np.array([theta]),
                                                np.array([PART_AFTERBODY]), m)[0]
                       for m in ms])
        step = np.max(np.abs(np.diff(cp)))          # grid is 0.01 in M
        worst = max(worst, step)
        ok_all &= step < 5e-3 and bool(np.all(np.isfinite(cp)))
    utils.report("V1d-continuity", "afterbody Cp continuous across M in [5,12]",
                 bool(ok_all), f"max |dCp| per 0.01 in M = {worst:.2e}")

    # --- 3. finiteness over the full production range ---------------------
    ok = True
    for M in np.linspace(3.0, 30.0, 271):
        th = np.linspace(-85.0, 85.0, 171)
        cp = lim_pressure_coefficient(th, np.full(th.shape, PART_AFTERBODY), M)
        ok &= bool(np.all(np.isfinite(cp)))
    utils.report("V1d-finite", "law finite over theta in [-85,85], M in [3,30]",
                 bool(ok))


if __name__ == "__main__":
    main()
