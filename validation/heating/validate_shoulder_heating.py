"""V3b - Shoulder-heating correlation vs thesis Table 4.1 source data.

Reproduces the thesis's section 4.3.2 check: the implemented linear
shoulder-to-stagnation flux ratio (Dirkx & Mooij fit) evaluated at the
Marvin & Sinclair / Jones / Wadham source conditions, with the
nose-sphere half-angle held at the Apollo value exactly as the thesis did.

Pass band: 12 % relative (the thesis reported deviations up to 11.7 %).
Trend checks: ratio rises with alpha, falls with shoulder radius.

Run from the repo root:
    python -m validation.heating.validate_shoulder_heating
"""

import numpy as np

import validation.apollo_case as apollo
from validation import utils
from capsule_opt.heating.heating import shoulder_heat_flux

REL_BAND = 0.12


def main():
    theta_sp1 = apollo.apollo_props()["theta_sp1"]      # Apollo nose-sphere
    d = utils.load_numeric("shoulder_ratio_reference.csv",
                           ["mach", "alpha_deg", "rs_over_rm", "reference_ratio"])

    computed = np.array([
        shoulder_heat_flux(1.0, M, a, theta_sp1, rs, 1.0)
        for M, a, rs in zip(d["mach"], d["alpha_deg"], d["rs_over_rm"])
    ])
    ref = d["reference_ratio"]
    rel = np.abs(computed - ref) / ref
    worst_i = int(np.argmax(rel))
    utils.report("V3b-shoulder-ratio",
                 "shoulder/stagnation ratio vs thesis Table 4.1 sources",
                 bool(np.all(rel <= REL_BAND)),
                 f"worst {100*rel[worst_i]:.1f}% at M={d['mach'][worst_i]:g}, "
                 f"alpha={d['alpha_deg'][worst_i]:g} deg")

    # trend checks (same-source families)
    j = (d["mach"] == 8.0)
    ms = np.isclose(d["mach"], 10.5) & np.isclose(d["alpha_deg"], 0.0)
    rising_alpha = computed[j][0] < computed[j][1] < computed[j][2]
    falling_rs = computed[ms][0] > computed[ms][1] > computed[ms][2]
    utils.report("V3b-trends", "ratio rises with alpha, falls with Rs",
                 bool(rising_alpha and falling_rs))

    for M, a, rs, r_ref, c in zip(d["mach"], d["alpha_deg"], d["rs_over_rm"],
                                  ref, computed):
        print(f"  M={M:<5g} alpha={a:<4g} Rs/Rm={rs:<5g}: "
              f"ref {r_ref:.2f}  computed {c:.3f}  ({100*(c-r_ref)/r_ref:+.1f}%)")


if __name__ == "__main__":
    main()
