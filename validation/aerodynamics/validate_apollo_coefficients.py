"""V2 - Apollo force/moment coefficients vs Dirkx (2017) Fig 7.7 (core check).

Exact-Apollo geometry + rcom CG (validation/apollo_case); alpha sweep 0-30
deg at M = 3 and 10, overlaid on the digitized LIM / wind-tunnel curves.
Conventions and the digitized y-scale caveat: validation/README.md sections
2 and 5. CD/CL are scored on curve shape after removing a best-fit constant
scale; L/D and Cm are compared directly.
Pass bands: shape |d| <= 0.05 (CD/CL), 0.01 (Cm), L/D <= 0.03 absolute.
Exception: CD-M3 uses 0.065 (same as V1b) -- documented beyond-detach method
difference (exact-TM saturation vs the Fig-3.5 correlation; detach boundary
lowest at M=3), not a digitization artifact.

Run from the repo root:
    python -m validation.aerodynamics.validate_apollo_coefficients
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import validation.apollo_case as apollo
from validation import utils
from capsule_opt.aerodynamics.solver import HypersonicAeroSolver

BANDS = {"CD": 0.05, "CL": 0.05, "Cm": 0.01, "L/D": 0.03}
# CD-M3: documented V1b beyond-detach method difference (see module docstring)
BAND_OVERRIDES = {("CD", 3.0): 0.065}
MACHS = (3.0, 10.0)


def digitized():
    """{coeff: {mach: {source: (alpha_deg, value)}}} from the reference CSV."""
    rows = utils.load_reference("apollo_flight_wind_tunnel_data.csv")
    out = {}
    for r in rows:
        c = r["coeff"]
        m = float(r["mach"])
        s = r["source"]
        out.setdefault(c, {}).setdefault(m, {}).setdefault(s, ([], []))
        out[c][m][s][0].append(float(r["alpha_deg"]))
        out[c][m][s][1].append(float(r["value"]))
    return {c: {m: {s: (np.array(a), np.array(v)) for s, (a, v) in d.items()}
                for m, d in mm.items()} for c, mm in out.items()}


def main():
    props = apollo.apollo_props()
    solver = HypersonicAeroSolver()
    scale = apollo.cm_ref_scale()
    ref = digitized()

    alphas = np.linspace(0.0, 30.0, 25)          # our (positive) convention
    utils.set_paper_style()

    for M in MACHS:
        res = [solver.coefficients(props, M, a) for a in alphas]
        ours = {"CD": np.array([r["CD"] for r in res]),
                "CL": np.array([r["CL"] for r in res]),
                "Cm": -np.array([r["Cm"] for r in res]) * scale,   # sign flip
                "L/D": np.array([r["L/D"] for r in res])}

        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        panels = [("CD", axes[0, 0]), ("CL", axes[0, 1]),
                  ("Cm", axes[1, 0]), ("L/D", axes[1, 1])]
        for coeff, ax in panels:
            a_dm, v_dm = ref[coeff][M]["lim"]        # alpha_DM negative
            interp = np.interp(-a_dm, alphas, ours[coeff])
            if coeff == "L/D":
                err = np.abs(interp - v_dm)
                detail = f"max|d|={err.max():.4f} (absolute)"
            else:
                k = float(np.median(interp / v_dm))          # fitted scale
                err = np.abs(interp / k - v_dm)
                detail = (f"shape max|d|={err.max():.4f}, fitted scale "
                          f"k={k:.3f}, abs max|d|={np.abs(interp - v_dm).max():.3f}")
            ok = bool(np.all(err <= BAND_OVERRIDES.get((coeff, M), BANDS[coeff])))
            utils.report(f"V2-{coeff}-M{M:g}",
                         f"Apollo {coeff} vs digitized D&M LIM (M={M:g})",
                         ok, detail)

            ax.plot(-a_dm, v_dm, "o", ms=4, label="D&M LIM (digitized)")
            w = ref[coeff][M].get("wind_tunnel")
            if w is not None:
                ax.plot(-w[0], w[1], "x", ms=5, label="wind tunnel")
            ax.plot(alphas, ours[coeff], "-", label="Present model")
            ax.set_xlabel(r"$\alpha$ [deg]")
            ax.set_ylabel(coeff)
            ax.legend(fontsize=8)
        fig.suptitle(f"Apollo aerodynamic coefficients, M = {M:g} (Fig. 7.7)")
        utils.save_fig(fig, f"v2_apollo_coefficients_M{M:g}.png")

    # method-behaviour check: CD over-prediction vs WT near alpha=0, M=10
    res0 = solver.coefficients(props, 10.0, 0.0)
    wt0 = float(ref["CD"][10.0]["wind_tunnel"][1][-1])
    utils.report("V2-CD-overprediction",
                 "M=10 CD(alpha~0) exceeds wind tunnel (method behaviour)",
                 res0["CD"] > wt0,
                 f"ours {res0['CD']:.3f} vs WT {wt0:.3f} "
                 f"(D&M LIM digitized 1.704)")


if __name__ == "__main__":
    main()
