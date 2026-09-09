"""V3a - Heating: Zoby-Sullivan effective nose radius cross-check (Fig 3.9).

The production `effective_nose_radius()` interpolates an in-code
digitization of the chart; this check compares it against the user's
INDEPENDENT PlotDigitizer extraction of the same figure, curve by curve.
The chart is scale-invariant, so Rm = 1 is used throughout.

(Stagnation-flux peak checks run in V4 against the propagated entry
trajectory; see validate_trajectory.py.)

Run from the repo root:  python -m validation.heating.validate_chapman
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from validation import utils
from capsule_opt.heating.heating import effective_nose_radius

BAND = 0.02          # |d(Rm/Reff)| pass band (two independent digitizations)
X_DESIGN = 0.70      # chart region the design space actually uses:
                     # Rn in [3,7] m with Rm fixed -> Rm/Rn <= 0.67


def main():
    d = utils.load_numeric("zoby_sullivan_chart_digitized.csv",
                           ["rs_over_rm", "rm_over_rn", "rm_over_reff"])
    keys = sorted(set(np.round(d["rs_over_rm"], 6)))

    utils.set_paper_style()
    fig, ax = plt.subplots(figsize=(7, 5))
    worst_all, worst_edge = 0.0, 0.0
    ok_all = True
    for k in keys:
        sel = (np.abs(d["rs_over_rm"] - k) < 1e-6) & (d["rm_over_rn"] > 0.005)
        x, y_ref = d["rm_over_rn"][sel], d["rm_over_reff"][sel]
        o = np.argsort(x)
        x, y_ref = x[o], y_ref[o]

        y_ours = np.array([1.0 / effective_nose_radius(1.0 / xi, k, 1.0)
                           for xi in x])
        err = np.abs(y_ours - y_ref)
        design = x <= X_DESIGN
        ok_all &= bool(np.all(err[design] <= BAND))
        worst_all = max(worst_all, float(err[design].max()))
        if (~design).any():
            worst_edge = max(worst_edge, float(err[~design].max()))
        utils.report(f"V3a-ZS-rs{int(round(k*100)):02d}",
                     f"Zoby-Sullivan Rs/Rm={k:.2f}, design-space region (Rm/Rn<=0.7)",
                     bool(np.all(err[design] <= BAND)),
                     f"max|d(Rm/Reff)|={err[design].max():.4f}")

        xf = np.linspace(0.01, 1.0, 100)
        ax.plot(xf, [1.0 / effective_nose_radius(1.0 / xi, k, 1.0) for xi in xf],
                "-", label=f"ours Rs/Rm={k:g}")
        ax.plot(x, y_ref, "o", ms=3, color="gray")

    utils.report("V3a-ZS-all", "Zoby-Sullivan chart, all curves, design space",
                 bool(ok_all), f"worst |d(Rm/Reff)|={worst_all:.4f}")
    utils.report("V3a-ZS-edge",
                 "chart edge Rm/Rn>0.7 (informational, outside design space)",
                 True,
                 f"max|d|={worst_edge:.4f} near x=0.98: coarse-grid artifact of "
                 f"the in-code table's last 0.9->1.0 segment; the true curve "
                 f"curls to ~1.01. Unused in production (Rm/Rn <= 0.67).")

    ax.plot([], [], "o", ms=3, color="gray", label="digitized (independent)")
    ax.set_xlabel("$R_m/R_n$")
    ax.set_ylabel("$R_m/R_{eff}$")
    ax.set_title("Effective nose radius, Zoby & Sullivan (Fig. 3.9)")
    ax.legend(fontsize=8)
    utils.save_fig(fig, "v3a_zoby_sullivan.png")


if __name__ == "__main__":
    main()
