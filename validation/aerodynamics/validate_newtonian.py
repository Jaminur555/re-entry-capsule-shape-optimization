"""V1a - Modified-Newtonian Cp(theta) vs the digitized Fig 3.5 curves.

Run from the repo root:  python -m validation.aerodynamics.validate_newtonian
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from validation import utils
from capsule_opt.aerodynamics.newtonian import cp_newtonian

BAND = 0.03          # |dCp| pass band (digitization noise dominates)


def curves():
    rows = utils.load_reference("cp_methods_dm2017_digitized.csv")
    out = {}
    for r in rows:
        if r["method"] != "modified_newtonian":
            continue
        key = r["mach"]
        out.setdefault(key, ([], []))
        out[key][0].append(float(r["theta_deg"]))
        out[key][1].append(float(r["cp"]))
    return {k: (np.array(a), np.array(b)) for k, (a, b) in out.items()}


def main():
    data = curves()
    order = sorted(data, key=lambda s: np.inf if s == "inf" else float(s))

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True)
    for ax, key in zip(axes.ravel(), order):
        M = 1.0e4 if key == "inf" else float(key)
        th, cp_ref = data[key]
        cp_ours = np.array([cp_newtonian(M, t) for t in th])
        err = np.abs(cp_ours - cp_ref)
        ok = bool(np.all(err <= BAND))
        utils.report(f"V1a-MN-M{key}", f"Modified Newtonian Cp vs Fig 3.5 (M={key})",
                     ok, f"max|dCp|={err.max():.4f} at theta={th[np.argmax(err)]:.1f} deg")

        thf = np.linspace(0, 90, 181)
        ax.plot(thf, [cp_newtonian(M, t) for t in thf], "-", label="ours")
        ax.plot(th, cp_ref, "o", ms=4, label="digitized D&M")
        ax.set_title(f"M = {key}   max|dCp| = {err.max():.3f}")
        ax.grid(alpha=0.3)
    axes[0, 0].set_ylabel("Cp")
    axes[1, 0].set_ylabel("Cp")
    for ax in axes[1]:
        ax.set_xlabel("local inclination [deg]")
    axes[0, 0].legend()
    fig.suptitle("V1a - Modified Newtonian vs Dirkx (2017) Fig 3.5")
    fig.tight_layout()
    out = utils.FIGURES_DIR / "v1a_newtonian.png"
    fig.savefig(out, dpi=150)
    print(f"  figure -> {out}")


if __name__ == "__main__":
    main()
