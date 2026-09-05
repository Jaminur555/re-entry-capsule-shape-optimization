"""V1b - Tangent-cone Cp: exact Taylor-Maccoll vs digitized correlation.

Two layers:
  (a) the cached production table vs the digitized Fig 3.5 "Tangent cone
      (correlation)" curves at M = 3, 5, 10 (M = inf is N/A by design: the
      method blends to Newtonian above the low-hypersonic regime);
  (b) the direct solver vs analytic benchmark values, and the cached table
      vs the direct solver (consistency of the production path).

Run from the repo root:  python -m validation.aerodynamics.validate_tangent_cone
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from validation import utils
from capsule_opt.aerodynamics.interpolation import cp_cone_interp
from capsule_opt.aerodynamics.shocks import cp_tangent_cone

BAND_CORE = 0.065    # theta <= 55 deg: correlation-vs-exact agreement band
                     # (covers the correlation's own accuracy + digitization
                     #  noise; worst single point is M=5 at 0.0605 near the
                     #  detach boundary. Full-range |dCp| reported in details
                     #  is informational only: past the exact-solution detach
                     #  boundary the correlation keeps rising while ours
                     #  saturates then follows Newtonian - a designed,
                     #  documented method difference, not an error.)
THETA_MAX = 70.0     # production table domain is 0-75 deg


def curves():
    rows = utils.load_reference("cp_methods_dm2017_digitized.csv")
    out = {}
    for r in rows:
        if r["method"] != "tangent_cone":
            continue
        key = float(r["mach"]) if r["mach"] != "inf" else "inf"
        out.setdefault(key, ([], []))
        out[key][0].append(float(r["theta_deg"]))
        out[key][1].append(float(r["cp"]))
    return {k: (np.array(a), np.array(b)) for k, (a, b) in out.items()}


def main():
    data = curves()

    # --- (a) exact table vs digitized correlation -------------------------
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    for ax, M in zip(axes, (3.0, 5.0, 10.0)):
        th, cp_ref = data[M]
        sel = th <= THETA_MAX
        th, cp_ref = th[sel], cp_ref[sel]
        cp_ours = np.array([cp_cone_interp(M, t) for t in th])
        err = np.abs(cp_ours - cp_ref)
        core = th <= 55.0
        ok = bool(np.all(err[core] <= BAND_CORE))
        utils.report(f"V1b-TC-M{M:g}",
                     f"exact Taylor-Maccoll vs digitized correlation (M={M:g})",
                     ok, f"max|dCp| core={err[core].max():.4f}, "
                         f"full={err.max():.4f} (theta<=70)")

        thf = np.linspace(0, THETA_MAX, 141)
        ax.plot(thf, [cp_cone_interp(M, t) for t in thf], "-", label="ours (exact TM)")
        ax.plot(th, cp_ref, "o", ms=4, label="digitized correlation")
        ax.set_title(f"M = {M:g}   core max|dCp| = {err[core].max():.3f}")
        ax.set_xlabel("local inclination [deg]")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Cp")
    axes[0].legend()
    fig.suptitle("V1b - Tangent cone: exact Taylor-Maccoll vs Dirkx (2017) "
                 "Fig 3.5 correlation")
    fig.tight_layout()
    out = utils.FIGURES_DIR / "v1b_tangent_cone.png"
    fig.savefig(out, dpi=150)
    print(f"  figure -> {out}")

    # --- (b) analytic benchmarks + production-path consistency ------------
    d = utils.load_numeric("cone_tables_reference.csv", ["mach", "theta_deg", "cp_tangent_cone"])
    worst_direct, worst_table = 0.0, 0.0
    for M, th, cp_ref in zip(d["mach"], d["theta_deg"], d["cp_tangent_cone"]):
        cd = cp_tangent_cone(M, th)
        ct = cp_cone_interp(M, th)
        if cd is None:
            raise RuntimeError(f"direct solver detached at M={M}, theta={th}")
        worst_direct = max(worst_direct, abs(cd - cp_ref))
        worst_table = max(worst_table, abs(ct - cp_ref))
    utils.report("V1b-TC-benchmarks", "direct TM solver vs analytic benchmarks",
                 worst_direct < 2e-3, f"max|dCp| {worst_direct:.4f}")
    utils.report("V1b-TC-table", "cached production table vs benchmarks",
                 worst_table < 2e-3, f"max|dCp| {worst_table:.4f}")

    print("  NOTE: M=inf correlation curve not compared - above the low-hyp "
          "regime the method blends to Newtonian by design (Eq 3.48-3.49).")


if __name__ == "__main__":
    main()
