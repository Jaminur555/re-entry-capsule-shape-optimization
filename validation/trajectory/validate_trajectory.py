"""V4b - Entry-trajectory profiles vs digitized Fig 7.17 (D&M LIM source).

Run from the repo root:  python -m validation.trajectory.validate_trajectory
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import validation.apollo_case as apollo
from validation import utils
from capsule_opt import config
from capsule_opt.atmosphere import atmosphere
from capsule_opt.optimization.metrics import heat_load, load_factor
from capsule_opt.trajectory.trajectory import bank_angle_sigma

BANDS = {              # absolute bands vs the digitized D&M LIM curves
    "altitude_km": 1.5,       # km
    "lift_to_drag": 0.02,
    "bank_angle_deg": 6.0,    # deg
    "load_factor_g": 0.25,    # g
}
HEAT_RATE_FRAC = 0.07         # heat-rate band as a fraction of the ref peak


def digitized(profile, source="lim", indep="mach"):
    rows = utils.load_reference("apollo_entry_reference_dm2017.csv")
    m, v = [], []
    for r in rows:
        if r["profile"] == profile and r["aero_source"] == source \
                and r["indep_name"] == indep:
            m.append(float(r["indep_value"]))
            v.append(float(r["value"]))
    o = np.argsort(m)
    return np.array(m)[o], np.array(v)[o]


def main():
    traj = apollo.run_entry()
    iCD, iCL, _iCM, _atr = apollo.apollo_database()
    props = apollo.apollo_props()
    m_kg = apollo.ENTRY_ICS["mass_kg"]
    A_ref = props["A_ref"]

    M = np.asarray(traj["Mach"], float)
    o = np.argsort(M)                       # ascending abscissa for np.interp
    Ms = M[o]

    # per-point CD/CL at the trim AoA (same look-ups as entry_eom)
    at = np.clip(np.asarray(traj["alpha_trim"], float),
                 config.aoa_nodes[0], config.aoa_nodes[-1])
    Mc = np.clip(M, config.mach_nodes[0], config.mach_nodes[-1])
    CD = np.array([float(iCD(np.array([[Mc[k], at[k]]]))[0]) for k in range(len(M))])
    CL = np.array([float(iCL(np.array([[Mc[k], at[k]]]))[0]) for k in range(len(M))])

    # bank angle: propagated state (rate-limited no-skip law) when available,
    # otherwise the raw post-hoc command (same law as entry_eom)
    if "bank_deg" in traj:
        sig = np.deg2rad(np.asarray(traj["bank_deg"], float))
    else:
        h = np.maximum(np.asarray(traj["h"], float), 0.0)
        V = np.asarray(traj["V"], float)
        rho = np.array([atmosphere(z)["rho"] for z in h])
        L = 0.5 * rho * V ** 2 * A_ref * np.abs(CL)
        sig = np.array([bank_angle_sigma(
            m_kg, L[k], V[k], h[k], np.deg2rad(traj["gamma_deg"][k]),
            np.deg2rad(traj["lat_deg"][k]), np.deg2rad(traj["chi_deg"][k]))
            for k in range(len(M))])

    heat = heat_load(traj, apollo.APOLLO_RN_NORM, apollo.APOLLO_RS_NORM,
                     apollo.APOLLO_RTHETA_NORM)
    nres = load_factor(traj, apollo.APOLLO_RN_NORM, apollo.APOLLO_RS_NORM,
                       apollo.APOLLO_RTHETA_NORM, m_kg,
                       shape_props=props, interp_CD=iCD)

    ours = {
        "altitude_km": np.asarray(traj["h"], float) / 1e3,
        "lift_to_drag": CL / CD,
        "bank_angle_deg": np.abs(np.rad2deg(sig)),
        "heat_rate_w_m2": heat["q_stag_arr"],
        "load_factor_g": nres["n_arr"],
    }

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for ax, name in zip(axes.ravel(), list(BANDS) + ["heat_rate_w_m2"]):
        # altitude is digitized vs time; every other profile vs Mach
        if name == "altitude_km":
            x_ref, v_ref = digitized(name, indep="time_s")
            x_ours = np.asarray(traj["t"], float)
            x_lab = "time [s]"
            y_ours = ours[name]
            sel = x_ref <= x_ours.max()
        else:
            x_ref, v_ref = digitized(name)
            x_lab = "Mach"
            x_ours, y_ours = Ms, ours[name][o]
            sel = (x_ref >= 3.0) & (x_ref <= Ms.max())
        if name == "lift_to_drag" and np.mean(v_ref) < 0:
            v_ref = -v_ref               # D&M alpha < 0 sign convention
        x_ref, v_ref = x_ref[sel], v_ref[sel]
        ours_ip = np.interp(x_ref, x_ours, y_ours)

        if name == "heat_rate_w_m2":
            band = HEAT_RATE_FRAC * v_ref.max()
        else:
            band = BANDS[name]
        err = np.abs(ours_ip - v_ref)
        utils.report(f"V4b-{name}", f"{name} vs digitized D&M LIM (Fig 7.17)",
                     bool(np.all(err <= band)),
                     f"max err {err.max():.3g} at {x_lab.split()[0]}="
                     f"{x_ref[np.argmax(err)]:.1f}, band {band:.3g}")

        ax.plot(x_ours, y_ours, label="ours")
        ax.plot(x_ref, v_ref, ".", ms=3, label="D&M LIM (digitized)")
        ax.set_xlabel(x_lab)
        ax.set_ylabel(name)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    axes[0, 0].set_title("V4b - entry profiles vs Fig 7.17 (LIM)")
    axes[1, 2].axis("off")
    fig.tight_layout()
    out = utils.FIGURES_DIR / "v4b_trajectory.png"
    fig.savefig(out, dpi=150)
    print(f"  figure -> {out}")


if __name__ == "__main__":
    main()
