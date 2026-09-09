"""V4c sensitivity - entry runs isolating the remaining-gap suspects.

BASE (shipped law) vs D&M: Qs +32%, Sg +16%, n_tot peak +22%, q-dot peak
0.7% (PASS). Variants: (a) low-M L/D bias (+0.023) - rescale CL by
`cl_frac` below Mach `cl_cut`; (b) bank modulation starts late (sigma = 0
until M~24 vs D&M onset ~25.5) - GATE variant activates modulation while
cos(sigma) > -eps with command arccos(clip(cos_s, 0, 1)); eps=0 reproduces
the shipped law exactly.

Run from the repo root:  python -m validation.system.sensitivity_bank_cl
"""

import numpy as np

import validation.apollo_case as apollo
from capsule_opt import config
from capsule_opt.trajectory import trajectory as tt
from capsule_opt.trajectory.trajectory import propagate_trajectory
from capsule_opt.optimization.metrics import heat_load, load_factor

# D&M reference values (Table 7.1 LIM column + digitized Fig 7.17 peaks)
QS_REF, SG_REF = 53.06e6, 1430.2e3     # [J m^-2], [m]
QDOT_REF, NTOT_REF = 0.616e6, 4.36     # [W m^-2], [g]

ORIGINAL_BANK = tt.bank_angle_sigma    # shipped law, restore after runs


def make_scaled_CL(iCL, frac, cut=8.0):
    """Wrap interp_CL: multiply CL by `frac` below Mach `cut`."""
    def wrapped(pts):
        pts = np.atleast_2d(np.asarray(pts, float))
        out = np.array(iCL(pts), float)
        out[pts[:, 0] < cut] *= frac
        return out
    return wrapped


def make_gated_bank(eps, cap_deg=None, sin_gamma=False):
    """Shipped no-skip law with an activation gate at cos(sigma) > -eps,
    optionally capping the commanded bank angle (D&M plateau ~82 deg).
    sin_gamma=True uses the Eq-2.40-as-printed omega term (sin gamma)."""
    def bank(m, L_force, V, h, gamma, lat, chi):
        r = config.r_earth + h
        g = config.go * (config.r_earth / r) ** 2
        om = np.sin(gamma) if sin_gamma else np.sin(chi)
        cos_s = (m / max(L_force, 1e-6)) * (
            g * (1 - V ** 2 / (g * r)) * np.cos(gamma)
            - 2 * config.omega_earth * V * np.cos(lat) * om)
        if cos_s >= 1.0 or cos_s <= -eps:
            return 0.0
        sig = float(np.arccos(np.clip(cos_s, 0.0, 1.0)))
        if cap_deg is not None:
            sig = min(sig, np.deg2rad(cap_deg))
        return sig
    return bank


def run_variant(cl_frac=1.0, cl_cut=8.0, gate_eps=0.0, cap_deg=None,
                sin_gamma=False):
    """One entry run; returns (traj, metrics dict). Monkey-patches restored."""
    iCD, iCL, _iCM, atr = apollo.apollo_database()
    if cl_frac != 1.0:
        iCL = make_scaled_CL(iCL, cl_frac, cl_cut)
    tt.bank_angle_sigma = (make_gated_bank(gate_eps, cap_deg, sin_gamma)
                           if gate_eps > 0.0 else ORIGINAL_BANK)
    try:
        traj = propagate_trajectory(
            apollo.APOLLO_RN_NORM, apollo.APOLLO_RS_NORM,
            apollo.APOLLO_RTHETA_NORM, apollo.ENTRY_ICS["mass_kg"],
            ho=apollo.ENTRY_ICS["ho"], Vo=apollo.ENTRY_ICS["Vo"],
            gamma0_deg=apollo.ENTRY_ICS["gamma0_deg"],
            lat0_deg=apollo.ENTRY_ICS["lat0_deg"],
            lon0_deg=apollo.ENTRY_ICS["lon0_deg"],
            chi0_deg=apollo.ENTRY_ICS["chi0_deg"],
            cg_params=apollo.APOLLO_CG_PARAMS, t_max=7200.0, dt_max=2.0,
            shape_props=apollo.apollo_props(),
            interp_CD=iCD, interp_CL=iCL, alpha_trim_arr=atr)
    finally:
        tt.bank_angle_sigma = ORIGINAL_BANK

    heat = heat_load(traj, apollo.APOLLO_RN_NORM, apollo.APOLLO_RS_NORM,
                     apollo.APOLLO_RTHETA_NORM)
    nres = load_factor(traj, apollo.APOLLO_RN_NORM, apollo.APOLLO_RS_NORM,
                       apollo.APOLLO_RTHETA_NORM, apollo.ENTRY_ICS["mass_kg"],
                       shape_props=apollo.apollo_props(), interp_CD=iCD)

    # D&M Eq 6.13 total-force peak (like V4c)
    Mc = np.clip(traj["Mach"], config.mach_nodes[0], config.mach_nodes[-1])
    Ac = np.clip(traj["alpha_trim"], config.aoa_nodes[0], config.aoa_nodes[-1])
    CDa = np.array([float(iCD(np.array([[Mc[k], Ac[k]]]))[0])
                    for k in range(Mc.size)])
    Cla = np.array([float(iCL(np.array([[Mc[k], Ac[k]]]))[0])
                    for k in range(Mc.size)])
    n_tot = nres["n_arr"] * np.sqrt(1.0 + (Cla / CDa) ** 2)

    bank = traj["bank_deg"]
    on = np.where(bank > 1.0)[0]
    met = {
        "Qs_MJ": heat["Qs"] / 1e6,
        "Sg_km": traj["sg"] / 1e3,
        "qdot_pk_MW": heat["q_stag_max"] / 1e6,
        "n_pk": nres["n_max"],
        "ntot_pk": float(n_tot.max()),
        "ntot_at_M": float(Mc[int(np.argmax(n_tot))]),
        "bank_onset_M": float(traj["Mach"][on[0]]) if on.size else np.nan,
        "bank_max_deg": float(bank.max()),
    }
    return traj, met


def main():
    variants = [
        ("BASE", {}),
        ("GATE.15cap82", dict(gate_eps=0.15, cap_deg=82.0)),
        ("GATE.30cap82", dict(gate_eps=0.30, cap_deg=82.0)),
        ("SING cap82", dict(gate_eps=1e-3, cap_deg=82.0, sin_gamma=True)),
    ]
    print(f"{'variant':<15} {'Qs MJ':>7} {'Sg km':>7} {'qdot MW':>8} "
          f"{'n_pk':>6} {'n_tot':>6} {'@M':>5} {'onset M':>8} {'sig max':>8}")
    print(f"{'D&M':<15} {QS_REF/1e6:>7.2f} {SG_REF/1e3:>7.1f} "
          f"{QDOT_REF/1e6:>8.3f} {'':>6} {NTOT_REF:>6.2f}")
    for label, kw in variants:
        _traj, m = run_variant(**kw)
        print(f"{label:<15} {m['Qs_MJ']:>7.2f} {m['Sg_km']:>7.1f} "
              f"{m['qdot_pk_MW']:>8.3f} {m['n_pk']:>6.2f} {m['ntot_pk']:>6.2f} "
              f"{m['ntot_at_M']:>5.1f} {m['bank_onset_M']:>8.1f} "
              f"{m['bank_max_deg']:>8.1f}")

    # DIAG: BASE n/q-dot/h profile around D&M's n-peak point (M=12.8) --
    # does our n curve *cross* their peak value (shape/timing issue) or sit
    # uniformly above it (structural scale)?
    traj, _m = run_variant()
    iCD, iCL, _iCM, _atr = apollo.apollo_database()
    Mc = np.clip(traj["Mach"], config.mach_nodes[0], config.mach_nodes[-1])
    Ac = np.clip(traj["alpha_trim"], config.aoa_nodes[0], config.aoa_nodes[-1])
    CDa = np.array([float(iCD(np.array([[Mc[k], Ac[k]]]))[0])
                    for k in range(Mc.size)])
    Cla = np.array([float(iCL(np.array([[Mc[k], Ac[k]]]))[0])
                    for k in range(Mc.size)])
    rho = 2 * traj["q_dyn"] / traj["V"] ** 2
    n_tot = (traj["q_dyn"] * apollo.apollo_props()["A_ref"] * CDa
             * np.sqrt(1 + (Cla / CDa) ** 2)
             / (apollo.ENTRY_ICS["mass_kg"] * config.go))
    # Sutton-Graves, same k as metrics.q_stag; Reff(4.694, 0.196 | Rm=1.956)
    # = 3.8988 m
    qdot = 1.83e-4 * np.sqrt(rho / 3.8988) * traj["V"] ** 3
    print("\nDIAG BASE: M / h_km / n_tot_g / qdot_MW / bank_deg "
          "(D&M npk 4.36 g at M=12.8)")
    for Mt in (20.0, 17.9, 15.0, 12.8, 10.0, 8.0, 5.0):
        k = int(np.argmin(np.abs(traj["Mach"] - Mt)))
        print(f"  M={traj['Mach'][k]:5.1f}  h={traj['h'][k]/1e3:6.1f}  "
              f"n={n_tot[k]:5.2f}  q={qdot[k]/1e6:5.3f}  "
              f"sig={traj['bank_deg'][k]:5.1f}")
    kq = int(np.argmax(qdot))
    print(f"  qdot peak at M={traj['Mach'][kq]:.1f}, h={traj['h'][kq]/1e3:.1f} km, "
          f"n={n_tot[kq]:.2f} g;  n(D-only) there="
          f"{(traj['q_dyn'][kq]*apollo.apollo_props()['A_ref']*CDa[kq]/(apollo.ENTRY_ICS['mass_kg']*config.go)):.2f}")


if __name__ == "__main__":
    main()
