"""Exact-Apollo validation case shared by all V-scripts (Dirkx 2017, Sec. 7.2).

Importing FIXES config.Rm_fixed / config.Lc_fixed at the Apollo values
(1.956 / 2.662 m) for the whole process; `restore_defaults()` undoes it.
Sign conventions (alpha flip, Cm scaling, CG z-sign): validation/README.md
section 2.

Run the self-test from the repo root:
    python -m validation.apollo_case
"""

from __future__ import annotations

import numpy as np

from capsule_opt import config
from capsule_opt.properties import get_capsule_properties
from capsule_opt.aerodynamics.aero_database import build_aero_database

# --- Apollo shape / reference data (book p. 135-136) -------------------
APOLLO_RN, APOLLO_RM, APOLLO_RS = 4.694, 1.956, 0.196   # [m]
APOLLO_THETA_C, APOLLO_LC = 33.0, 2.662                  # [deg], [m]
APOLLO_CREF = 3.9116                                     # [m] = 2 * Rm
APOLLO_RCOM = (1.0367, 0.0, 0.1369)                      # [m], D&M axes
# Book prints Sref = 39.441 m2 -- typo, inconsistent with Fig 7.7 (README sec 2)
APOLLO_SREF = np.pi * APOLLO_RM ** 2                     # = 12.0173 m2

# --- normalized design vars (inverse of geometry.cap_params) ------------
APOLLO_RN_NORM = (APOLLO_RN - config.RN_MIN) / (config.RN_MAX - config.RN_MIN)
APOLLO_RS_NORM = (APOLLO_RS - config.RS_MIN) / (config.RS_MAX - config.RS_MIN)
APOLLO_RTHETA_NORM = ((APOLLO_THETA_C - config.THETA_C_MIN)
                      / (config.THETA_C_MAX - config.THETA_C_MIN))

_DEFAULTS = {"Rm_fixed": config.Rm_fixed, "Lc_fixed": config.Lc_fixed}


def apply_apollo_overrides():
    """Fix the shared geometry at the Apollo values (idempotent)."""
    config.Rm_fixed, config.Lc_fixed = APOLLO_RM, APOLLO_LC


def restore_defaults():
    """Undo the Apollo overrides (validation scripts never need this)."""
    config.Rm_fixed, config.Lc_fixed = _DEFAULTS["Rm_fixed"], _DEFAULTS["Lc_fixed"]


apply_apollo_overrides()      # module import == the process runs the Apollo case


def solve_cg_params():
    """Invert the CG placement so (x_cg, z_cg) lands on D&M's rcom.

    Both mappings are linear, so the inversion is direct; z_cg =
    -dz_over_h * h_local, and the Apollo offset needs z_cg = -0.1369.
    """
    base = get_capsule_properties(APOLLO_RN_NORM, APOLLO_RS_NORM,
                                  APOLLO_RTHETA_NORM, cg_params=(0.5, 0.7))
    # default dx_over_L = normalized_param(0.5) = MIN + 0.5*(MAX-MIN)
    dx_over_l_default = config.DX_OVER_L_MIN + 0.5 * (config.DX_OVER_L_MAX - config.DX_OVER_L_MIN)
    x_centroid = base["X_cg"] - dx_over_l_default * base["L_total"]
    dx_over_l = (APOLLO_RCOM[0] - x_centroid) / base["L_total"]
    h_local = float(np.interp(APOLLO_RCOM[0], base["X"], base["Y"]))
    dz_over_h = APOLLO_RCOM[2] / h_local          # magnitude; z_cg comes out -ve

    dx_norm = (dx_over_l - config.DX_OVER_L_MIN) / (config.DX_OVER_L_MAX - config.DX_OVER_L_MIN)
    dz_norm = (dz_over_h - config.DZ_OVER_H_MIN) / (config.DZ_OVER_H_MAX - config.DZ_OVER_H_MIN)
    if not (0.0 <= dx_norm <= 1.0 and 0.0 <= dz_norm <= 1.0):
        raise RuntimeError(
            f"Apollo CG outside parameterized bounds: dx_norm={dx_norm:.3f}, "
            f"dz_norm={dz_norm:.3f} (bounds are [0, 1] by construction)")
    return dx_norm, dz_norm


APOLLO_CG_PARAMS = solve_cg_params()

# --- D&M validation-entry initial conditions (book p. 142) --------------
ENTRY_ICS = dict(ho=120000.0, Vo=7830.0, gamma0_deg=-4.0,
                 lat0_deg=-23.75, lon0_deg=225.5, chi0_deg=49.6,
                 mass_kg=4532.0)

_cache: dict = {}


def apollo_props():
    """Property dict of the exact-Apollo capsule with the rcom CG."""
    if "props" not in _cache:
        _cache["props"] = get_capsule_properties(
            APOLLO_RN_NORM, APOLLO_RS_NORM, APOLLO_RTHETA_NORM,
            cg_params=APOLLO_CG_PARAMS)
    return _cache["props"]


def cm_ref_scale():
    """Factor converting our Cm to the D&M normalization (cref = 3.9116 m)."""
    if "cm_scale" not in _cache:
        _cache["cm_scale"] = apollo_props()["L_total"] / APOLLO_CREF
    return _cache["cm_scale"]


def apollo_database(rebuild=False):
    """Cached (interp_CD, interp_CL, interp_CM, alpha_trim) for the case."""
    if rebuild or "db" not in _cache:
        _cache["db"] = build_aero_database(
            APOLLO_RN_NORM, APOLLO_RS_NORM, APOLLO_RTHETA_NORM,
            cg_params=APOLLO_CG_PARAMS)
    return _cache["db"]


def run_entry(dt_max=2.0, t_max=7200.0):
    """Propagate the D&M validation entry (Fig 7.17 conditions); cached.

    gamma0 = -4 deg (the validation entry; the -2 deg default is the
    optimization case of book Sect. 7.3.2).
    """
    if "traj" not in _cache:
        from capsule_opt.trajectory.trajectory import propagate_trajectory
        iCD, iCL, _iCM, atr = apollo_database()
        _cache["traj"] = propagate_trajectory(
            APOLLO_RN_NORM, APOLLO_RS_NORM, APOLLO_RTHETA_NORM,
            ENTRY_ICS["mass_kg"],
            ho=ENTRY_ICS["ho"], Vo=ENTRY_ICS["Vo"],
            gamma0_deg=ENTRY_ICS["gamma0_deg"],
            lat0_deg=ENTRY_ICS["lat0_deg"], lon0_deg=ENTRY_ICS["lon0_deg"],
            chi0_deg=ENTRY_ICS["chi0_deg"],
            cg_params=APOLLO_CG_PARAMS,
            t_max=t_max, dt_max=dt_max,
            shape_props=apollo_props(),
            interp_CD=iCD, interp_CL=iCL, alpha_trim_arr=atr)
    return _cache["traj"]


if __name__ == "__main__":
    props = apollo_props()

    # geometry / CG fidelity
    assert abs(props["R_max"] - APOLLO_RM) < 1e-4, props["R_max"]
    assert abs(props["A_ref"] - APOLLO_SREF) < 0.05, props["A_ref"]
    assert abs(props["X_cg"] - APOLLO_RCOM[0]) < 1e-9, props["X_cg"]
    assert abs(props["Z_cg"] + APOLLO_RCOM[2]) < 1e-9, props["Z_cg"]
    print(f"A_ref = {props['A_ref']:.4f} m2 (target {APOLLO_SREF:.4f})")
    print(f"L_total = {props['L_total']:.4f} m, volume = {props['volume']:.2f} m3")
    print(f"CG = ({props['X_cg']:.4f}, {props['Z_cg']:+.4f}) m "
          f"[D&M rcom = (1.0367, +0.1369), our z sign flipped]")

    # trim behaviour vs the digitized Fig 7.7 reference
    iCD, iCL, iCM, atr = apollo_database()
    idx10 = int(np.argmin(np.abs(config.mach_nodes - 10.0)))
    idx3 = int(np.argmin(np.abs(config.mach_nodes - 3.0)))
    for label, k, lo, hi, dm in (("M=10", idx10, 19.0, 26.0, -22.2),
                                 ("M=3", idx3, 21.0, 28.0, -24.8)):
        assert lo < atr[k] < hi, f"trim {label}: {atr[k]:.2f} outside ({lo},{hi})"
        print(f"trim {label}: {atr[k]:.2f} deg (digitized D&M LIM: {dm} deg)")

    a10 = float(atr[idx10])
    ld = float(iCL([[10.0, a10]])[0]) / float(iCD([[10.0, a10]])[0])
    assert 0.28 < ld < 0.38, ld
    print(f"L/D at trim (M=10) = {ld:.3f} (digitized D&M LIM ~ 0.335)")
    print("\napollo_case self-test PASS")
