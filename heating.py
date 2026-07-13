"""Convective heating: effective nose radius and stagnation / shoulder heat flux."""

import numpy as np

from . import config
from .parameterization import cap_params
from .properties import get_capsule_properties


#=======================================================================================
#                 Digitalized Zoby & Sullivan chart (effective nose radius)
#=======================================================================================
# Axes: x = Rm/Rn; y = Rm/Reff; family parameter = Rs/Rm

ZS_Rm_over_Rn = np.linspace(0, 1, 11)

ZS_Table = {           # Rs/Rm: Rm/Reff at each Rm/Rn
    0.00: np.array([0.286, 0.310, 0.355, 0.411, 0.478, 0.556, 0.639, 0.728, 0.824, 0.916, 1.0]),
    0.10: np.array([0.297, 0.327, 0.372, 0.427, 0.490, 0.560, 0.639, 0.728, 0.824, 0.916, 1.0]),
    0.20: np.array([0.318, 0.339, 0.378, 0.427, 0.490, 0.560, 0.639, 0.728, 0.824, 0.916, 1.0]),
    0.25: np.array([0.351, 0.362, 0.394, 0.441, 0.500, 0.568, 0.639, 0.728, 0.824, 0.916, 1.0]),
    0.30: np.array([0.368, 0.376, 0.406, 0.445, 0.500, 0.568, 0.639, 0.728, 0.824, 0.916, 1.0])
}

ZS_Rs_over_Rm_keys = np.array(sorted(ZS_Table.keys()))


def effective_nose_radius(Rn, Rs, Rm=config.Rm_fixed):
    """Effective nose radius "Reff" from the Zoby & Sullivan chart.

    Parameters
    ----------
    Rn, Rs, Rm : float
        Nose, shoulder and maximum body radii [m].

    Returns
    -------
    float
        Effective nose radius "Reff" [m] (falls back to "Rn" if undefined).
    """
    if Rm <= 0 or Rn <= 0 or Rs < 0:
        raise ValueError("Rn, Rs and Rm must be non-negative")
    if Rm > Rn:
        raise ValueError("Constraint Rm < Rn violated")

    Rm_over_Rn = np.clip(Rm / Rn, 0.0, 1.0)
    Rs_over_Rm = np.clip(Rs / Rm, 0.0, 0.3)  # chart only goes to 0.30

    # for each Rs/Rm key, interpolate in the Rm/Rn direction
    Rm_over_Reff_at_keys = np.array([np.interp(Rm_over_Rn, ZS_Rm_over_Rn, ZS_Table[k])
                                     for k in ZS_Rs_over_Rm_keys])

    # interpolate across the Rs/Rm family
    Rm_over_Reff = np.interp(Rs_over_Rm, ZS_Rs_over_Rm_keys, Rm_over_Reff_at_keys)

    if Rm_over_Reff < config.eps:
        return Rn
    Reff = Rm / Rm_over_Reff
    return Reff


def stagnation_heat_flux(rho, V, Rn, Rs=None, Rm=config.Rm_fixed, Tw_by_Taw=0, use_effective_rn=True):
    """Stagnation-point convective heat flux.

    Parameters
    ----------
    rho : float
        Freestream density [kg m^-3].
    V : float
        Freestream velocity [m s^-1].
    Rn : float
        Nose radius [m].
    Rs : float, optional
        Shoulder radius [m] (used only when "use_effective_rn" is True).
    Rm : float, optional
        Maximum body radius [m].
    Tw_by_Taw : float, optional
        Wall-to-adiabatic-wall temperature ratio in [0, 1].
    use_effective_rn : bool, optional
        Use the Zoby-Sullivan effective nose radius instead of "Rn".

    Returns
    -------
    float
        Stagnation heat flux "q_stag" [W m^-2].
    """
    if rho < 0:
        raise ValueError("Density must be non-negative.")
    if V < 0:
        raise ValueError("Velocity must be non-negative.")
    if Rn <= 0:
        raise ValueError("Nose Radius must be positive")
    if not (0.0 <= Tw_by_Taw <= 1):
        raise ValueError("Tw_by_Taw must be in [0,1].")

    if use_effective_rn and Rs is not None:
        R_for_k = effective_nose_radius(Rn, Rs, Rm)
    else:
        R_for_k = Rn

    k = (1.83e-4 / np.sqrt(R_for_k)) * (1.0 - Tw_by_Taw)
    q_stag = k * (rho ** .5) * (V ** 3)

    return q_stag


def shoulder_heat_flux(q_stag, M, AoA_deg, theta_sp_rad, Rs, Rm=config.Rm_fixed):
    """Shoulder heat flux as a correlated fraction of the stagnation value.

    Uses an empirical linear correlation in Mach, AoA, "Rs/Rm" and the
    nose-sphere half-angle "theta_sp_rad".
    """
    c1, c2, c3, c4, c5 = -0.0006, 0.0185, -0.5321, -0.2939, 1.3630

    rs_over_rm = Rs / Rm

    ratio = (c1 * M + c2 * AoA_deg + c3 * rs_over_rm + c4 * theta_sp_rad + c5)
    ratio = max(ratio, 0)

    return ratio * q_stag


def capsule_heating(rho, V, M, AoA_deg, rn, rs, r_theta, Tw_by_Taw):
    """Stagnation and shoulder heat flux for a given capsule shape and state.

    Parameters
    ----------
    rho, V, M, AoA_deg : float
        Freestream density [kg m^-3], velocity [m s^-1], Mach and angle of
        attack [deg].
    rn, rs, r_theta : float
        Normalized design parameters in [0, 1].
    Tw_by_Taw : float
        Wall-to-adiabatic-wall temperature ratio in [0, 1].

    Returns
    -------
    dict
        Stagnation flux, shoulder flux, effective nose radius, shoulder/stag
        ratio and the physical "Rn", "Rs".
    """
    props = get_capsule_properties(rn, rs, r_theta)
    Rn, Rs, theta_c_deg = cap_params(rn, rs, r_theta)

    theta_sp1 = props['theta_sp1']  # nose-sphere half-angle

    R_eff = effective_nose_radius(Rn, Rs)
    q_stag = stagnation_heat_flux(rho, V, Rn=Rn, Rs=Rs, Tw_by_Taw=Tw_by_Taw, use_effective_rn=True)
    q_shldr = shoulder_heat_flux(q_stag, M, AoA_deg, theta_sp1, Rs)

    return {
        'q_stag': q_stag,
        'q_shoulder': q_shldr,
        'R_eff': R_eff,
        'ratio' : q_shldr / q_stag if q_stag > 0 else 0.0,
        'Rn': Rn,
        'Rs': Rs
    }
