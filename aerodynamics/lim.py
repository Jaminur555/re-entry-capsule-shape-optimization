import numpy as np
from .newtonian import cp_newtonian
from .interpolation import cp_cone_interp, cp_prandtl_interp
from ..properties import PART_NOSE, PART_AFTERBODY
from .. import config


def high_hyp_weight(M):
    """Blend weight on the high-hyp (Newtonian) law, continuous form of Eq 3.48-3.49.
       0 at M<=5 (pure cone), 1 at M>=12 (pure Newtonian); smoothstep in between."""
    s = np.clip((M - config.M_HIGH_MIN) / (config.M_LOW_MAX - config.M_HIGH_MIN), 0.0, 1.0)   # (M - 5) / 7
    return 3.0 * s * s - 2.0 * s ** 3


def lim_pressure_coefficient(theta_deg, panel_part, M):
    """Dirkx & Mooij local-inclination pressure coefficient (Table 3.2 / Sec 7.3.1).

    theta_deg  : local surface inclination vs freestream (deg); >0 windward, <0 lee.
    panel_part : per-panel tag (PART_NOSE blunt forebody, PART_AFTERBODY round),
                 broadcastable to theta_deg.
    M          : freestream Mach (scalar).

    Windward nose/shoulder (blunt) -> modified Newtonian
    Windward afterbody (round)     -> cone (low-hyp) cubic-blended to Newtonian (high-hyp)
    Leeward (theta < 0)            -> Prandtl-Meyer expansion
    """
    theta_deg = np.asarray(theta_deg, dtype=float)
    part = np.asarray(panel_part)
    if part.shape != theta_deg.shape:                                  # broadcast (N,) -> (N, P)
        part = part.reshape(part.shape + (1,) * (theta_deg.ndim - part.ndim))

    wind     = theta_deg >= 0.0
    theta_w  = np.where(wind, theta_deg, 0.0)                          # safe windward arg (0 where lee)
    cp_newt  = cp_newtonian(M, theta_w)                                # blunt / high-hyp law
    cp_cone  = cp_cone_interp(M, theta_w)                              # round afterbody, low-hyp law
    w_hi     = high_hyp_weight(M)                                      # 0 at M<=5, 1 at M>=12
    cp_after = (1.0 - w_hi) * cp_cone + w_hi * cp_newt
    cp_lee   = cp_prandtl_interp(M, np.where(wind, 0.0, -theta_deg))   # |theta| on lee

    is_nose  = (part == PART_NOSE)
    is_after = (part == PART_AFTERBODY)
    return np.where(wind & is_nose, cp_newt,
           np.where(wind & is_after, cp_after, cp_lee))

