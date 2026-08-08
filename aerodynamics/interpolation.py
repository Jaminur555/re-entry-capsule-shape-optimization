from .. import config
import numpy as np
from .cp_tables import load_cone_table, build_prandtl_table

# ===========================================================================
#          Table interpolation helper (bilinear in M, linear in theta)
# ===========================================================================

def interp_cp(table, machs, thetas_deg, M, theta_deg):
    """Interpolate 'table' (shape (nM, nTheta)) at a single Mach M and an array
       of panel deflection angles theta_deg (any shape)."""
    machs      = np.asarray(machs)
    thetas_deg = np.asarray(thetas_deg)
    theta_deg  = np.asarray(theta_deg)

    M_c     = np.clip(M, machs[0], machs[-1])
    iM      = np.clip(np.searchsorted(machs, M_c) - 1, 0, machs.size - 2)
    wM      = (M_c - machs[iM]) / (machs[iM + 1] - machs[iM])
    row     = table[iM, :] + wM * (table[iM + 1] - table[iM, :])
    theta_c = np.clip(theta_deg, thetas_deg[0], thetas_deg[-1])
    return np.interp(theta_c, thetas_deg, row)


# ===========================================================================
# Local-inclination law databases (Dirkx & Mooij Sec 3.3.2, 7.3.1)
#   windward nose/shoulder (blunt)  -> modified Newtonian   (analytic)
#   windward afterbody (round)      -> tangent cone (low-hyp) blended to
#                                      modified Newtonian (high-hyp)
#   leeward / lee (theta < 0)       -> Prandtl-Meyer        (analytic/table)
# Only the cone law is expensive, so it alone is tabulated and cached.
# ===========================================================================

M_LOW_MAX, M_HIGH_MIN   = config.M_LOW_MAX, config.M_HIGH_MIN
CONE_MACHS, CONE_THETAS = config.CONE_MACHS, config.CONE_THETAS
CONE_CACHE = config.CONE_CACHE

def cp_cone_interp(M, theta_deg):
    """Windward afterbody Cp from the cached cone table (bilinear in M, theta)."""
    return interp_cp(CP_CONE, CONE_MACHS, CONE_THETAS, M, theta_deg)

# Eager build/load at import (one-time build, then instant CSV load).
CP_CONE = load_cone_table()


# ==============================================================================
# Prandtl-Meyer lee-side table (cheap, built eagerly) + LIM pressure dispatcher
# ==============================================================================

PRANDTL_MACHS, PRANDTL_THETAS = config.PRANDTL_MACHS, config.PRANDTL_THETAS

def cp_prandtl_interp(M, theta_turn_deg):
    """Lee-side Cp from the Prandtl-Meyer table (theta_turn = |theta| >= 0, deg)."""
    return interp_cp(CP_PRANDTL, PRANDTL_MACHS, PRANDTL_THETAS, M, theta_turn_deg)
CP_PRANDTL = build_prandtl_table()
