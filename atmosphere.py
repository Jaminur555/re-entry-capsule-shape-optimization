"""Atmospheric model: US Standard Atmosphere 1976 (0-120 km).

Below 86 km the original US-76 geopotential-layer model is used.

Above 86 km (the homopause) the atmosphere is no longer compositionally mixed:
the mean molecular weight falls with altitude, so density cannot be recovered
from P and T with the constant sea-level gas constant. The 86-120 km region is
therefore interpolated directly from the tabulated T, P, rho in
"config.H_EXT / T_EXT / P_EXT / RHO_EXT" (USA1976 Tabulated data).
The table's first row reproduces the US-76 86 km boundary, so the two pieces
join continuously.

Reference: NOAA/NASA, "US Standard Atmosphere 1976", US Government Printing Office.
"""

import numpy as np

from . import config

# Pre-compute the logarithms of the 86-120 km pressure and density columns once
# at import. Pressure and density decay ~exponentially with altitude, so
# log-linear interpolation (np.interp in log-space, then exp) is far more
# accurate than linear interpolation across the orders of magnitude involved.

LNP_EXT   = np.log(config.P_EXT)
LNRHO_EXT = np.log(config.RHO_EXT)
H_MESO    = config.H_EXT[0]      # 86000 m -- join between model and table


def atmosphere(h):
    """Return atmospheric properties at geometric altitude "h" [m].
    Parameters
    ----------
    h : float
        Geometric altitude [m] (clipped to [0, 120 000]).

        0-86 km   : US-76 geopotential-layer model.
        86-120 km : direct (log-linear) interpolation of the tabulated
                      T, P, rho (density is not derived from the ideal-gas law,
                      because the mean molecular weight varies above the
                      homopause).
    Returns
    -------
    dict
        Density 'rho' [kg m^-3], pressure 'P' [Pa], temperature 'T' [K],
        speed of sound 'a' [m s^-1], viscosity 'mu' [Pa s],
        geometric ('h') and geopotential ('H') altitudes.
    """
    
    h = float(np.clip(h, 0.0, config.H_ATM_MAX))

    # geopotential altitude 
    H = (config.r_earth * h) / (config.r_earth + h)

    if h <= H_MESO:
        #  0-86 km: US-76 geopotential-layer model 
        i = int(np.searchsorted(config.H_b, H, side='right') - 1)
        i = int(np.clip(i, 0, 7))

        dH = H - config.H_b[i]
        T = config.T_b[i] + config.L_b[i] * dH

        if abs(config.L_b[i]) < config.eps:               # isothermal layer
            P = config.P_b[i] * np.exp(-config.go * dH / (config.R_air * config.T_b[i]))
        else:                                              # gradient layer
            P = config.P_b[i] * (T / config.T_b[i]) ** (config.go / (-config.L_b[i] * config.R_air))

        rho = P / (config.R_air * T)

    else:
        # 86-120 km: direct table interpolation
        T   = float(np.interp(h, config.H_EXT, config.T_EXT))
        P   = float(np.exp(np.interp(h, config.H_EXT, LNP_EXT)))
        rho = float(np.exp(np.interp(h, config.H_EXT, LNRHO_EXT)))

    a  = np.sqrt(config.gamma * config.R_air * T)         # speed of sound
    mu = 1.458e-6 * T ** 1.5 / (T + 110.4)                # Sutherland's viscosity law

    return {'rho': rho, 'P': P, 'T': T, 'a': a, 'mu': mu, 'h': h, 'H': H}


def mach_from_velocity(v, h):
    """Mach number for velocity "v" [m s^-1] at altitude "h" [m]."""
    return v / atmosphere(h)['a']
