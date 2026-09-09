"""US Standard Atmosphere 1976, 0-120 km.

Below 86 km: US-76 geopotential-layer model. Above the homopause the variable
mean molecular weight means rho cannot be recovered from P and T with the
sea-level gas constant, so the 86-120 km region is interpolated directly
from the tabulated T, P, rho in config (which joins the layer model
continuously at 86 km).

Reference: NOAA/NASA, "US Standard Atmosphere 1976".
"""

import numpy as np

from .. import config

# Log-space interpolation of the 86-120 km P and rho columns: both decay
# ~exponentially with altitude, far beyond linear-interpolation accuracy.

LNP_EXT   = np.log(config.P_EXT)
LNRHO_EXT = np.log(config.RHO_EXT)
H_MESO    = config.H_EXT[0]      # 86000 m -- join between model and table


def atmosphere(h):
    """Atmospheric properties at geometric altitude "h" [m] (clipped to
    [0, 120 000]).

    Returns
    -------
    dict
        'rho' [kg m^-3], 'P' [Pa], 'T' [K], 'a' [m s^-1], 'mu' [Pa s],
        geometric 'h' and geopotential 'H' [m].
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
