from .. import config
import numpy as np

def cp_max(M, gamma=config.gamma):
    """Stagnation (maximum) pressure coefficient via the Rayleigh-Pitot relation."""
    g, squ_M = gamma, M ** 2
    term1    = (((g + 1) ** 2) * squ_M) / ((4 * g * squ_M) - (2 * (g - 1)))
    term2    = g / (g - 1)
    term3    = (1 - g + (2 * g * squ_M)) / (g + 1)

    return (2 / (g * squ_M)) * (((term1 ** term2) * term3) - 1)

def cp_newtonian(M, theta_deg, gamma=config.gamma):
    """Blunt / high-hyp windward Cp: modified Newtonian (theta in degrees)."""
    
    return cp_max(M, gamma) * np.sin(np.deg2rad(theta_deg)) ** 2