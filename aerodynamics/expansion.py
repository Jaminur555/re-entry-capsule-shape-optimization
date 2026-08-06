from .. import config
import numpy as np


def nu(M, gamma=config.gamma):
    """Prandtl-Meyer function [rad]."""
    if M <= 1.0:
        return 0.0
    k = np.sqrt((gamma + 1) / (gamma - 1))
    return k * np.arctan(np.sqrt((M ** 2 - 1) / (k ** 2))) - np.arctan(np.sqrt(M ** 2 - 1.0))


def nu_inv(nu_val, gamma=config.gamma):
    """Inverse Prandtl-Meyer: M from nu, by bisection on M in [1, 60]."""
    if nu_val <= 0.0:
        return 1.0
    lo, hi = 1.0, 60.0
    f_lo   = nu(lo, gamma) - nu_val
    if nu(hi, gamma) - nu_val < 0.0:
        return hi
    for _ in range(100):
        mid   = 0.5 * (lo + hi)
        f_mid = nu(mid, gamma) - nu_val
        if abs(f_mid) < config.eps:
            break
        if f_lo * f_mid < 0.0:
            hi = mid
        else:
            lo, f_lo = mid, f_mid
    return mid


def cp_prandtl_meyer(M, theta_turn, gamma=config.gamma):
    """Lee-side Cp for an expansion of magnitude theta_turn >= 0 (rad)."""
    theta_turn = abs(theta_turn)
    M_loc      = nu_inv(nu(M, gamma) + theta_turn, gamma)
    p_ratio    = ((1 + 0.5 * (gamma - 1) * M ** 2) /
               (1 + 0.5 * (gamma - 1) * M_loc ** 2)) ** (gamma / (gamma - 1))
    Cp         = (2.0 / (gamma * M ** 2)) * (p_ratio - 1.0)
    return max(Cp, -2.0 / (gamma * M ** 2))     # vacuum floor
