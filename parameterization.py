"""Normalized -> physical design-parameter mapping for the capsule."""

import numpy as np
from . import config


def normalized_param(x_max, x_min, r):
    """Map a normalized coordinate "r" in [0, 1] to the range [x_min, x_max].

    Parameters
    ----------
    x_max, x_min : float
        Endpoints of the physical range; 'x_max' must be >= 'x_min'.
    r : float
        Normalized coordinate (clipped to [0, 1]).

    Returns
    -------
    float
        Physical value "x_min + r * (x_max - x_min)".
    """
    r = np.clip(r, 0.0, 1.0)
    if x_max < x_min:
        raise ValueError("The first argument should be the maximum value")
    return x_min + r * (x_max - x_min)


def cap_params(rn, rs, r_theta):
    """Convert normalized design parameters to physical capsule dimensions.

    Parameters
    ----------
    rn, rs, r_theta : float
        Normalized parameters in [0, 1] for nose radius, shoulder radius and
        afterbody cone half-angle respectively.

    Returns
    -------
    tuple of float
        "(Rn, Rs, theta_c)"  nose radius [m], shoulder radius [m] and
        afterbody cone half-angle [deg].
    """
    Rn = normalized_param(config.RN_MAX, config.RN_MIN, rn)
    Rs = normalized_param(config.RS_MAX, config.RS_MIN, rs)
    theta_c = normalized_param(config.THETA_C_MAX, config.THETA_C_MIN, r_theta)
    return Rn, Rs, theta_c
