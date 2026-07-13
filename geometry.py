"""Meridional (side-profile) generation of the axisymmetric capsule."""

import numpy as np

from . import config
from .parameterization import cap_params


def capsule_profile(rn, rs, r_theta):
    """Build the capsule meridional profile "(X, Y)" for normalized parameters.

    The profile is assembled from four analytical segments: nose sphere,
    toroidal shoulder, conical frustum and rear sphere.

    Parameters
    ----------
    rn, rs, r_theta : float
        Normalized design parameters in [0, 1].

    Returns
    -------
    X, Y : ndarray
        Axial and radial coordinates of the meridional profile [m].

    Raises
    ------
    ValueError
        If the rear-sphere radius is non-positive or the geometry contains
        non-finite values (infeasible shape).
    """
    Lc = config.Lc_fixed
    Rn, Rs, theta_c = cap_params(rn, rs, r_theta)
    theta_c = np.deg2rad(theta_c)

    # ========== nose sphere ============
    # center = (Rn, 0)

    arg = (config.Rm_fixed - Rs) / (Rn - Rs)
    arg = np.clip(arg, 0, 1)
    theta_ns_max = np.arcsin(arg)
    theta_ns = np.linspace(0, theta_ns_max, 200)

    X_ns, Y_ns = Rn * (1 - np.cos(theta_ns)), Rn * np.sin(theta_ns)

    # ========== toroidal segment ==========

    Xt_c, Yt_c = Rn - (Rn - Rs) * np.cos(theta_ns_max), (Rn - Rs) * np.sin(theta_ns_max)
    theta_t_start, theta_t_final = np.pi - theta_ns_max, np.pi / 2 - theta_c
    theta_t = np.linspace(theta_t_start, theta_t_final, 100)

    X_t, Y_t = Xt_c + Rs * np.cos(theta_t), Yt_c + Rs * np.sin(theta_t)

    # ========== conical frustum ==========

    X_t_end, Y_t_end = X_t[-1], Y_t[-1]
    Lc_conical = np.linspace(0, Lc, 500)

    X_c, Y_c = X_t_end + Lc_conical * np.cos(theta_c), Y_t_end - Lc_conical * np.sin(theta_c)

    # ========== rear sphere ==========

    X_c_end, Y_c_end = X_c[-1], Y_c[-1]

    R_rs = Y_c_end / np.cos(theta_c)  # radius

    if R_rs <= 0:
        raise ValueError(f"Infeasible Geometry: Rear sphere radius R_rs={R_rs:.2f} is non-positive.")

    Xrs_c = X_c_end - R_rs * np.sin(theta_c)  # center of rear sphere (x_rs, 0)
    theta_rs = np.linspace(np.pi / 2 - theta_c, 0, 200)

    X_rs, Y_rs = Xrs_c + R_rs * np.cos(theta_rs), R_rs * np.sin(theta_rs)

    # ========== output ==========
    
    X = np.concatenate([X_ns, X_t, X_c, X_rs])
    Y = np.concatenate([Y_ns, Y_t, Y_c, Y_rs])

    if not np.all(np.isfinite(X)) or not np.all(np.isfinite(Y)):
        raise ValueError("Geometry contains non-finite values.")

    return X, Y
