"""Derived geometric and inertial properties of the capsule."""
import numpy as np
import scipy as sp

from . import config
from .geometry import capsule_profile
from .parameterization import cap_params, normalized_param

# Local-inclination part tags
PART_NOSE = 0   # blunt forebody (nose sphere + shoulder) -> Modified Newtonian 
PART_AFTERBODY = 1 # low-inclination afterbody (cone + rear) -> tangent-cone/Newtonian

def get_capsule_properties(rn, rs, r_theta, cg_params=(0.5, 0.7)):
    """Compute the capsule geometric properties and center-of-gravity model.
    Parameters
    ----------
    rn, rs, r_theta : float
        Normalized design parameters in [0, 1].
    cg_params : tuple of float, optional
        Normalized "(dx, dz)" offsets placing the CG relative to the volume
        centroid (axial offset over length, radial offset over local height).
    Returns
    -------
    dict
        Lengths, reference area, surface normals, panel discretization,
        volume, surface area, center of gravity and nose-sphere half-angle.
    """
    Lc = config.Lc_fixed
    Rn, Rs, theta_c = cap_params(rn, rs, r_theta)
    theta_c = np.deg2rad(theta_c)

    # ========== total length ==========

    arg = (config.Rm_fixed - Rs) / (Rn - Rs)
    arg = np.clip(arg, 0, 1)
    theta_ns_max = np.arcsin(arg)

    Yt_c    = (Rn - Rs) * np.sin(theta_ns_max)         # toroid centre y
    Y_t_end = Yt_c + Rs * np.cos(theta_c)              # toroid end y
    Y_c_end = Y_t_end - Lc * np.sin(theta_c)           # cone end y
    R_rs    = Y_c_end / np.cos(theta_c)                # rear sphere radius

    Lsp1 = Rn * (1 - np.cos(theta_ns_max))
    Lt = Rs * (np.cos(theta_ns_max) + np.sin(theta_c))
    Lsp2 = R_rs * (1 - np.sin(theta_c))
    L_total = Lsp1 + Lt + Lc + Lsp2
    L_ref = L_total

    # ========== geometry ================

    X, Y, breaks = capsule_profile(rn, rs, r_theta)

    # ========== maximum body radius ==========

    R_max = np.max(Y)

    # ========== reference frontal area ==========

    A_ref = np.pi * (R_max ** 2)

    # ========== surface panels ==========

    dx, dz = np.diff(X), np.diff(Y)
    ds = np.sqrt(dx ** 2 + dz ** 2)

    X_mid = 0.5 * (X[:-1] + X[1:])
    Z_mid = 0.5 * (Y[:-1] + Y[1:])

   #=========== local-inclination part tag =========
   # nose spher + shoulder = blunt forebody; cone + rear sphere = afterbody.

    panel_part = np.where(np.arange(X_mid.size) < breaks['cone_start'], 
                        PART_NOSE, PART_AFTERBODY)
 
   # ========== outwards unit normal ================

    ds_safe = np.where(ds < config.eps, config.eps, ds)
    nx, nz = -dz / ds_safe, dx / ds_safe

    # ==========Surface area ==========

    surf_area = 2.0 * np.pi * np.sum(Z_mid * ds)   # Pappus' centroid theorem

    # ========== center of gravity (offsets from the volume centroid) ==========

    sort_idx = np.argsort(X)
    X_s, Y_s = X[sort_idx], Y[sort_idx]  # ensure X is monotonic for integration

    volume = np.pi * sp.integrate.simpson(y=Y_s ** 2, x=X_s)
    x_centroid = (np.pi * sp.integrate.simpson(y=Y_s ** 2 * X_s, x=X_s)) / volume

    dx_over_L = normalized_param(config.DX_OVER_L_MAX, config.DX_OVER_L_MIN, cg_params[0])
    dx_com = dx_over_L * L_total
    x_cg = x_centroid + dx_com

    h_local = np.interp(x_cg, X, Y)
    dz_over_h = normalized_param(config.DZ_OVER_H_MAX, config.DZ_OVER_H_MIN, cg_params[1])
    z_cg = -dz_over_h * h_local

    return {
        "L_total": L_total,
        "R_max"  : R_max,
        "X": X,
        "Y": Y,
        "A_ref": A_ref,
        "X_mid": X_mid,
        "Z_mid": Z_mid,
        "ds" : ds,
        "n_x": nx,
        "n_z": nz,
        "X_cg" : x_cg,
        "Z_cg" : z_cg,
        "L_ref": L_ref,
        "theta_sp1" : theta_ns_max,
        "surf_area" : surf_area,
        "volume"    : volume,
        "panel_part": panel_part        
    }
