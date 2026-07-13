"""Aerodynamic coefficient database (Mach x AoA grid) and interpolators."""

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from . import config
from .aerodynamics import HypersonicAeroSolver
from .properties import get_capsule_properties
from .trim_solver import trim_alpha


def build_aero_database(rn, rs, r_theta, cg_params=(0.5, 0.7)):
    """Build the Mach x AoA aerodynamic coefficient database for a capsule.

    Evaluates the modified-Newtonian solver over the configured Mach/AoA grid,
    builds "RegularGridInterpolator" objects for CD, CL and Cm, and computes
    the trim AoA at each Mach node.

    Parameters
    ----------
    rn, rs, r_theta : float
        Normalized design parameters in [0, 1].
    cg_params : tuple, optional
        CG placement offsets.

    Returns
    -------
    tuple
        "(interp_CD, interp_CL, interp_CM, alpha_trim_arr)" -- three
        interpolators and the array of trim AoA per Mach node [deg].
    """
    props_local = get_capsule_properties(rn, rs, r_theta, cg_params=cg_params)
    solver_local = HypersonicAeroSolver()

    DB_CD_local = np.zeros((len(config.mach_nodes), len(config.aoa_nodes)))
    DB_CL_local = np.zeros_like(DB_CD_local)
    DB_CM_local = np.zeros_like(DB_CD_local)

    for i, M in enumerate(config.mach_nodes):
        for j, aoa in enumerate(config.aoa_nodes):
            r = solver_local.coefficients(props_local, M, aoa)
            DB_CD_local[i, j] = r['CD']
            DB_CL_local[i, j] = r['CL']
            DB_CM_local[i, j] = r['Cm']

    interp_CD = RegularGridInterpolator((config.mach_nodes, config.aoa_nodes), DB_CD_local, method='linear', bounds_error=False, fill_value=None)
    interp_CL = RegularGridInterpolator((config.mach_nodes, config.aoa_nodes), DB_CL_local, method='linear', bounds_error=False, fill_value=None)
    interp_CM = RegularGridInterpolator((config.mach_nodes, config.aoa_nodes), DB_CM_local, method='linear', bounds_error=False, fill_value=None)

    alpha_trim_arr = np.zeros(len(config.mach_nodes))
    for i, M in enumerate(config.mach_nodes):
        at = trim_alpha(M, interp_CM)
        alpha_trim_arr[i] = at if at is not None else 20.0

    return interp_CD, interp_CL, interp_CM, alpha_trim_arr
