"""Performance metrics: volumetric efficiency, integrated heat load, load factor."""

import numpy as np
import scipy as sp

from . import config
from .atmosphere import atmosphere
from .heating import effective_nose_radius
from .parameterization import cap_params
from .properties import get_capsule_properties


def volumetric_efficiency(rn, rs, r_theta, cg_params=(0.5, 0.7)):
    """Dimensionless volumetric efficiency "eta_V = 6*sqrt(pi)*V/S^1.5".

    Equals 1 for a sphere; smaller values mean more wetted area per enclosed
    volume.

    Parameters
    ----------
    rn, rs, r_theta : float
        Normalized design parameters in [0, 1].
    cg_params : tuple, optional
        CG placement offsets.

    Returns
    -------
    float
        Volumetric efficiency "eta_V" in (0, 1].
    """
    props = get_capsule_properties(rn, rs, r_theta, cg_params=cg_params)
    S, V = props['surf_area'], props['volume']

    return 6.0 * np.sqrt(np.pi) * (V / (S ** 1.5))


def heat_load(traj, rn, rs, r_theta, Tw_by_Taw=0):
    """Integrated stagnation heat load and peak stagnation / shoulder fluxes.

    Parameters
    ----------
    traj : dict
        Trajectory history (from
        :func:'capsule_opt.trajectory.propagate_trajectory').
    rn, rs, r_theta : float
        Normalized design parameters in [0, 1].
    Tw_by_Taw : float, optional
        Wall-to-adiabatic-wall temperature ratio in [0, 1].

    Returns
    -------
    dict
        Integrated heat load "Qs" [J m^-2], peak stagnation / shoulder fluxes
        [W m^-2] and the time-resolved flux arrays.
    """
    props = get_capsule_properties(rn, rs, r_theta, cg_params=(0.5, 0.7))
    Rn, Rs, _ = cap_params(rn, rs, r_theta)
    theta_sp1 = props['theta_sp1']

    t = traj['t']
    h = traj['h']
    V = traj['V']
    M = traj['Mach']
    at = traj['alpha_trim']

    h_clipped = np.maximum(h.astype(float), 0.0)
    rho_arr   = np.array([atmosphere(h)['rho'] for h in h_clipped])

    # stagnation heat flux
    R_for_k    = effective_nose_radius(Rn, Rs)
    K_val      = (1.83e-4 / np.sqrt(R_for_k))
    q_stag_arr = K_val * (rho_arr ** 0.5) * (V.astype(float) ** 3)

    # shoulder heat flux
    c1, c2, c3, c4, c5 = -0.0006, 0.0185, -0.5321, -0.2939, 1.3630
    rs_over_rm = Rs / config.Rm_fixed
    ratio_arr = (c1 * M + c2 * at + c3 * rs_over_rm + c4 * theta_sp1 + c5)
    ratio_arr = np.maximum(ratio_arr, 0.0)
    q_shldr_arr = ratio_arr * q_stag_arr

    Qs = sp.integrate.trapezoid(q_stag_arr, t)

    return {
        'Qs'         : float(Qs),
        'q_stag_max' : float(np.max(q_stag_arr)),
        'q_shldr_max': float(np.max(q_shldr_arr)),
        'q_stag_arr' : q_stag_arr,
        'q_shldr_arr': q_shldr_arr
    }


def load_factor(traj, rn, rs, r_theta, m,
                cg_params=(0.5, 0.7),
                shape_props=None, interp_CD=None):
    """Peak aerodynamic load factor "n" (in g) along the trajectory.

    Parameters
    ----------
    traj : dict
        Trajectory history.
    rn, rs, r_theta : float
        Normalized design parameters in [0, 1].
    m : float
        Vehicle mass [kg].
    cg_params : tuple, optional
        CG placement offsets.
    shape_props : dict, optional
        Capsule properties (provides "A_ref").
    interp_CD : callable, optional
        Drag interpolator from the aero database.

    Returns
    -------
    dict
        Peak load factor "n_max" [g], the time-resolved array and the
        time/altitude at which the peak occurs.
    """
    h = traj['h']
    V = traj['V']
    M = traj['Mach']
    at = traj['alpha_trim']

    A_ref = shape_props['A_ref']

    h_clipped = np.maximum(h.astype(float), 0.0)
    rho_arr   = np.array([atmosphere(h)['rho'] for h in h_clipped])

    M_clipped = np.clip(M.astype(float), config.mach_nodes[0], config.mach_nodes[-1])
    at_clipped = np.clip(at.astype(float), config.aoa_nodes[0], config.aoa_nodes[-1])

    CD_arr = np.array([float(interp_CD(np.array([[M_clipped[k], at_clipped[k]]]))[0]) for k in range(len(h))])
    n_arr  = (0.5 * rho_arr * V.astype(float) ** 2 * CD_arr * A_ref) / (m * config.go)

    idx = int(np.argmax(n_arr))

    return {
        'n_max'  : float(n_arr[idx]),
        'n_arr'  : n_arr,
        't_n_max': float(traj['t'][idx]),
        'h_n_max': float(traj['h'][idx])
    }
