"""Constraint evaluation and static-stability / trim margins."""

import numpy as np

from . import config
from .metrics import heat_load, load_factor
from .properties import get_capsule_properties
from .trim_solver import trim_alpha


def pitch_stability_margin(rn, rs, r_theta, cg_params,
                           mach_check=10,
                           delta_alpha=0.5,
                           interp_CM=None):
    """Pitch-stability derivative "dCm/da" at trim (must be <= 0 to be stable).

    Parameters
    ----------
    rn, rs, r_theta : float
        Normalized design parameters in [0, 1].
    cg_params : tuple
        CG placement offsets.
    mach_check : float, optional
        Mach number at which to check stability.
    delta_alpha : float, optional
        Finite-difference increment [deg].
    interp_CM : callable, optional
        Pitching-moment interpolator from the aero database.

    Returns
    -------
    float
        "dCm/da" at trim; "1.0" (treated as unstable) if no trim exists.
    """
    at = trim_alpha(mach_check, interp_CM)
    if at is None:
        return 1.0  # no trim -> treat as unstable (positive violation)
    Mc   = float(np.clip(mach_check, config.mach_nodes[0], config.mach_nodes[-1]))
    cm_p = float(interp_CM(np.array([[Mc, at + delta_alpha]]))[0])
    cm_n = float(interp_CM(np.array([[Mc, at - delta_alpha]]))[0])

    return (cm_p - cm_n) / (2.0 * delta_alpha)


def trim_on_nose_margin(rn, rs, r_theta, cg_params, interp_CM=None):
    """Margin between the trim AoA and the nose-sphere half-angle [deg].

    A negative value means trim occurs before the nose, i.e. on-nose trim
    (infeasible). Returns "90.0" (infeasible) if no trim exists.
    """
    at = trim_alpha(float(np.clip(10.0, config.mach_nodes[0], config.mach_nodes[-1])), interp_CM)

    if at is None:
        return 90.0  # no trim -> infeasible
    props = get_capsule_properties(rn, rs, r_theta, cg_params=cg_params)
    theta_N_deg = np.rad2deg(props['theta_sp1'])
    return float(at) - theta_N_deg


def evaluate_constraints(traj, rn, rs, r_theta, m,
                         cg_params=(0.5, 0.7),
                         shape_props=None, interp_CD=None, interp_CM=None):
    """Evaluate all design constraints for a trajectory and capsule shape.

    Constraints are written as "g(x) <= 0" (feasible). They bound the peak
    stagnation / shoulder heat flux, the load factor, the pitch stability and
    the trim-on-nose margin.

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
    interp_CD, interp_CM : callable, optional
        Drag / pitching-moment interpolators from the aero database.

    Returns
    -------
    dict
        "feasible" flag, the five normalized constraint values and a
        "details" dict with the underlying peak quantities.
    """
    hl = heat_load(traj, rn, rs, r_theta)
    lf = load_factor(traj, rn, rs, r_theta, m=m, cg_params=cg_params,
                     shape_props=shape_props, interp_CD=interp_CD)

    g_q_stag  = (hl['q_stag_max'] - config.C_Q_STAG_MAX) / config.C_Q_STAG_MAX
    g_q_shldr = (hl['q_shldr_max'] - config.C_Q_SHLDR_MAX) / config.C_Q_SHLDR_MAX
    g_n_max   = (lf['n_max'] - config.C_N_MAX) / config.C_N_MAX
    g_pitch   = pitch_stability_margin(rn, rs, r_theta, cg_params, interp_CM=interp_CM)
    g_nose    = trim_on_nose_margin(rn, rs, r_theta, cg_params, interp_CM=interp_CM)

    feasible = (g_q_stag  <= 0.0 and
                g_q_shldr <= 0.0 and
                g_n_max   <= 0.0 and
                g_pitch   <= config.C_CMA_MAX and
                g_nose    <= 0.0)

    return {
        'feasible'        : feasible,
        'g_q_stag'        : g_q_stag,
        'g_q_shldr'       : g_q_shldr,
        'g_n_max'         : g_n_max,
        'g_pitch_stable'  : g_pitch,
        'g_trim_on_nose'  : g_nose,
        'details'         : {
            'q_stag_max'  : hl['q_stag_max'],
            'q_shldr_max' : hl['q_shldr_max'],
            'n_max'       : lf['n_max'],
            'Qs'          : hl['Qs'],
            'sg'          : traj['sg']
        }
    }
