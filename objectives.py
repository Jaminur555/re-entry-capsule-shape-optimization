"""Objective integration: evaluate one capsule design end-to-end.

This module is the integration point that wires together the property, aero,
trajectory, metric and constraint modules to score a single candidate shape.
"""

from . import config
from .aero_database import build_aero_database
from .constraints import evaluate_constraints
from .metrics import heat_load, volumetric_efficiency
from .parameterization import cap_params
from .properties import get_capsule_properties
from .trajectory import propagate_trajectory


def evaluate_shape(rn, rs, r_theta,
                   cg_params=(0.5, 0.7),
                   m=None,
                   entry_conditions=None):
    """Evaluate objectives and constraints for one capsule design.

    Assembles the full analysis pipeline for a candidate shape: geometric
    properties -> aero database -> trajectory -> performance metrics ->
    constraints.

    Parameters
    ----------
    rn, rs, r_theta : float
        Normalized design parameters in [0, 1].
    cg_params : tuple, optional
        CG placement offsets.
    m : float or None, optional
        Vehicle mass [kg]. If "None" (default), the mass is derived from the
        capsule volume via the Apollo-based constant density
        "config.CAPSULE_DENSITY". Pass a float
        to force a fixed mass (e.g. for a sensitivity study).
    entry_conditions : dict, optional
        Overrides for the entry-interface initial conditions (any of "ho",
        "Vo", "gamma0_deg", "lat0_deg", "lon0_deg", "chi0_deg").

    Returns
    -------
    dict
        "objectives" ('eta_V', 'Qs', 'sg'), "constraints", the
        "feasible" flag, the full "traj" history and a "shape" summary.
    """

    ec = {'ho': 120000.0, 'Vo': 7830.0, 'gamma0_deg': -2.0,
          'lat0_deg': -23.75, 'lon0_deg': 225.5, 'chi0_deg': 49.59}
    if entry_conditions is not None:
        ec.update(entry_conditions)

    shape_props = get_capsule_properties(rn, rs, r_theta, cg_params=cg_params)
    interp_CD, interp_CL, interp_CM, alpha_trim_arr = build_aero_database(rn, rs, r_theta, cg_params)

    # vehicle mass: Apollo-based constant-density model
    # unless a fixed mass is supplied for a sensitivity study.
    if m is None:
        m = config.CAPSULE_DENSITY * shape_props['volume']

    # =================== trajectory ===================

    traj = propagate_trajectory(
        rn, rs, r_theta, m,
        ho=ec['ho'], Vo=ec['Vo'], gamma0_deg=ec['gamma0_deg'],
        lat0_deg=ec['lat0_deg'], lon0_deg=ec['lon0_deg'], chi0_deg=ec['chi0_deg'],
        cg_params=cg_params,
        shape_props=shape_props,
        interp_CD=interp_CD, interp_CL=interp_CL, alpha_trim_arr=alpha_trim_arr
    )

    # =================== objectives ===================

    eta_V = volumetric_efficiency(rn, rs, r_theta, cg_params=cg_params)
    hl    = heat_load(traj, rn, rs, r_theta)
    sg    = traj['sg']

    # =================== constraints ===================
    con = evaluate_constraints(traj, rn, rs, r_theta, m=m, cg_params=cg_params,
                               shape_props=shape_props, interp_CD=interp_CD, interp_CM=interp_CM)

    # =================== shape summary ===================
    Rn, Rs, theta_c = cap_params(rn, rs, r_theta)

    return {
        'objectives': {
            'eta_V': eta_V,
            'Qs'   : hl['Qs'],
            'sg'   : sg
        },
        'constraints': con,
        'feasible'   : con['feasible'],
        'traj'       : traj,
        'shape'      : {
            'rn': rn, 'rs': rs, 'r_theta': r_theta,
            'Rn': Rn, 'Rs': Rs, 'theta_c': theta_c,
            'L_total': shape_props['L_total'],
            'A_ref'  : shape_props['A_ref'],
            'volume' : shape_props['volume'],
            'm'      : m
        }
    }
