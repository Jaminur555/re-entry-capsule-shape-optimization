"""Trajectory propagation of the entry capsule through the atmosphere."""

import numpy as np
from scipy.integrate import solve_ivp
from . import config
from .atmosphere import atmosphere, mach_from_velocity


def bank_angle_sigma(m, L_force, V, h, gamma, lat, chi):
    """Bank angle 'sigma' that cancels the out-of-plane acceleration (no-skip).

    Parameters
    ----------
    m : float
        Vehicle mass [kg].
    L_force : float
        Aerodynamic lift force magnitude [N].
    V, h, gamma, lat, chi : float
        Speed [m s^-1], altitude [m], flight-path angle, latitude and heading
        angle [rad].

    Returns
    -------
    float
        Bank angle [rad].
    """
    r     = config.r_earth + h
    g     = config.go * (config.r_earth / r) ** 2
    Vc2   = g * r
    cos_s = (m / max(L_force, 1e-6)) * (g * (1 - V ** 2 / Vc2) * np.cos(gamma) -
                                 2 * config.omega_earth * V * np.cos(lat) * np.sin(chi))
    return float(np.arccos(np.clip(cos_s, -1, 1)))


def entry_eom(t, state, m, A_ref, interp_CD, interp_CL, alpha_trim_arr):
    """Point-mass entry equations of motion over a spherical, rotating Earth.

    State vector "state = [h, V, gamma, lat, lon, chi]". Aerodynamic
    coefficients are read from the supplied interpolators at the trim AoA.
    """
    h, V, gamma, lat, lon, chi = state
    h = max(float(h), 0.0)

    # atmospheric state
    atm = atmosphere(h)
    rho = atm['rho']
    q   = 0.5 * rho * V ** 2          # dynamic pressure [Pa]
    M   = float(np.clip(V / atm['a'], config.mach_nodes[0], config.mach_nodes[-1]))

    # trim AoA and aerodynamic coefficients (from database)
    at = float(np.clip(np.interp(M, config.mach_nodes, alpha_trim_arr), config.aoa_nodes[0], config.aoa_nodes[-1]))
    CD = float(interp_CD(np.array([[M, at]]))[0])
    CL = float(interp_CL(np.array([[M, at]]))[0])

    D = q * A_ref * CD                # drag force [N]
    L = q * A_ref * abs(CL)           # lift force magnitude [N]

    # bank angle (no-skip condition)
    sigma = bank_angle_sigma(m, L, V, h, gamma, lat, chi)

    # altitude-corrected gravity
    r     = config.r_earth + h
    g     = config.go * (config.r_earth / r) ** 2

    # equations of motion
    dh     = V * np.sin(gamma)
    dV     = - D / m - g * np.sin(gamma)
    dlat   = V * np.cos(gamma) * np.cos(chi) / r
    dlon   = V * np.sin(chi) * np.cos(gamma) / (r * np.cos(lat))
    dchi   = ((L * np.sin(sigma) / (m * V * np.cos(gamma))) + (V * np.cos(gamma) * np.sin(chi) * np.tan(lat) / r)
            - 2 * config.omega_earth * (np.cos(chi) * np.tan(gamma) * np.cos(lat) - np.sin(lat)))
    dgamma = ((L * np.cos(sigma) / (m * V)) - (g / V - V / r) * np.cos(gamma) + 2.0 * config.omega_earth * np.cos(lat) * np.sin(chi))

    return [dh, dV, dgamma, dlat, dlon, dchi]


def mach3_event(t, state, *args):
    """Event: Mach-3 crossing (terminal, decreasing) -- ends the integration."""
    return mach_from_velocity(state[1], max(state[0], 0.0)) - 3

mach3_event.terminal = True
mach3_event.direction = -1


def propagate_trajectory(rn, rs, r_theta, m,
                         ho=120000.0,
                         Vo=7830.0,
                         gamma0_deg=-2.0,
                         lat0_deg=-23.75,
                         lon0_deg=225.5,
                         chi0_deg=49.59,
                         cg_params=(0.5, 0.7),
                         t_max=7200.0,
                         dt_max=2.0,
                         shape_props=None, interp_CD=None, interp_CL=None, alpha_trim_arr=None):
    """Integrate the entry trajectory from the entry interface to the Mach-3 event.

    Parameters
    ----------
    rn, rs, r_theta : float
        Normalized design parameters in [0, 1].
    m : float
        Vehicle mass [kg].
    ho, Vo, gamma0_deg, lat0_deg, lon0_deg, chi0_deg : float, optional
        Entry-interface initial conditions (altitude [m], speed [m s^-1] and
        angles [deg]).
    cg_params : tuple, optional
        CG placement offsets.
    t_max, dt_max : float, optional
        Maximum integration time [s] and max step [s].
    shape_props : dict, optional
        Capsule properties (provides "A_ref").
    interp_CD, interp_CL : callable, optional
        Drag / lift interpolators from the aero database.
    alpha_trim_arr : ndarray, optional
        Trim AoA per Mach node.

    Returns
    -------
    dict
        Time histories of the state variables plus Mach, dynamic pressure, trim
        AoA, great-circle ground-track range "sg" [m] and the solver status.
    """
    state0 = [ho, Vo,
              np.deg2rad(gamma0_deg), np.deg2rad(lat0_deg),
              np.deg2rad(lon0_deg), np.deg2rad(chi0_deg)]

    sol = solve_ivp(
        entry_eom, [0.0, t_max], state0,
        args=(m, shape_props['A_ref'], interp_CD, interp_CL, alpha_trim_arr),
        method='RK45',
        events=mach3_event,
        max_step=dt_max,
        rtol=1e-6, atol=1e-8,
        dense_output=False
    )

    h_arr   = sol.y[0]
    V_arr   = sol.y[1]
    gam_arr = sol.y[2]
    lat_arr = sol.y[3]
    lon_arr = sol.y[4]
    chi_arr = sol.y[5]

    h_clipped = np.maximum(h_arr, 0.0)
    rho_arr   = np.array([atmosphere(h)['rho'] for h in h_clipped])
    q_arr     = 0.5 * rho_arr * V_arr ** 2
    Mach_arr  = np.array([mach_from_velocity(V_arr[k], h_clipped[k]) for k in range(len(sol.t))])

    M_clipped = np.clip(Mach_arr, config.mach_nodes[0], config.mach_nodes[-1])
    at_arr    = np.interp(M_clipped, config.mach_nodes, alpha_trim_arr)
    at_arr    = np.clip(at_arr, config.aoa_nodes[0], config.aoa_nodes[-1])

    # great-circle ground track
    phi0, phi1 = lat_arr[0], lat_arr[-1]
    dl = lon_arr[-1] - lon_arr[0]
    sg = config.r_earth * np.arccos(np.clip(np.sin(phi0) * np.sin(phi1) +
                                     np.cos(phi0) * np.cos(phi1) * np.cos(dl), -1, 1))

    return {
        't'         : sol.t,
        'h'         : h_arr,
        'V'         : V_arr,
        'gamma_deg' : np.rad2deg(gam_arr),
        'lat_deg'   : np.rad2deg(lat_arr),
        'lon_deg'   : np.rad2deg(lon_arr),
        'chi_deg'   : np.rad2deg(chi_arr),
        'Mach'      : Mach_arr,
        'q_dyn'     : q_arr,
        'alpha_trim': at_arr,
        'sg'        : sg,
        'status'    : sol.message
    }
