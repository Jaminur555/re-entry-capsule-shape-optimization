import numpy as np
from .. import config
from .lim import lim_pressure_coefficient
from .newtonian import cp_max

class HypersonicAeroSolver:
    """Local-inclination-method aerodynamic coefficient solver for axisymmetric bodies.

    Parameters
    ----------
    gamma : float, optional
        Ratio of specific heats (default from :mod:`capsule_opt.config`).
    """

    def __init__(self, gamma=config.gamma):
        self.gamma = gamma

    def coefficients(self, props, M, AoA):
        """Integrate axial, normal and moment coefficients over the body panels.

        Parameters
        ----------
        props : dict
            Capsule property dictionary (from
            :func:`capsule_opt.properties.get_capsule_properties`).
        M : float
            Freestream Mach number (must be >= 3).
        AoA : float
            Angle of attack [deg].

        Returns
        -------
        dict
            Drag "CD", lift "CL" and pitching moment "Cm" about the CG,
            plus axial/normal coefficients, Mach, AoA, "Cps" (Cp_max) and "L/D".
        """
        if M < 3:
            raise ValueError("Aerodynamic model is not valid below Mach 3.")

        alpha = np.deg2rad(AoA)
        Cp_max = cp_max(M)

        # capsule properties
        x_mid, z_mid = props["X_mid"], props["Z_mid"]
        L_ref, A_ref = props["L_ref"], props["A_ref"]
        nx, nz = props['n_x'], props['n_z']
        x_cg, z_cg = props['X_cg'], props['Z_cg']
        ds = props['ds']

        phi = np.linspace(0, np.pi * 2, 181)[:-1]       # azimuthal discretization
        d_phi = phi[1] - phi[0]

        # local flow inclination: sin(theta) = -Vinf.nhat
        v_dot_n = nx[:, None] * np.cos(alpha) + nz[:, None] * np.sin(alpha) * np.cos(phi)
        sin_delta = -v_dot_n
        theta_deg = np.rad2deg(np.arcsin(np.clip(sin_delta, -1.0, 1.0)))

        # Local-inclination Cp per panel
        Cp = lim_pressure_coefficient(theta_deg, props['panel_part'], M)

        # ring surface area formed by revolving a panel about the symmetry axis
        dA = z_mid[:, None] * d_phi * ds[:, None]

        # axial / normal force contributions (body axes)
        dCa = -Cp * nx[:, None] * dA
        dCn = -Cp * nz[:, None] * np.cos(phi) * dA
        Ca, Cn = np.sum(dCa), np.sum(dCn)

        # pitching moment about the CG
        x_rltv = x_mid[:, None] - x_cg
        z_rltv = z_mid[:, None] * np.cos(phi) - z_cg
        Cm = -np.sum((x_rltv * dCn) - (z_rltv * dCa))

        Ca, Cn, Cm = Ca / A_ref, Cn / A_ref, Cm / (A_ref * L_ref)

        # Body -> wind axes.the body x-axis points downstream (nose at the origin,
        # freestream Vinf = +x), so Ca is the downstream (drag-direction) axial force.
        # The lift perpendicular is therefore (sin a, -cos a): CL = Ca sin a - Cn cos a
        # -> positive L/D at positive AoA (verified ~0.30 at the trim AoA ~20 deg, matching Apollo).
        #  Drag is the Vinf projection as usual. 
        CD = Ca * np.cos(alpha) + Cn * np.sin(alpha)
        CL = Ca * np.sin(alpha) - Cn * np.cos(alpha)

        l_by_d = CL / CD if np.abs(CD) > config.eps else 0.00

        return {'CD': CD,
                'CL': CL,
                'Cm': Cm,
                'Cn': Cn,
                'Ca': Ca,
                'Ma': M,
                'Aoa': AoA,
                'Cps': Cp_max,
                'L/D': l_by_d}
