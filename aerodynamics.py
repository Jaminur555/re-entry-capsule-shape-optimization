"""Hypersonic aerodynamics via the modified Newtonian pressure model."""

import numpy as np
from . import config


class HypersonicAeroSolver:
    """Modified-Newtonian aerodynamic coefficient solver for axisymmetric bodies.
    Parameters
    ----------
    gamma : float, optional
        Ratio of specific heats (default from :mod:`capsule_opt.config`).
    """

    def __init__(self, gamma=config.gamma):
        self.gamma = gamma

    def cp_max(self, M):
        """Stagnation (maximum) pressure coefficient via the Rayleigh-Pitot relation.
        Parameters
        ----------
        M : float
            Freestream Mach number (valid for M >= 3).
        Returns
        -------
        float
            Stagnation pressure coefficient "Cp_max".
        """
        
        g, squ_M = self.gamma, M ** 2
        term1 = (((g + 1) ** 2) * squ_M) / ((4 * g * squ_M) - (2 * (g - 1)))
        term2 = g / (g - 1)
        term3 = (1 - g + (2 * g * squ_M)) / (g + 1)
        Cps = (2 / (g * squ_M)) * (((term1 ** term2) * term3) - 1)

        return Cps

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
            plus axial/normal coefficients, Mach, AoA, "Cp_max" and "L/D".
        """
        if M < 3:
            raise ValueError("Aerodynamic model is not valid below Mach 3.")

        alpha = np.deg2rad(AoA)
        Cp_max = self.cp_max(M)

        # capsule properties
        x_mid, z_mid = props["X_mid"], props["Z_mid"]
        L_ref, A_ref = props["L_ref"], props["A_ref"]
        nx, nz = props['n_x'], props['n_z']
        x_cg, z_cg = props['X_cg'], props['Z_cg']
        ds = props['ds']

        phi = np.linspace(0, np.pi * 2, 181)[:-1]  # azimuthal discretization (axisymmetric body)
        d_phi = phi[1] - phi[0]

        # dot product between freestream direction and local outward normal
        v_dot_n = nx[:, None] * np.cos(alpha) + nz[:, None] * np.sin(alpha) * np.cos(phi)
        sin_delta = -v_dot_n  # sign depends on adopted AoA / body-axis convention

        # modified Newtonian pressure applied only on windward surface panels
        Cp = Cp_max * (np.maximum(0, sin_delta) ** 2)

        # ring surface area formed by revolving a panel about the symmetry axis
        dA = z_mid[:, None] * d_phi * ds[:, None]

        # axial and normal force contributions in body-axis coordinates
        dCa = -Cp * nx[:, None] * dA
        dCn = -Cp * nz[:, None] * np.cos(phi) * dA

        Ca, Cn = np.sum(dCa), np.sum(dCn)

        # pitching moment about the center of gravity
        x_rltv = x_mid[:, None] - x_cg
        z_rltv = z_mid[:, None] * np.cos(phi) - z_cg
        dMy = (x_rltv * dCn) - (z_rltv * dCa)

        Cm = -np.sum(dMy)

        Ca, Cn, Cm = Ca / A_ref, Cn / A_ref, Cm / (A_ref * L_ref)

        # convert to drag and lift
        CD = Ca * np.cos(alpha) + Cn * np.sin(alpha)
        CL = -Ca * np.sin(alpha) + Cn * np.cos(alpha)

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
