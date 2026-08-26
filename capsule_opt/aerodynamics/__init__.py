"""
Hypersonic aerodynamics package based on the
Local Inclination Method (LIM).

Modules
-------
solver : Aerodynamic force and moment integration.
newtonian : Modified-Newtonian pressure model (cp_max, cp_newtonian).
shocks : Tangent-cone (Taylor-Maccoll) pressure model.
expansion : Prandtl-Meyer expansion model.
cp_tables : Cone / Prandtl-Meyer lookup-table build, load, and interpolation.
lim : Local-inclination pressure-law dispatcher.
"""

from .solver import HypersonicAeroSolver

__all__ = ["HypersonicAeroSolver"]