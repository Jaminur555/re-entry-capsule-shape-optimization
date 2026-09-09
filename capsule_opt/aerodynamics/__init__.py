"""Hypersonic aerodynamics based on the Local Inclination Method (LIM):
pressure laws (modified Newtonian, tangent-cone, Prandtl-Meyer), lookup
tables, and force/moment integration."""

from .solver import HypersonicAeroSolver

__all__ = ["HypersonicAeroSolver"]