"""Capsule shape multi-objective optimization package.

A modular reorganization of the monolithic "run_optimization.py" into focused
modules, one per concern:

    config           : centralized tunable constants (single source of truth)
    parameterization ; normalized -> physical design-parameter mapping
    geometry         : meridional (side-profile) capsule generation
    properties       : derived geometric / inertial properties
    atmosphere       : US Standard Atmosphere 1976 (0-86 km) extended to 120 km
    aerodynamics     : modified-Newtonian hypersonic coefficient solver
    heating          : effective nose radius and stagnation / shoulder heating
    trim_solver      : trim angle-of-attack solver
    aero_database    : Mach x AoA coefficient grid, interpolators
    trajectory       : entry equations of motion and trajectory propagation
    metrics          : volumetric efficiency, heat load, load factor
    constraints      : constraint evaluator, stability / trim margins
    objectives       : single-design evaluation (the integration point)
    optimization     : pymoo CapsuleOptimization problem
    plotting         : Pareto-front visualization and trade-off reporting
    main             : run-the-optimization orchestrator

Typical usage::

    from capsule_opt import main
    main()
"""

from .main import main
from .objectives import evaluate_shape
from .optimization import CapsuleOptimization

__all__ = ["main", "CapsuleOptimization", "evaluate_shape"]
