"""Uncertain-variable definitions for the robust-optimization study (uq/PLAN.md)."""

import numpy as np
from scipy import stats

# (name, frozen scipy distribution); model-form scales sized from V2 validation.
UNCERTAIN_VARS = [
    ("gamma0_deg", stats.norm(-2.0, 0.5)),    # entry FPA [deg]
    ("V0",         stats.norm(7830.0, 50.0)),  # entry speed [m/s]
    ("k_CD",       stats.norm(1.00, 0.05)),    # drag model-form scale
    ("k_CL",       stats.norm(1.00, 0.05)),    # lift model-form scale
]


def sample(n, seed=42):
    """(n, len(UNCERTAIN_VARS)) samples via independent inverse-CDF of a unit LHS."""
    u = stats.qmc.LatinHypercube(d=len(UNCERTAIN_VARS), seed=seed).random(n)
    return np.column_stack([dist.ppf(u[:, j]) for j, (_, dist) in enumerate(UNCERTAIN_VARS)])


if __name__ == "__main__":
    x = sample(10000)
    for j, (name, dist) in enumerate(UNCERTAIN_VARS):
        print(f"{name:10s} mean {x[:, j].mean():10.4f} (nom {dist.mean():10.4f}) "
              f"std {x[:, j].std(ddof=1):8.4f} (nom {dist.std():8.4f}) "
              f"range [{x[:, j].min():.4f}, {x[:, j].max():.4f}]")
