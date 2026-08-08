from .. import config
import numpy as np
from .shocks import cp_tangent_cone
from .newtonian import cp_newtonian
from .expansion import cp_prandtl_meyer



M_LOW_MAX, M_HIGH_MIN   = config.M_LOW_MAX, config.M_HIGH_MIN
CONE_MACHS, CONE_THETAS = config.CONE_MACHS, config.CONE_THETAS
CONE_CACHE = config.CONE_CACHE


def build_cone_table(gamma=config.gamma, quiet=False):
    """Tabulate cp_tangent_cone over CONE_MACHS x CONE_THETAS (theta >= 0, deg).
       Detached entries (None) are left as NaN; clean_cone_table fills them."""
    table = np.full((CONE_MACHS.size, CONE_THETAS.size), np.nan)
    for i, M in enumerate(CONE_MACHS):
        for j, th in enumerate(CONE_THETAS):
            cp = cp_tangent_cone(M, th, gamma)
            if cp is not None:
                table[i, j] = cp
        if not quiet:
            print(f"  cone table M={M:4.1f}  ({i + 1}/{CONE_MACHS.size})")
    return table


def clean_cone_table(table, gamma=config.gamma):
    """Post-process the raw cone table to be NaN-free and monotone in theta:
       (1) detached cells (NaN) <- modified-Newtonian Cp (Dirkx's detached-shock
       fallback); (2) per Mach row, walk from high theta downwards and clamp
       Cp[j] = min(Cp[j], Cp[j+1]) to remove low-theta spikes (spurious
       strong-shock-branch roots near the Mach angle). Correct cells untouched."""
    cleaned = np.array(table, dtype=float)
    for i, M in enumerate(CONE_MACHS):
        row      = cleaned[i]
        nan_mask = np.isnan(row)
        if nan_mask.any():
            row[nan_mask] = cp_newtonian(M, CONE_THETAS[nan_mask], gamma)
        for j in range(CONE_THETAS.size - 2, -1, -1):
            if row[j] > row[j + 1]:
                row[j] = row[j + 1]
        cleaned[i] = row
    return cleaned


def load_cone_table():
    """Load the cached cone table; build + save it on first use (one-time, ~3-6 min).

    The cache is a human-readable CSV (rows = Mach, columns = theta-deg, with
    header labels) so it can be inspected or edited in any spreadsheet. Mach and
    theta grids still come from config; only the Cp grid is stored. Resolved
    against the PARENT package dir (capsule_opt/) so a single copy is reused."""
    import os
    path = os.path.join(os.path.dirname(os.path.dirname(__file__)), CONE_CACHE)
    if os.path.exists(path):
        return _read_cone_csv(path)
    print(f"Building tangent-cone table (one-time, ~3-6 min) -> {path}")
    table = clean_cone_table(build_cone_table())
    _write_cone_csv(path, table)
    return table


def _read_cone_csv(path):
    """Read a labeled cone-table CSV -> (nMach, nTheta) Cp array."""
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return data[:, 1:]                        # drop the Mach-label first column


def _write_cone_csv(path, table):
    """Write the cone table as a labeled CSV (Mach rows, theta-deg columns)."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("Mach\\theta_deg," + ",".join(f"{t:g}" for t in CONE_THETAS) + "\n")
        for i, M in enumerate(CONE_MACHS):
            f.write(f"{M:g}," + ",".join(f"{v:.8f}" for v in table[i]) + "\n")



PRANDTL_MACHS, PRANDTL_THETAS = config.PRANDTL_MACHS, config.PRANDTL_THETAS

def build_prandtl_table(gamma=config.gamma):
    """Tabulate cp_prandtl_meyer over PRANDTL_MACHS x PRANDTL_THETAS (|theta| in deg).
       Cheap (~1 s), so built eagerly at import -- no disk cache needed."""
    table = np.zeros((PRANDTL_MACHS.size, PRANDTL_THETAS.size))
    for i, M in enumerate(PRANDTL_MACHS):
        for j, th in enumerate(PRANDTL_THETAS):
            table[i, j] = cp_prandtl_meyer(M, np.deg2rad(th), gamma)
    return table
