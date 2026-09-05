import numpy as np
import os

from .. import config
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
    """Post-process the raw cone table to be NaN-free and monotone in theta.

    (1) Remove low-theta spikes among the attached cells (spurious
        strong-shock-branch roots near the Mach angle): walking right-to-left,
        no attached cell may exceed the next attached cell to its right.
    (2) Fill detached cells (NaN) with the LARGER of the modified-Newtonian
        fallback (Dirkx's detached-shock treatment) and the last attached
        cone Cp. Taking the max keeps Cp from DROPPING when the shock
        detaches -- a plain Newtonian fill sits below the last attached
        value at high theta (e.g. M=10: fill ~1.30 vs attached 1.465 at
        55 deg), which used to drag valid cells down via the monotone
        enforcement and biased the afterbody Cp low.
    """
    cleaned = np.array(table, dtype=float)
    n = CONE_THETAS.size
    for i, M in enumerate(CONE_MACHS):
        row = cleaned[i]

        # (1) spike removal among attached cells (NaN cells skipped)
        for j in range(n - 2, -1, -1):
            if np.isnan(row[j]):
                continue
            k = j + 1
            while k < n and np.isnan(row[k]):
                k += 1
            if k < n and row[j] > row[k]:
                row[j] = row[k]

        # (2) detached fill: max(Newtonian, last attached value)
        run_max = -np.inf
        for j in range(n):
            if np.isnan(row[j]):
                fill = cp_newtonian(M, CONE_THETAS[j], gamma)
                row[j] = fill if not np.isfinite(run_max) else max(fill, run_max)
            else:
                run_max = max(run_max, row[j])
        cleaned[i] = row
    return cleaned


def load_cone_table():
    """Load the cached cone table; build + save it on first use (one-time, ~3-6 min).

    The cache is a human-readable CSV (rows = Mach, columns = theta-deg, with
    header labels) so it can be inspected or edited in any spreadsheet. Mach and
    theta grids still come from config; only the Cp grid is stored. Resolved
    against the PARENT package dir (capsule_opt/) so a single copy is reused."""

    path = config.CONE_CACHE

    if os.path.exists(path):
        table = read_cone_csv(path)
        if table.shape == (CONE_MACHS.size, CONE_THETAS.size):
            return table
        print(f"Cahched cone table at {path} has the wrong shape"
              f"{table.shape} (expected {(CONE_MACHS.size, CONE_THETAS.size)}); rebuilding.")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    print(f"Building tangent-cone table (one-time, ~3-6 min) -> {path}")
    table = clean_cone_table(build_cone_table())
    write_cone_csv(path, table)
    return table


def read_cone_csv(path):
    """Read a labeled cone-table CSV -> (nMach, nTheta) Cp array."""
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return data[:, 1:]                        # drop the Mach-label first column


def write_cone_csv(path, table):
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
