"""Trim angle-of-attack solver (zero pitching moment) via bisection."""

import numpy as np
from . import config


def Cm_at(M_clip: float, aoa: float, interp_CM) -> float:
    """Pitching-moment coefficient at "(M, aoa)" from the aero-database interpolator."""
    return float(interp_CM(np.array([[M_clip, aoa]]))[0])


def trim_alpha(M: float, interp_CM,
               aoa_min: float = -5,
               aoa_max: float = 35,
               n_braket: int = 200) -> float | None:
    """Find the trim angle of attack where the pitching moment crosses zero.

    Brackets a sign change of "Cm" on a coarse grid, then refines with
    bisection, keeping the root where "Cm"  is locally decreasing (stable).

    Parameters
    ----------
    M : float
        Mach number (clipped to the aero-database range).
    interp_CM : callable
        Pitching-moment interpolator from
        :func:`capsule_opt.aero_database.build_aero_database`.
    aoa_min, aoa_max : float, optional
        Angle-of-attack search bounds [deg].
    n_braket : int, optional
        Number of grid points used to bracket the root.

    Returns
    -------
    float or None
        Trim angle of attack [deg], or "None" if no stable root is found.
    """
    M_c = float(np.clip(M, config.mach_nodes[0], config.mach_nodes[-1]))

    aoas = np.linspace(aoa_min, aoa_max, n_braket)
    Cms = np.array([Cm_at(M_c, a, interp_CM) for a in aoas])

    result = None
    for idx in np.where(np.diff(np.sign(Cms)))[0]:
        Co, C1 = Cms[idx], Cms[idx + 1]
        if C1 < Co:
            ao, a1 = aoas[idx], aoas[idx + 1]

            for _ in range(10):  # bisection refinement
                amid = 0.5 * (ao + a1)
                cmid = Cm_at(M_c, amid, interp_CM)
                if cmid * Co < 0:
                    a1, C1 = amid, cmid
                else:
                    ao, Co = amid, cmid

            result = float(0.5 * (ao + a1))
            break
    return result
