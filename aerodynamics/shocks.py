from .. import config
import numpy as np
from scipy.optimize import brentq

def shock_post(M, beta, gamma=config.gamma):
    """State immediately behind an oblique shock at angle beta."""

    Mn1   = M * np.sin(beta)
    Mn2sq = (1 + (gamma - 1) / 2 * Mn1 ** 2) / (gamma * Mn1 ** 2 - (gamma - 1) / 2)
    delta = np.arctan(2 * (1 / np.tan(beta)) * (M ** 2 * np.sin(beta) ** 2 - 1)
                      / (M ** 2 * (gamma + np.cos(2 * beta)) + 2))
    M2    = np.sqrt(Mn2sq) / np.sin(beta - delta)

    return Mn1, delta, M2 / np.sqrt(1 + (gamma - 1) / 2 * M2 ** 2)


def tm_rhs(psi, u, v, gamma=config.gamma):
    """Taylor-Maccoll ODE right-hand side (velocity normalized by a0)."""

    A = 1 - (gamma - 1) / 2 * (u * u + v * v)
    return v, (u*v*v - 2*u*A - v*A/np.tan(psi)) / (A - v * v)


def tm_integrate(u0, v0, beta, thc, gamma=config.gamma, n=200):
    """Fixed-step RK4 of Taylor-Maccoll from shock (psi=beta) down to surface (psi=thc).
       Returns (u, v) at the surface, or None if the sonic singularity is hit.
       np.errstate silences the expected divide-by-zero at the singularity."""
    
    h        = (thc - beta) / n                       # thc < beta -> h negative
    psi, u, v = beta, u0, v0
    with np.errstate(divide='ignore', invalid='ignore'):
        for _ in range(n):
            k1u, k1v = tm_rhs(psi, u, v, gamma)
            k2u, k2v = tm_rhs(psi + h / 2, u + h * k1u / 2, v + h * k1v / 2, gamma)
            k3u, k3v = tm_rhs(psi + h / 2, u + h * k2u / 2, v + h * k2v / 2, gamma)
            k4u, k4v = tm_rhs(psi + h, u + h * k3u, v + h * k3v, gamma)

            u   += (h / 6) * (k1u + 2 * k2u + 2 * k3u + k4u)
            v   += (h / 6) * (k1v + 2 * k2v + 2 * k3v + k4v)
            psi += h

            if not (np.isfinite(u) and np.isfinite(v)):
                return None
            
        return u, v


def cp_tangent_cone(M, theta_deg, gamma=config.gamma):
    """Tangent-cone Cp (conical shock).
       theta_deg = local windward surface inclination (deg >= 0).
       Returns None if the shock is detached (no attached solution)."""
    
    thc = np.deg2rad(theta_deg)
    mu  = np.arcsin(min(1 / M, 1.0))

    def resid(beta):
        if beta <= thc:                        # shock must stand outside the cone
            return np.nan
        Mn1, delta, V2 = shock_post(M, beta, gamma)
        if Mn1 <= 1:
            return np.nan
        s = tm_integrate(V2 * np.cos(beta - delta),
                         -V2 * np.sin(beta - delta),
                         beta, thc, gamma)
        return s[1] if s is not None else np.nan

    bs = np.linspace(max(mu, thc) + 1e-3, np.pi / 2 - 1e-3, 90)
    r  = np.array([resid(b) for b in bs])

    # weakest-shock (smallest-beta) sign change among finite adjacent samples
    pairs = [(bs[k], bs[k + 1]) for k in range(len(bs) - 1)
             if np.isfinite(r[k]) and np.isfinite(r[k + 1]) and r[k] * r[k + 1] < 0]
    if not pairs:
        return None

    def safe_root(a, b):
        # shrink [a,b] to a NaN-free sign-changing sub-interval, then brentq
        sub = np.linspace(a, b, 41)
        rs  = np.array([resid(x) for x in sub])
        clean = [(sub[k], sub[k + 1]) for k in range(len(sub) - 1)
                 if np.isfinite(rs[k]) and np.isfinite(rs[k + 1]) and rs[k] * rs[k + 1] < 0]
        if not clean:
            return None
        a, b = clean[0]
        return brentq(resid, a, b, xtol=1e-9, maxiter=50)

    try:
        beta = brentq(resid, pairs[0][0], pairs[0][1], xtol=1e-9, maxiter=50)
    except ValueError:
        # bracket held a NaN notch (erratic RK4 near the Mach-angle limit)
        beta = safe_root(pairs[0][0], pairs[0][1])
        if beta is None:
            return None

    Mn1, delta, V2 = shock_post(M, beta, gamma)
    # Final high-resolution surface integration. tm_integrate returns None if the
    # RK4 hits the sonic singularity before reaching the surface -- a marginal /
    # near-detached root. Treat that as "no attached solution" (contract: None),
    # matching the other None-returning paths instead of crashing on unpack.
    surf = tm_integrate(V2 * np.cos(beta - delta),
                        -V2 * np.sin(beta - delta),
                        beta, thc, gamma, n=300)
    if surf is None:
        return None
    u_s, _ = surf

    As  = 1 - (gamma - 1) / 2 * u_s ** 2
    Ms2 = u_s ** 2 / As

    p_s_p02 = (1 + (gamma - 1) / 2 * Ms2) ** (-gamma / (gamma - 1))
    p02_p01 = ((((gamma + 1) * Mn1 ** 2) / ((gamma - 1) * Mn1 ** 2 + 2)) ** (gamma / (gamma - 1))
               * ((gamma + 1) / (2 * gamma * Mn1 ** 2 - (gamma - 1))) ** (1 / (gamma - 1)))
    p01_pin = (1 + (gamma - 1) / 2 * M ** 2) ** (gamma / (gamma - 1))

    return (2.0 / (gamma * M ** 2)) * (p_s_p02 * p02_p01 * p01_pin - 1)