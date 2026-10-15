"""Heliocentric N-body integration: Sun + Earth(-Moon barycentre) + Jupiter +
Saturn + massless test particle(s) (e.g. Ceres), for the perturbed-asteroid
epicycle study.

Method: *barycentric* Newtonian point-mass integration of all bodies including
the Sun as a moving body (no indirect terms needed); heliocentric output is
obtained by subtracting the Sun's barycentric position.  Units: AU, day,
masses in solar masses (GM_sun = k^2 = kepler.MU_SUN, so G*m_i = MU_SUN * m_i).
Bodies with mass 0 are test particles: they feel all massive bodies and exert
no force.  Integrator: scipy DOP853, rtol=1e-12, atol=1e-14 (AU, AU/day).

Masses (Sun/planet mass ratios, IAU 2012 / JPL DE430 values, Folkner et al.
2014, IPN PR 42-196; Standish 1998 for the planets): Earth+Moon 328900.5614,
Jupiter system 1047.3486 (planet+satellites), Saturn system 3497.898.  Ceres
(GM ~ 62.6 km^3/s^2 ~ 4.7e-10 Msun) is treated as massless by default.

Initial conditions: each body's osculating Keplerian elements (bodies.py; the
planet sets are JPL "approximate positions" at J2000 with mu = MU_SUN) are
propagated with Kepler's equation to the common start epoch t0 and converted to
heliocentric position/velocity; the Sun is then placed so that the barycentre
is at rest at the origin.  Note the heliocentric velocities come from two-body
Kepler with mu = k^2 (not k^2 (1 + m)), so the planets' osculating orbits at t0
differ from the element sets at the 1e-3 (Jupiter) level in mean motion;
irrelevant for this study, which compares against the same integration.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .kepler import MU_SUN, Elements, position, velocity

#: mass / M_sun for the planets (see module docstring)
MASS_RATIOS = {
    "Earth": 1.0 / 328900.5614,
    "Jupiter": 1.0 / 1047.3486,
    "Saturn": 1.0 / 3497.898,
}


def _accel(pos: np.ndarray, m: np.ndarray) -> np.ndarray:
    """pos (n,3), m (n,) masses in Msun -> accelerations (n,3) [AU/day^2]."""
    d = pos[None, :, :] - pos[:, None, :]            # d[i,j] = r_j - r_i
    r2 = np.sum(d * d, axis=-1)
    np.fill_diagonal(r2, 1.0)
    w = m[None, :] / (r2 * np.sqrt(r2))
    np.fill_diagonal(w, 0.0)
    return MU_SUN * np.einsum("ij,ijk->ik", w, d)


def energy(pos: np.ndarray, vel: np.ndarray, m: np.ndarray) -> np.ndarray:
    """Total energy [AU^2/day^2 in units where G*M_sun = MU_SUN; per M_sun].

    pos, vel: (..., n, 3); m: (n,).  Returns (...,).  Massless bodies carry
    no energy here (they are not conserved subsystems by themselves).
    """
    kin = 0.5 * np.einsum("i,...ik,...ik->...", m, vel, vel)
    d = pos[..., :, None, :] - pos[..., None, :, :]
    r = np.sqrt(np.sum(d * d, axis=-1))
    n = len(m)
    iu = np.triu_indices(n, 1)
    pot = -MU_SUN * np.sum((m[:, None] * m[None, :])[iu] / r[..., iu[0], iu[1]], axis=-1)
    return kin + pot


def integrate_state(bodies: dict[str, Elements], t0: float, t_span_days: float,
                    dt_sample: float, masses: dict[str, float] | None = None,
                    planets_massless: bool = False, rtol: float = 1e-12,
                    atol: float = 1e-14):
    """Integrate; return (t, names, m, pos_bary, vel_bary).

    names[0] is "Sun".  pos_bary/vel_bary: (N, n, 3) barycentric.  t: (N,) JD,
    t = t0 + k*dt_sample up to and including t0 + t_span_days (when divisible).
    ``masses``: name -> mass ratio (default MASS_RATIOS, others 0);
    ``planets_massless=True`` sets all non-Sun masses to 0 (pure Kepler check).
    """
    mm = dict(MASS_RATIOS if masses is None else masses)
    names = ["Sun"] + list(bodies)
    m = np.array([1.0] + [0.0 if planets_massless else mm.get(k, 0.0) for k in bodies])
    rh = np.array([position(e, t0) for e in bodies.values()])
    vh = np.array([velocity(e, t0) for e in bodies.values()])
    msum = m.sum()
    rs = -(m[1:, None] * rh).sum(0) / msum
    vs = -(m[1:, None] * vh).sum(0) / msum
    r0 = np.vstack([rs, rh + rs])
    v0 = np.vstack([vs, vh + vs])
    n = len(names)

    def rhs(_t, y):
        p, v = y[:3 * n].reshape(n, 3), y[3 * n:].reshape(n, 3)
        return np.concatenate([v.ravel(), _accel(p, m).ravel()])

    nstep = int(np.floor(t_span_days / dt_sample + 1e-9))
    te = t0 + dt_sample * np.arange(nstep + 1)
    sol = solve_ivp(rhs, (t0, te[-1]), np.concatenate([r0.ravel(), v0.ravel()]),
                    method="DOP853", t_eval=te, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    Y = sol.y.T
    return te, names, m, Y[:, :3 * n].reshape(-1, n, 3), Y[:, 3 * n:].reshape(-1, n, 3)


def integrate(bodies: dict[str, Elements], t0: float, t_span_days: float, dt_sample: float,
              masses: dict[str, float] | None = None, planets_massless: bool = False,
              **kw) -> dict:
    """Heliocentric trajectories: dict with "t" (N,) JD and, per body name, an
    (N,3) heliocentric ecliptic J2000 position array [AU] ("Sun" excluded).
    See :func:`integrate_state` for arguments."""
    t, names, _m, p, _v = integrate_state(bodies, t0, t_span_days, dt_sample, masses,
                                          planets_massless, **kw)
    out = {"t": t}
    for i, k in enumerate(names):
        if i:
            out[k] = p[:, i] - p[:, 0]
    return out


def geocentric(traj: dict, body: str = "Ceres", earth: str = "Earth") -> np.ndarray:
    """Geocentric position (N,3) = r_body - r_earth from an :func:`integrate` result."""
    return traj[body] - traj[earth]
