"""Adiabatic temperature profiles for the deep atmosphere.

Two interchangeable models, both mapping a pressure grid onto temperatures
given a single anchor point:

``parmentier_adiabat``
    The Parmentier et al. (2015) analytic fit, whose gradient is
    ``d ln T / d ln P = a - b*T``.  Cheap, closed form, and accurate to better
    than 2% below ~1500 K.

``tabulated_adiabat``
    Integration of the tabulated ``\\nabla_ad(T, P)`` for an H/He mixture (see
    ``data/README.md``).  Above ~2000 K, where H2 dissociation flattens the
    real adiabat to less than half the analytic value, this is the only one of
    the two that is right.

They take the same arguments and return the same shape, so either can be
dropped in for the other.  ``Carma.adiabat`` selects between them by name.
"""

import os
import warnings
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp

from .constants import (BAR_TO_BARYE, PARMENTIER_A_COEFF, PARMENTIER_B_COEFF)

_DATA = os.path.join(os.path.dirname(__file__), "data")

#: The bundled gradient table. Passed to the Fortran engine by path rather
#: than copied into each run directory.
TABLE_FILE = os.path.join(_DATA, "adiabat_grad.txt")

ADIABAT_PARMENTIER = "parmentier"
ADIABAT_TABLE = "table"


@lru_cache(maxsize=1)
def _load_grad_table():
    """Read the bundled H/He adiabatic gradient table.

    Returns
    -------
    tuple
        ``(log10 T [K], log10 P [bar], grad)`` with shapes ``(nt,)``,
        ``(np,)`` and ``(nt, np)``.  The specific heat block that follows the
        gradient in the file is not returned; nothing needs it yet.
    """
    with open(TABLE_FILE) as f:
        rows = [line for line in f if not line.lstrip().startswith("#")]

    nt, npres = (int(tok) for tok in rows[0].split())
    values = np.array([float(tok) for row in rows[1:] for tok in row.split()])

    expected = nt + npres + 2 * nt * npres
    if values.size != expected:
        raise ValueError(f"{TABLE_FILE} holds {values.size} values, expected "
                         f"{expected} for a {nt}x{npres} table")

    t_axis = values[:nt]
    p_axis = values[nt:nt + npres]
    grad = values[nt + npres:nt + npres + nt * npres].reshape(nt, npres)

    return t_axis, p_axis, grad


def _interpolate_grad(log_t, log_p, excursion=None):
    """Bilinearly interpolate the gradient table, clamping at its edges.

    Clamping rather than extrapolating matches how PICASO reads the same
    table.  ``excursion`` is an optional dict that accumulates how far outside
    the grid the requests reached, so a caller driving this in a loop can warn
    once at the end instead of on every evaluation.
    """
    t_axis, p_axis, grad = _load_grad_table()

    if excursion is not None:
        for name, value, axis in (("temperature", log_t, t_axis),
                                  ("pressure", log_p, p_axis)):
            lo, hi = excursion.get(name, (0.0, 0.0))
            excursion[name] = (min(lo, float(np.min(value)) - axis[0]),
                               max(hi, float(np.max(value)) - axis[-1]))

    log_t = np.clip(log_t, t_axis[0], t_axis[-1])
    log_p = np.clip(log_p, p_axis[0], p_axis[-1])

    it = np.clip(np.searchsorted(t_axis, log_t) - 1, 0, t_axis.size - 2)
    ip = np.clip(np.searchsorted(p_axis, log_p) - 1, 0, p_axis.size - 2)

    ft = (log_t - t_axis[it]) / (t_axis[it + 1] - t_axis[it])
    fp = (log_p - p_axis[ip]) / (p_axis[ip + 1] - p_axis[ip])

    return ((1 - ft) * (1 - fp) * grad[it, ip]
            + ft * (1 - fp) * grad[it + 1, ip]
            + ft * fp * grad[it + 1, ip + 1]
            + (1 - ft) * fp * grad[it, ip + 1])


def _warn_excursion(excursion):
    """Emit one warning describing how far a request left the table."""
    t_axis, p_axis, _ = _load_grad_table()

    messages = []
    for name, unit, axis in (("temperature", "K", t_axis),
                             ("pressure", "bar", p_axis)):
        below, above = excursion.get(name, (0.0, 0.0))
        if below < 0:
            messages.append(f"{name} down to {10**(axis[0] + below):.3g} {unit}"
                            f" (table starts at {10**axis[0]:.3g} {unit})")
        if above > 0:
            messages.append(f"{name} up to {10**(axis[-1] + above):.3g} {unit}"
                            f" (table ends at {10**axis[-1]:.3g} {unit})")

    if messages:
        warnings.warn("adiabat gradient table clamped at its edges: "
                      + "; ".join(messages), UserWarning, stacklevel=3)


def grad_ad(T, P):
    """Adiabatic gradient ``d ln T / d ln P`` of an H/He mixture.

    Bilinear in ``(log10 T, log10 P)`` over the bundled table, clamped at the
    table edges.  Useful on its own for testing a profile's actual gradient
    against the adiabatic one.

    Parameters
    ----------
    T : ArrayLike
        Temperature [K]
    P : ArrayLike
        Pressure [barye]

    Returns
    -------
    np.ndarray
        The gradient, dimensionless, broadcast to the shape of ``T`` and ``P``
    """
    log_t = np.log10(np.asarray(T, dtype=float))
    log_p = np.log10(np.asarray(P, dtype=float) / BAR_TO_BARYE)

    excursion = {}
    out = _interpolate_grad(log_t, log_p, excursion)
    _warn_excursion(excursion)

    return out


def parmentier_adiabat(P, t0, p0):
    """Temperature along a Parmentier et al. (2015) adiabat.

    A fit to the Saumon et al. (1995) equation of state, with a gradient of
    ``d ln T / d ln P = PARMENTIER_A_COEFF - PARMENTIER_B_COEFF * T``.

    Parameters
    ----------
    P : ArrayLike
        Pressures to evaluate at [barye], in any order
    t0 : float
        Anchor temperature [K]
    p0 : float
        Anchor pressure [barye]

    Returns
    -------
    np.ndarray
        Temperature [K] at each pressure, shaped like ``P``

    References
    ----------
    .. [1] Parmentier, V., Guillot, T., Fortney, J. J., & Marley,
       M. S. 2015, A&A, 574, A35
    """
    P = np.asarray(P, dtype=float)

    kappa = (t0 / (PARMENTIER_A_COEFF - PARMENTIER_B_COEFF * t0)
             * (P / p0) ** PARMENTIER_A_COEFF)

    return PARMENTIER_A_COEFF * kappa / (1 + PARMENTIER_B_COEFF * kappa)


def tabulated_adiabat(P, t0, p0):
    """Temperature along the tabulated H/He adiabat.

    Integrates ``d ln T / d ln P = grad_ad(T, P)`` outward from the anchor.
    Unlike the analytic fit this has no closed form, so the whole pressure grid
    is integrated in one sweep -- pass all the pressures at once rather than
    calling this per level.

    Parameters
    ----------
    P : ArrayLike
        Pressures to evaluate at [barye], in any order
    t0 : float
        Anchor temperature [K]
    p0 : float
        Anchor pressure [barye]

    Returns
    -------
    np.ndarray
        Temperature [K] at each pressure, shaped like ``P``

    Notes
    -----
    The table covers ``T`` in [10, 3981] K and ``P`` in [1e-2, 1e3] bar.
    Beyond it the gradient is held at its edge value and a warning is raised;
    a profile that leaves the grid is being extrapolated, not tabulated.
    """
    P = np.asarray(P, dtype=float)

    x = np.log(np.atleast_1d(P).ravel())
    x0 = np.log(p0)

    # Points at the anchor keep t0 exactly; the rest are integrated away from
    # it in both directions, each branch starting from the anchor so the
    # integration error never accumulates across it.
    log_t = np.full(x.shape, np.log(t0))
    excursion = {}

    def gradient(xi, y):
        return _interpolate_grad(y / np.log(10),
                                 (xi - np.log(BAR_TO_BARYE)) / np.log(10),
                                 excursion)

    for direction in (1, -1):
        branch = np.nonzero(direction * (x - x0) > 0)[0]
        if branch.size == 0:
            continue

        branch = branch[np.argsort(direction * x[branch])]

        # rtol=1e-8 leaves a few parts in 1e6 of integration error, which is
        # enough to show up against the Fortran port. The whole integration
        # costs milliseconds, so there is no reason not to ask for more.
        solution = solve_ivp(gradient, (x0, x[branch[-1]]), [np.log(t0)],
                             t_eval=x[branch], rtol=1e-10, atol=1e-12)
        if not solution.success:
            raise RuntimeError("adiabat integration failed between "
                               f"{np.exp(x0):.4g} and "
                               f"{np.exp(x[branch[-1]]):.4g} barye: "
                               f"{solution.message}")

        log_t[branch] = solution.y[0]

    _warn_excursion(excursion)

    # exp(log(t0)) is not always t0, and the anchor is worth having exact.
    temperature = np.exp(log_t)
    temperature[x == x0] = t0

    return temperature.reshape(P.shape)


#: The selectable adiabat models, by the name ``Carma.adiabat`` takes.
ADIABAT_MODELS = {ADIABAT_PARMENTIER: parmentier_adiabat,
                  ADIABAT_TABLE: tabulated_adiabat}

#: How each model is encoded in the Fortran namelist.
ADIABAT_CODES = {ADIABAT_PARMENTIER: 0, ADIABAT_TABLE: 1}


def get_adiabat(name):
    """Look up an adiabat model by name.

    Parameters
    ----------
    name : str
        One of the keys of ``ADIABAT_MODELS``

    Returns
    -------
    callable
        ``f(P, t0, p0) -> np.ndarray``
    """
    try:
        return ADIABAT_MODELS[name]
    except KeyError:
        raise ValueError(f"unknown adiabat {name!r}, expected one of "
                         f"{sorted(ADIABAT_MODELS)}") from None
