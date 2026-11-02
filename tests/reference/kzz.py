"""Mixing length Kzz, a Python port of PICASO's ``get_kzz``.

``test_rce.py`` checks ``rce_kzz`` in ``carma_rce.F90`` against this module, and
``test_kzz_picaso.py`` checks this module against PICASO.

Units are cgs (pressure in barye, Kzz in cm^2/s) except ``fnet_levels``, which
is W/m^2 as the Fortran and ``Results.flux_net`` carry it.
"""

import numpy as np
from numpy.typing import ArrayLike

from carmapy.adiabat import ADIABAT_PARMENTIER, ADIABAT_TABLE, grad_ad
from carmapy.constants import PARMENTIER_A_COEFF, PARMENTIER_B_COEFF

__all__ = ["mixing_length_kzz", "mixing_length_kzz_layers"]

#: Scale factor on the mixing length velocity.
KZZ_SCALEF = 1.0 / 3.0
#: Floor on the mixing length as a fraction of the scale height.  This is what
#: keeps Kzz finite in a stably stratified layer.
KZZ_MIXL_MIN = 0.1
#: How fast the convective flux may fall with height, as a fraction of the
#: pressure ratio across a pair of layers.  Represents convective overshoot;
#: the value is arbitrary.
KZZ_CHF_FALL = 1.0 / 3.0
#: Floor on the convective flux, as a fraction of the effective temperature.
#: Enters as the fourth power.
KZZ_TEFF_MIN = 0.05
#: Specific heat used by the parameterisation, as a multiple of the specific
#: gas constant (diatomic ideal gas).  Independent of the run's own Cp, as in
#: PICASO.
KZZ_CP_OVER_R = 7.0 / 2.0

#: Universal gas constant [erg/mol/K], written as PICASO writes it so the two
#: agree to the digit.
RGAS_CGS = 8.3143e7
#: Stefan-Boltzmann constant [erg/cm^2/s/K^4], likewise.
SIGMA_SB_CGS = 0.56687e-4
#: W/m^2 -> erg/cm^2/s
FLUX_SI2CGS = 1e3

#: Floor on the adiabatic gradient, matching ``GRAD_AD_MIN`` in
#: ``carma_rce.F90``.  The Parmentier fit crosses zero at 9600 K, outside the
#: domain it was fitted on; the tabulated gradient is positive everywhere, so
#: this binds only on the fit.
GRAD_AD_MIN = 1e-3


def _grad_ad(T, P, adiabat):
    """The adiabatic gradient of the named model, floored, as in ``rce_grad``."""
    if adiabat == ADIABAT_TABLE:
        grad = grad_ad(T, P)
    elif adiabat == ADIABAT_PARMENTIER:
        grad = PARMENTIER_A_COEFF - PARMENTIER_B_COEFF * np.asarray(T,
                                                                    dtype=float)
    else:
        raise ValueError(f"unknown adiabat {adiabat!r}, expected one of "
                         f"{[ADIABAT_PARMENTIER, ADIABAT_TABLE]}")

    return np.maximum(grad, GRAD_AD_MIN)


def _lapse_rate(T, P):
    """``d ln T / d ln P`` of a profile, on its own grid.

    Centred where there are neighbours on both sides, one-sided at the ends.
    """
    log_t = np.log(np.asarray(T, dtype=float))
    log_p = np.log(np.asarray(P, dtype=float))

    if log_t.size == 1:
        return np.zeros(1)

    dtdp = np.empty_like(log_t)
    dtdp[1:-1] = ((log_t[2:] - log_t[:-2]) / (log_p[2:] - log_p[:-2]))
    dtdp[0] = (log_t[1] - log_t[0]) / (log_p[1] - log_p[0])
    dtdp[-1] = (log_t[-1] - log_t[-2]) / (log_p[-1] - log_p[-2])

    return dtdp


def mixing_length_kzz(P_centers: ArrayLike,
                      T_centers: ArrayLike,
                      P_levels: ArrayLike,
                      fnet_levels: ArrayLike,
                      t_int: float,
                      wt_mol: float,
                      surface_grav: float,
                      adiabat: str = ADIABAT_PARMENTIER,
                      mixl_scale: float = 1.0,
                      return_rescale: bool = False):
    """Eddy diffusion implied by a radiative-convective profile.

    The mixing length velocity is set by the convective heat flux, and the
    mixing length itself by how close the layer's lapse rate sits to the
    adiabatic one::

        Kzz = scalef * H * (l/H)**(4/3) * (R * F_conv / (rho * Cp))**(1/3)

    with ``l = mixl_scale * max(0.1, min(1, dlnT/dlnP / grad_ad)) * H``.

    Evaluated in every layer; the floor on ``l/H`` leaves stable layers a small
    residual Kzz.

    The convective flux is the emergent flux less the layer's net radiative
    flux, then limited by an overshoot floor on how fast it may fall with
    height, rescaled so the internal flux passes through the base, and floored
    at ``sigma*(0.05*t_int)**4``.

    Parameters
    ----------
    P_centers : ArrayLike
        Layer pressures [barye], bottom to top, ``(NZ,)``
    T_centers : ArrayLike
        Layer temperatures [K], ``(NZ,)``
    P_levels : ArrayLike
        Level pressures [barye], ``(NZ+1,)``
    fnet_levels : ArrayLike
        Net *upward* flux at levels [W/m^2], ``(NZ+1,)``.  This is what
        ``Results.flux_net`` holds.
    t_int : float
        Internal temperature [K]
    wt_mol : float
        Mean molecular weight [g/mol]
    surface_grav : float
        Surface gravity [cm/s^2]
    adiabat : str, optional
        Which adiabat the lapse ratio is measured against.  See
        ``carmapy.adiabat``.  By default the Parmentier fit.
    mixl_scale : float, optional
        Multiplier on the mixing length ``l``, applied after the ``0.1*H``
        floor.  Kzz goes as ``l**(4/3)``, so a scale of ``f`` moves Kzz by
        ``f**(4/3)``.  Must be positive.  By default 1, as in PICASO.
    return_rescale : bool, optional
        Also return the factor the convective flux was scaled by at the base
        of the column, ``sigma*t_int**4 / F_TOA`` (1 in radiative-convective
        equilibrium).  By default False.

    Returns
    -------
    np.ndarray
        ``(NZ+1,)`` eddy diffusion at levels [cm^2/s], ready for
        ``Carma.add_kzz``.  With ``return_rescale``, a ``(kzz, rescale)``
        tuple.

    Notes
    -----
    With an incident beam this departs from PICASO, which uses the infrared
    net flux alone; ``fnet`` here is infrared plus shortwave.
    """
    p = np.asarray(P_centers, dtype=float)
    pl = np.asarray(P_levels, dtype=float)

    kz_lay, rescale = mixing_length_kzz_layers(
        P_centers, T_centers, P_levels, fnet_levels,
        t_int=t_int, wt_mol=wt_mol, surface_grav=surface_grav,
        adiabat=adiabat, mixl_scale=mixl_scale)

    # Linear in (ln Kzz, ln P), held at the end layers' values on the two
    # boundary levels.
    kzz = np.empty(p.size + 1)
    kzz[0] = kz_lay[0]
    kzz[-1] = kz_lay[-1]

    w = np.log(pl[1:-1] / p[:-1]) / np.log(p[1:] / p[:-1])
    kzz[1:-1] = np.exp((1 - w) * np.log(kz_lay[:-1]) + w * np.log(kz_lay[1:]))

    if return_rescale:
        return kzz, rescale

    return kzz


def mixing_length_kzz_layers(P_centers: ArrayLike,
                             T_centers: ArrayLike,
                             P_levels: ArrayLike,
                             fnet_levels: ArrayLike,
                             t_int: float,
                             wt_mol: float,
                             surface_grav: float,
                             adiabat: str = ADIABAT_PARMENTIER,
                             mixl_scale: float = 1.0):
    """The mixing length Kzz on the layers, as PICASO's ``get_kzz`` returns it.

    :func:`mixing_length_kzz` is this plus interpolation onto levels.
    Arguments are :func:`mixing_length_kzz`'s.

    Returns
    -------
    tuple[np.ndarray, float]
        ``(NZ,)`` eddy diffusion at layer centres [cm^2/s], and the factor the
        convective flux was rescaled by at the base of the column.
    """
    p = np.asarray(P_centers, dtype=float)
    t = np.asarray(T_centers, dtype=float)
    pl = np.asarray(P_levels, dtype=float)
    fnet = np.asarray(fnet_levels, dtype=float) * FLUX_SI2CGS

    if mixl_scale <= 0:
        raise ValueError(f"mixl_scale must be positive, got {mixl_scale}")

    nz = p.size
    if t.size != nz:
        raise ValueError(f"P_centers has {nz} entries but T_centers has "
                         f"{t.size}")
    if pl.size != nz + 1 or fnet.size != nz + 1:
        raise ValueError(f"P_levels and fnet_levels must have {nz + 1} "
                         f"entries, got {pl.size} and {fnet.size}")

    r_atmos = RGAS_CGS / wt_mol
    cp = KZZ_CP_OVER_R * r_atmos

    # Emergent flux against the target; their ratio rescales the convective
    # flux below.
    f_sum = fnet[-1]
    f_target = SIGMA_SB_CGS * t_int**4
    flx_min = SIGMA_SB_CGS * (KZZ_TEFF_MIN * t_int) ** 4

    # ---- convective flux --------------------------------------------------
    # What the radiation is not carrying, with fnet averaged onto layers.  The
    # deepest layer carries the whole flux convectively.
    chf = f_sum - 0.5 * (fnet[:-1] + fnet[1:])
    chf[0] = f_sum

    # Convective overshoot: the flux may not fall faster with height than a
    # fixed fraction of the pressure ratio.  Swept upward, since each layer is
    # floored against the one below.
    for iz in range(1, nz):
        chf[iz] = max(chf[iz], KZZ_CHF_FALL * (p[iz] / p[iz - 1]) * chf[iz - 1])

    # Put the imposed internal flux through the base of the column, correcting
    # the layers above it by the same factor.
    rescale = f_target / chf[0] if chf[0] > 0 else 1.0
    chf = np.maximum(chf * rescale, flx_min)

    # ---- mixing length ----------------------------------------------------
    lapse = np.minimum(1.0, _lapse_rate(t, p) / _grad_ad(t, p, adiabat))

    scale_h = r_atmos * t / surface_grav
    # The scale multiplies the mixing length after the floor.
    mixl = mixl_scale * np.maximum(KZZ_MIXL_MIN, lapse) * scale_h
    rho = p / (r_atmos * t)

    kz_lay = (KZZ_SCALEF * scale_h * (mixl / scale_h) ** (4 / 3)
              * (r_atmos * chf / (rho * cp)) ** (1 / 3))

    return kz_lay, rescale
