"""Check ``carmapy.kzz`` against the PICASO routine it is a port of.

``tests/integration/test_rce.py`` pins the Fortran to ``carmapy.kzz``; this
pins ``carmapy.kzz`` to ``picaso.climate.get_kzz``.  Only the two together
support the claim that a CARMA column and a PICASO column with the same profile
mix the same way -- the first alone would be satisfied by two identical
mistakes.

PICASO's reference data is *not* needed.  ``picaso.climate`` imports without it,
and the adiabatic gradient table ``get_kzz`` interpolates is the same Saumon
table CARMApy bundles as ``data/adiabat_grad.txt``, so the bundle can be
assembled here.  Skipped when picaso is not installed.

Run with ``CARMAPY_RUN_INTEGRATION=1 pytest tests/integration/test_kzz_picaso.py``.
"""
from collections import namedtuple

import numpy as np
import pytest

pytestmark = pytest.mark.integration

SIGMA_SB_CGS = 0.56687e-4


@pytest.fixture(scope="module")
def picaso_get_kzz():
    """PICASO's ``get_kzz``, as the plain Python function behind the jit.

    ``get_kzz`` is ``@jit(nopython=True)``, and numba type-infers the
    ``moist=True`` branch even though nothing here takes it -- that branch
    reads condensate fields off ``Atmosphere`` which a dry column has no
    business carrying.  ``.py_func`` is the same source without the
    compilation, which is what is under test here anyway.
    """
    climate = pytest.importorskip("picaso.climate",
                                  reason="picaso not installed")
    return climate.get_kzz.py_func


def _adiabat_bundle():
    """PICASO's ``AdiabatBundle``, built from the table CARMApy bundles.

    ``did_grad_cp`` indexes the grid with hard-coded bounds for 53 temperatures
    and 26 pressures, which is the shape of this table; a differently shaped
    one would read out of bounds rather than fail cleanly, so the shape is
    asserted here.
    """
    from carmapy.adiabat import TABLE_FILE

    with open(TABLE_FILE) as f:
        rows = [line for line in f if not line.lstrip().startswith("#")]

    nt, npres = (int(tok) for tok in rows[0].split())
    assert (nt, npres) == (53, 26), "did_grad_cp assumes a 53x26 grid"

    values = np.array([float(tok) for row in rows[1:] for tok in row.split()])
    t_table = values[:nt]
    p_table = values[nt:nt + npres]
    grad = values[nt + npres:nt + npres + nt * npres].reshape(nt, npres)
    cp = values[nt + npres + nt * npres:].reshape(nt, npres)

    Bundle = namedtuple("AdiabatBundle", ["t_table", "p_table", "grad", "cp"])
    return Bundle(t_table, p_table, grad, cp)


def _column(nz=40, t_int=1100.0, wt_mol=2.2, grav=1e5):
    """A profile both sides can be handed.

    Built from *levels*, because PICASO derives its layer grid from them --
    layer pressures as the geometric mean and layer temperatures as the
    arithmetic one -- while CARMApy takes the layer grid directly.  Starting
    from levels and deriving the layers the same way is what makes the two
    evaluate the same column rather than two similar ones.
    """
    from carmapy.adiabat import tabulated_adiabat

    # Bottom to top, inside the gradient table's own range so neither side is
    # comparing its clamping behaviour rather than its interpolation.
    pl_bar = np.geomspace(3e2, 1e-2, nz + 1)[::-1]

    # An adiabatic interior under a shallower radiative gradient, joined at
    # 5 bar.  A real structure rather than a smooth curve, and the point of it:
    # the lapse ratio saturates at 1 below the join and falls away above it, so
    # the ``min(1, ...)`` and the ``max(0.1, ...)`` both do work.
    p_rcb = 5.0
    tl = tabulated_adiabat(pl_bar * 1e6, 2200.0, 3e8)
    above = pl_bar < p_rcb
    t_rcb = np.interp(np.log(p_rcb), np.log(pl_bar), tl)
    tl[above] = t_rcb * (pl_bar[above] / p_rcb) ** 0.08

    p_bar = np.sqrt(pl_bar[1:] * pl_bar[:-1])
    t = 0.5 * (tl[1:] + tl[:-1])

    # Net upward flux at levels [W/m^2]: nothing at the base, where convection
    # carries everything, rising to a little over the internal flux at the top.
    # The bump partway up is what a cloud deck does, and it is here on purpose
    # -- it carries more than the column emits, so the convective flux beneath
    # it goes negative.  That is what puts the overshoot sweep and the absolute
    # floor to work instead of leaving them decorative, and a negative flux
    # under a cube root is exactly the case the floors exist to survive.
    f_int = SIGMA_SB_CGS * t_int**4 * 1e-3
    x = np.linspace(0.0, 1.0, nz + 1)
    fnet = 1.05 * f_int * x ** 0.6 * (1 + 0.9 * np.exp(-((x - 0.72) / 0.1)**2))

    return dict(nz=nz, t_int=t_int, wt_mol=wt_mol, grav=grav,
                pl_bar=pl_bar, tl=tl, p_bar=p_bar, t=t, fnet=fnet)


def _picaso_kzz(get_kzz, case, dtdp):
    """``get_kzz`` on ``case``, returned bottom-to-top.

    PICASO's arrays run top-to-bottom, so everything is reversed on the way in
    and the answer reversed back on the way out.  ``dtdp`` is passed in rather
    than recomputed because PICASO takes it from its climate solver instead of
    differencing the profile; handing it the lapse rate CARMApy computed is
    what isolates the parameterisation from how the two get their gradients.
    """
    nz = case["nz"]

    Atmosphere = namedtuple("Atmosphere",
                            ["p_level", "t_level", "mmw_layer", "dtdp"])
    atm = Atmosphere(p_level=case["pl_bar"][::-1],
                     t_level=case["tl"][::-1],
                     mmw_layer=np.full(nz, case["wt_mol"]),
                     dtdp=dtdp[::-1])

    # sigma*T_int^4 in cgs, as the flux the column is asked to carry.  PICASO
    # reads it off tidal[0] and takes its absolute value.
    tidal = np.full(nz + 1, -SIGMA_SB_CGS * case["t_int"] ** 4)

    # The layer net flux PICASO differences against, and the emergent flux it
    # sums for the total.  Both cgs.
    fnet_cgs = case["fnet"] * 1e3
    flux_net_ir_layer = (0.5 * (fnet_cgs[:-1] + fnet_cgs[1:]))[::-1]
    flux_plus_ir_attop = np.array([fnet_cgs[-1]])

    kz = get_kzz(case["grav"] * 1e-2, tidal, flux_net_ir_layer,
                 flux_plus_ir_attop, _adiabat_bundle(),
                 np.zeros(6, dtype=np.int64), atm, False)

    # get_kzz pads its layer array to level length by repeating the last entry;
    # the layers are what the two share.
    return np.asarray(kz)[:nz][::-1]


def test_kzz_matches_picaso(picaso_get_kzz):
    """``mixing_length_kzz_layers`` against ``picaso.climate.get_kzz``.

    Both evaluate the same closed-form expression on the same column against
    the same gradient table, so the agreement should be at round-off.  What
    would show up here and nowhere else is a coefficient copied wrong, a floor
    applied in the wrong order, or the convective flux sweep run in the wrong
    direction -- the index order is reversed between the two, and getting that
    backwards leaves a profile that still looks entirely plausible.
    """
    from carmapy.kzz import _lapse_rate, mixing_length_kzz_layers

    case = _column()

    got, rescale = mixing_length_kzz_layers(
        case["p_bar"] * 1e6, case["t"], case["pl_bar"] * 1e6, case["fnet"],
        t_int=case["t_int"], wt_mol=case["wt_mol"],
        surface_grav=case["grav"], adiabat="table")

    ref = _picaso_kzz(picaso_get_kzz, case,
                      _lapse_rate(case["t"], case["p_bar"]))

    np.testing.assert_allclose(got, ref, rtol=1e-10)

    # The column radiates more than it was asked to, so the rescale is live
    # rather than an identity that would agree by default.  Both sides derive
    # it the same way, from the target flux over the emergent one.
    f_target = SIGMA_SB_CGS * case["t_int"] ** 4
    f_sum = case["fnet"][-1] * 1e3
    assert rescale == pytest.approx(f_target / f_sum, rel=1e-12)
    assert abs(rescale - 1.0) > 0.02
