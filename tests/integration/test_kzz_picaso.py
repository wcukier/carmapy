"""``tests.reference.kzz`` against ``picaso.climate.get_kzz``.

``test_rce.py`` checks the Fortran against the reference module; this checks
the reference module against PICASO.  PICASO's reference data is not needed:
the gradient table ``get_kzz`` interpolates is the one CARMApy bundles as
``data/adiabat_grad.txt``.  Skipped when picaso is not installed.
"""
from collections import namedtuple

import numpy as np
import pytest

pytestmark = pytest.mark.integration

SIGMA_SB_CGS = 0.56687e-4


@pytest.fixture(scope="module")
def picaso_get_kzz():
    """PICASO's ``get_kzz`` without the numba jit.

    numba type-infers the ``moist=True`` branch, which reads condensate fields
    a dry ``Atmosphere`` does not carry; ``.py_func`` skips compilation.
    """
    climate = pytest.importorskip("picaso.climate",
                                  reason="picaso not installed")
    return climate.get_kzz.py_func


def _adiabat_bundle():
    """PICASO's ``AdiabatBundle``, built from the table CARMApy bundles.

    ``did_grad_cp`` hard-codes a 53x26 grid, so the shape is asserted.
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
    """A column built from levels, with layers derived as PICASO derives them
    (geometric-mean pressure, arithmetic-mean temperature)."""
    from carmapy.adiabat import tabulated_adiabat

    # Bottom to top, inside the gradient table's range.
    pl_bar = np.geomspace(3e2, 1e-2, nz + 1)[::-1]

    # An adiabatic interior under a shallower radiative gradient, joined at
    # 5 bar, so both the ``min(1, ...)`` and ``max(0.1, ...)`` clips bind.
    p_rcb = 5.0
    tl = tabulated_adiabat(pl_bar * 1e6, 2200.0, 3e8)
    above = pl_bar < p_rcb
    t_rcb = np.interp(np.log(p_rcb), np.log(pl_bar), tl)
    tl[above] = t_rcb * (pl_bar[above] / p_rcb) ** 0.08

    p_bar = np.sqrt(pl_bar[1:] * pl_bar[:-1])
    t = 0.5 * (tl[1:] + tl[:-1])

    # Net upward flux at levels [W/m^2]: zero at the base, rising to a little
    # over the internal flux at the top.  The bump partway up carries more than
    # the column emits, so the convective flux beneath it goes negative and the
    # overshoot sweep and absolute floor both bind.
    f_int = SIGMA_SB_CGS * t_int**4 * 1e-3
    x = np.linspace(0.0, 1.0, nz + 1)
    fnet = 1.05 * f_int * x ** 0.6 * (1 + 0.9 * np.exp(-((x - 0.72) / 0.1)**2))

    return dict(nz=nz, t_int=t_int, wt_mol=wt_mol, grav=grav,
                pl_bar=pl_bar, tl=tl, p_bar=p_bar, t=t, fnet=fnet)


def _picaso_kzz(get_kzz, case, dtdp):
    """``get_kzz`` on ``case``, returned bottom-to-top.

    PICASO's arrays run top-to-bottom.  ``dtdp`` is passed in because PICASO
    takes it from its climate solver; using CARMApy's lapse rate isolates the
    parameterisation.
    """
    nz = case["nz"]

    Atmosphere = namedtuple("Atmosphere",
                            ["p_level", "t_level", "mmw_layer", "dtdp"])
    atm = Atmosphere(p_level=case["pl_bar"][::-1],
                     t_level=case["tl"][::-1],
                     mmw_layer=np.full(nz, case["wt_mol"]),
                     dtdp=dtdp[::-1])

    # sigma*T_int^4 in cgs; PICASO reads abs(tidal[0]).
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
    """``mixing_length_kzz_layers`` against ``picaso.climate.get_kzz``, to
    round-off."""
    from tests.reference.kzz import _lapse_rate, mixing_length_kzz_layers

    case = _column()

    got, rescale = mixing_length_kzz_layers(
        case["p_bar"] * 1e6, case["t"], case["pl_bar"] * 1e6, case["fnet"],
        t_int=case["t_int"], wt_mol=case["wt_mol"],
        surface_grav=case["grav"], adiabat="table")

    ref = _picaso_kzz(picaso_get_kzz, case,
                      _lapse_rate(case["t"], case["p_bar"]))

    np.testing.assert_allclose(got, ref, rtol=1e-10)

    # The column radiates more than its target, so the rescale is not 1.
    f_target = SIGMA_SB_CGS * case["t_int"] ** 4
    f_sum = case["fnet"][-1] * 1e3
    assert rescale == pytest.approx(f_target / f_sum, rel=1e-12)
    assert abs(rescale - 1.0) > 0.02
