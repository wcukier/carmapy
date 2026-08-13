"""The two adiabat models and the gradient table behind one of them.

``parmentier_adiabat`` is checked against the closed form it replaced, so the
refactor that moved it out of ``extend_atmosphere`` cannot change a number.
``tabulated_adiabat`` has no closed form to check against, so it is pinned
three ways instead: it must reproduce the table at its own nodes, agree with
the fit where the two physically should agree, and disagree with it in the
specific direction H2 dissociation implies.
"""
import warnings

import numpy as np
import pytest

from carmapy.adiabat import (ADIABAT_MODELS, ADIABAT_PARMENTIER, ADIABAT_TABLE,
                             get_adiabat, grad_ad, parmentier_adiabat,
                             tabulated_adiabat)
from carmapy.adiabat import _load_grad_table
from carmapy.constants import (BAR_TO_BARYE, PARMENTIER_A_COEFF,
                               PARMENTIER_B_COEFF)

pytestmark = pytest.mark.unit

# Pressures spanning the table, 1e-2 to 1e3 bar.
P_IN_TABLE = np.geomspace(1e4, 1e9, 25)


def _closed_form(P, t0, p0):
    """The expression extend_atmosphere used before carmapy.adiabat existed."""
    a, b = PARMENTIER_A_COEFF, PARMENTIER_B_COEFF
    k = t0 / (a - b * t0) * (P / p0) ** a
    return a * k / (1 + b * k)


# ---------------------------------------------------------------------------
# the analytic model
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("t0", [500.0, 1200.0, 3000.0])
def test_parmentier_matches_closed_form(t0):
    """Bit-for-bit with the pre-refactor expression."""
    got = parmentier_adiabat(P_IN_TABLE, t0, 1e7)
    assert np.array_equal(got, _closed_form(P_IN_TABLE, t0, 1e7))


def test_parmentier_gradient_is_a_minus_bt():
    """d ln T / d ln P of the fit must be a - b*T, the identity everything
    else in this module is reasoned from."""
    P = np.geomspace(1e6, 1e8, 400)
    T = parmentier_adiabat(P, 1500.0, 1e7)

    gradient = np.gradient(np.log(T), np.log(P))
    expected = PARMENTIER_A_COEFF - PARMENTIER_B_COEFF * T

    # Endpoints excluded: np.gradient is one-sided there.
    assert np.allclose(gradient[1:-1], expected[1:-1], rtol=1e-5)


# ---------------------------------------------------------------------------
# the gradient table
# ---------------------------------------------------------------------------

def test_grad_ad_reproduces_table_nodes():
    """Bilinear interpolation must be exact at the grid points."""
    t_axis, p_axis, grad = _load_grad_table()

    for it in (0, 17, 34, len(t_axis) - 1):
        for ip in (0, 9, 18, len(p_axis) - 1):
            got = grad_ad(10 ** t_axis[it], 10 ** p_axis[ip] * BAR_TO_BARYE)
            assert got == pytest.approx(grad[it, ip], rel=1e-12)


def test_grad_ad_shows_the_dissociation_dip():
    """The feature the analytic fit cannot represent: at low pressure the
    gradient collapses well below the fit near 2500 K."""
    T = 2512.0
    tabulated = grad_ad(T, 1e-2 * BAR_TO_BARYE)
    fit = PARMENTIER_A_COEFF - PARMENTIER_B_COEFF * T

    assert tabulated < 0.5 * fit


def test_grad_ad_is_pressure_independent_when_cool():
    """Below ~1500 K the mixture is ideal and the table has no pressure
    dependence, which is why the fit works there at all."""
    values = [grad_ad(800.0, p * BAR_TO_BARYE)
              for p in (1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3)]

    assert np.allclose(values, values[0], rtol=1e-12)


def test_grad_ad_clamps_and_warns_outside_the_table():
    t_axis, p_axis, grad = _load_grad_table()

    with pytest.warns(UserWarning, match="clamped"):
        hot = grad_ad(50_000.0, 1e6 * BAR_TO_BARYE)

    assert hot == pytest.approx(grad[-1, -1], rel=1e-12)


def test_grad_ad_does_not_warn_inside_the_table():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        grad_ad(1000.0, 1.0 * BAR_TO_BARYE)


# ---------------------------------------------------------------------------
# the tabulated model
# ---------------------------------------------------------------------------

def test_tabulated_returns_the_anchor_exactly():
    assert tabulated_adiabat(np.array([1e7]), 1800.0, 1e7)[0] == 1800.0


def test_tabulated_integrates_its_own_gradient():
    """The profile it produces must have grad_ad as its actual gradient."""
    P = np.geomspace(1e5, 1e8, 300)
    T = tabulated_adiabat(P, 1200.0, 1e7)

    gradient = np.gradient(np.log(T), np.log(P))
    expected = grad_ad(T, P)

    # Endpoints excluded: np.gradient is one-sided there.
    assert np.allclose(gradient[2:-2], expected[2:-2], rtol=2e-3)


def test_tabulated_agrees_with_the_fit_when_cool():
    """Below ~1500 K the two models describe the same physics.

    The gradients agree to better than 2% there, but integrating over three
    decades of pressure accumulates that into a slightly larger difference in
    temperature, hence 3% rather than 2%.
    """
    P = np.geomspace(1e5, 1e8, 20)
    fit = parmentier_adiabat(P, 800.0, 1e7)
    table = tabulated_adiabat(P, 800.0, 1e7)

    cool = fit < 1500.0
    assert cool.sum() > 5
    assert np.allclose(table[cool], fit[cool], rtol=0.03)


def test_tabulated_runs_cooler_than_the_fit_when_hot():
    """Dissociation flattens the real adiabat, so integrating through it must
    land below the fit. The fit's error is what motivates this module."""
    # 1 to 20 bar from a 2000 K anchor stays inside the table, so this is the
    # table's own statement rather than an extrapolation.
    P = np.geomspace(1e6, 2e7, 40)
    fit = parmentier_adiabat(P, 2000.0, 1e6)
    table = tabulated_adiabat(P, 2000.0, 1e6)

    assert table[-1] < 0.95 * fit[-1]


def test_both_models_are_order_independent():
    """Callers pass pressure grids in whatever order they have them in."""
    shuffled = np.array([7, 0, 19, 3, 12])

    for model in (parmentier_adiabat, tabulated_adiabat):
        full = model(P_IN_TABLE, 1200.0, 1e7)
        subset = model(P_IN_TABLE[shuffled], 1200.0, 1e7)
        assert np.allclose(subset, full[shuffled], rtol=1e-5)


@pytest.mark.parametrize("model", [parmentier_adiabat, tabulated_adiabat])
def test_both_models_increase_with_pressure(model):
    T = model(P_IN_TABLE, 1200.0, 1e7)
    assert np.all(np.diff(T) > 0)


@pytest.mark.parametrize("model", [parmentier_adiabat, tabulated_adiabat])
def test_both_models_preserve_shape(model):
    P = P_IN_TABLE.reshape(5, 5)
    assert model(P, 1200.0, 1e7).shape == (5, 5)


# ---------------------------------------------------------------------------
# selection
# ---------------------------------------------------------------------------

def test_get_adiabat_returns_the_named_model():
    assert get_adiabat(ADIABAT_PARMENTIER) is parmentier_adiabat
    assert get_adiabat(ADIABAT_TABLE) is tabulated_adiabat


def test_get_adiabat_rejects_unknown_names():
    with pytest.raises(ValueError, match="unknown adiabat"):
        get_adiabat("saumon")


def test_models_are_interchangeable():
    """Same signature, same return shape -- the drop-in property."""
    results = [model(P_IN_TABLE, 1200.0, 1e7)
               for model in ADIABAT_MODELS.values()]

    assert all(r.shape == P_IN_TABLE.shape for r in results)
