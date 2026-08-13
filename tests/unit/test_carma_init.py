import pytest
from carmapy.carmapy import Carma
from carmapy.constants import I_CART, I_LOGP


@pytest.mark.unit
def test_default_has_h2o():
    c = Carma("_t")
    assert "H2O" in c.gases


@pytest.mark.unit
def test_default_igridv_1d():
    c = Carma("_t")
    assert c.igridv == I_CART
    assert c.is_2d is False


@pytest.mark.unit
def test_2d_igridv():
    c = Carma("_t", is_2d=True)
    assert c.igridv == I_LOGP
    assert c.is_2d is True
