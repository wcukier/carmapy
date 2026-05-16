import pytest
from carmapy.carmapy import Carma
from carmapy.constants import I_CART, I_LOGP


@pytest.mark.unit
def test_default_nbin():
    c = Carma("_t")
    assert c.NBIN == 80


@pytest.mark.unit
def test_default_has_h2o():
    c = Carma("_t")
    assert "H2O" in c.gases


@pytest.mark.unit
def test_default_igridv_1d():
    c = Carma("_t")
    assert c.igridv == I_CART


@pytest.mark.unit
def test_2d_igridv():
    c = Carma("_t", is_2d=True)
    assert c.igridv == I_LOGP


@pytest.mark.unit
def test_default_groups_empty():
    c = Carma("_t")
    assert len(c.groups) == 0


@pytest.mark.unit
def test_default_elems_empty():
    c = Carma("_t")
    assert len(c.elems) == 0


@pytest.mark.unit
def test_default_nucs_empty():
    c = Carma("_t")
    assert len(c.nucs) == 0


@pytest.mark.unit
def test_default_growth_empty():
    c = Carma("_t")
    assert len(c.growth) == 0


@pytest.mark.unit
def test_default_coags_empty():
    c = Carma("_t")
    assert len(c.coags) == 0


@pytest.mark.unit
def test_name_stored():
    c = Carma("my_simulation")
    assert c.name == "my_simulation"


@pytest.mark.unit
def test_default_nz_zero():
    c = Carma("_t")
    assert c.NZ == 0


@pytest.mark.unit
def test_default_not_2d():
    c = Carma("_t")
    assert c.is_2d is False


@pytest.mark.unit
def test_2d_flag_set():
    c = Carma("_t", is_2d=True)
    assert c.is_2d is True


@pytest.mark.unit
def test_default_restart_false():
    c = Carma("_t")
    assert c.restart is False


@pytest.mark.unit
def test_default_output_gap():
    c = Carma("_t")
    assert c.output_gap == 1000


@pytest.mark.unit
def test_default_n_tstep():
    c = Carma("_t")
    assert c.n_tstep == 1_000_000
