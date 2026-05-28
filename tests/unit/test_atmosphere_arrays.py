import pytest
import numpy as np
from carmapy.carmapy import Carma


@pytest.mark.unit
def test_add_P_sets_nz():
    c = Carma("_t")
    P = np.linspace(1e9, 1e6, 11)  # 11 levels → NZ=10
    c.add_P(P)
    assert c.NZ == 10


@pytest.mark.unit
def test_add_P_computes_centers():
    c = Carma("_t")
    P = np.array([300.0, 200.0, 100.0])
    c.add_P(P)
    expected = np.array([250.0, 150.0])
    np.testing.assert_allclose(c.P_centers, expected)


@pytest.mark.unit
def test_add_P_wrong_shape_raises():
    c = Carma("_t")
    P = np.linspace(1e9, 1e6, 11)
    c.add_P(P)
    with pytest.raises(ValueError):
        c.add_P(np.linspace(1e9, 1e6, 7))  # wrong length


@pytest.mark.unit
def test_add_T_sets_nz():
    c = Carma("_t")
    T = np.linspace(2000, 1000, 11)
    c.add_T(T)
    assert c.NZ == 10


@pytest.mark.unit
def test_add_T_computes_centers():
    c = Carma("_t")
    T = np.array([2000.0, 1500.0, 1000.0])
    c.add_T(T)
    expected = np.array([1750.0, 1250.0])
    np.testing.assert_allclose(c.T_centers, expected)


@pytest.mark.unit
def test_add_T_wrong_shape_raises():
    c = Carma("_t")
    P = np.linspace(1e9, 1e6, 11)
    c.add_P(P)  # NZ=10, so T must have 11 elements
    with pytest.raises(ValueError):
        c.add_T(np.linspace(2000, 1000, 8))


@pytest.mark.unit
def test_add_T_2d_array_raises_in_1d_mode():
    c = Carma("_t")
    T_2d = np.ones((11, 4))
    with pytest.raises(ValueError):
        c.add_T(T_2d)


@pytest.mark.unit
def test_add_T_1d_array_raises_in_2d_mode():
    c = Carma("_t", is_2d=True)
    with pytest.raises(ValueError):
        c.add_T(np.linspace(2000, 1000, 11))


@pytest.mark.unit
def test_add_kzz_sets_nz():
    c = Carma("_t")
    kzz = np.ones(11) * 1e8
    c.add_kzz(kzz)
    assert c.NZ == 10


@pytest.mark.unit
def test_add_kzz_wrong_shape_raises():
    c = Carma("_t")
    P = np.linspace(1e9, 1e6, 11)
    c.add_P(P)
    with pytest.raises(ValueError):
        c.add_kzz(np.ones(5) * 1e8)


@pytest.mark.unit
def test_add_z_sets_nz():
    c = Carma("_t")
    z = np.linspace(0, 1e8, 11)
    c.add_z(z)
    assert c.NZ == 10


@pytest.mark.unit
def test_add_z_computes_centers():
    c = Carma("_t")
    z = np.array([0.0, 1e6, 2e6])
    c.add_z(z)
    expected = np.array([0.5e6, 1.5e6])
    np.testing.assert_allclose(c.z_centers, expected)


@pytest.mark.unit
def test_add_z_wrong_shape_raises():
    c = Carma("_t")
    P = np.linspace(1e9, 1e6, 11)
    c.add_P(P)
    with pytest.raises(ValueError):
        c.add_z(np.linspace(0, 1e8, 5))


@pytest.mark.unit
def test_arrays_must_be_compatible():
    """Confirm NZ is set from first array and validated against subsequent ones."""
    c = Carma("_t")
    c.add_P(np.linspace(1e9, 1e6, 11))   # NZ = 10
    c.add_T(np.linspace(2000, 1000, 11))  # must be 11 elements — OK
    assert c.NZ == 10
