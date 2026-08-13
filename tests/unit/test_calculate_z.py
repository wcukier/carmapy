import pytest
import numpy as np
from carmapy.carmapy import Carma
from carmapy.constants import k_B, PROTON_MASS


@pytest.fixture
def isothermal_carma():
    """10-layer isothermal atmosphere at T=1000K, g=980, mu=2."""
    c = Carma("_t")
    c.set_physical_params(surface_grav=980.0, wt_mol=2.0)
    P = np.logspace(9, 6, 11)   # bottom (index 0) = highest P
    T = np.ones(11) * 1000.0
    kzz = np.ones(11) * 1e8
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    return c


@pytest.mark.unit
def test_calculate_z_z_levels_monotonic(isothermal_carma):
    c = isothermal_carma
    c.calculate_z()
    assert np.all(np.diff(c.z_levels) > 0), "z_levels must increase upward"


@pytest.mark.unit
def test_calculate_z_bottom_is_zero(isothermal_carma):
    isothermal_carma.calculate_z()
    assert isothermal_carma.z_levels[0] == pytest.approx(0.0)


@pytest.mark.unit
def test_calculate_z_centers_between_levels(isothermal_carma):
    c = isothermal_carma
    c.calculate_z()
    expected = (c.z_levels[:-1] + c.z_levels[1:]) / 2
    np.testing.assert_allclose(c.z_centers, expected, rtol=1e-10)


@pytest.mark.unit
def test_calculate_z_isothermal_scale_height(isothermal_carma):
    """Each dz must equal H * ln(P_bot / P_top) for an isothermal atmosphere."""
    c = isothermal_carma
    c.calculate_z()
    H = k_B * 1000.0 / (2.0 * PROTON_MASS * 980.0)
    expected_dz = H * np.log(c.P_levels[:-1] / c.P_levels[1:])
    actual_dz = np.diff(c.z_levels)
    np.testing.assert_allclose(actual_dz, expected_dz, rtol=1e-6)


@pytest.mark.unit
def test_calculate_z_uses_stored_wt_mol(isothermal_carma):
    """calculate_z() with no argument should use self.wt_mol."""
    c = isothermal_carma
    c.calculate_z()       # uses self.wt_mol = 2.0
    z1 = c.z_levels.copy()
    c.calculate_z(wt_mol=2.0)   # explicit same value
    np.testing.assert_allclose(c.z_levels, z1, rtol=1e-12)


@pytest.mark.unit
def test_calculate_z_raises_without_wt_mol():
    c = Carma("_t")
    c.set_physical_params(surface_grav=980.0)
    P = np.logspace(9, 6, 11)
    T = np.ones(11) * 1000.0
    kzz = np.ones(11) * 1e8
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    with pytest.raises(RuntimeError):
        c.calculate_z()


@pytest.mark.unit
def test_calculate_z_raises_without_T():
    c = Carma("_t")
    c.set_physical_params(surface_grav=980.0, wt_mol=2.0)
    c.add_P(np.logspace(9, 6, 11))
    with pytest.raises(RuntimeError):
        c.calculate_z()


@pytest.mark.unit
def test_calculate_z_raises_without_surface_grav():
    c = Carma("_t")
    c.set_physical_params(wt_mol=2.0)
    P = np.logspace(9, 6, 11)
    T = np.ones(11) * 1000.0
    kzz = np.ones(11) * 1e8
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    with pytest.raises(RuntimeError):
        c.calculate_z()


@pytest.mark.unit
def test_extend_atmosphere_increases_nz(isothermal_carma):
    c = isothermal_carma
    c.calculate_z()
    original_nz = c.NZ
    c.extend_atmosphere(max_P=1e12)
    assert c.NZ > original_nz


@pytest.mark.unit
def test_extend_atmosphere_noop_when_shallow(isothermal_carma):
    c = isothermal_carma
    c.calculate_z()
    original_nz = c.NZ
    c.extend_atmosphere(max_P=1e3)  # shallower than current bottom
    assert c.NZ == original_nz


@pytest.mark.unit
def test_extend_atmosphere_pressure_monotonic(isothermal_carma):
    c = isothermal_carma
    c.calculate_z()
    c.extend_atmosphere(max_P=1e12)
    assert np.all(np.diff(c.P_levels) < 0), "P_levels must decrease bottom→top"
