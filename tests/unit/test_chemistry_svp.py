import pytest
import numpy as np
from carmapy.chemistry import saturation_vapor_pressure
from carmapy.constants import gas_dict, TCOEFF_KCL


@pytest.mark.unit
def test_svp_kcl_analytic():
    """KCl SVP formula: 1e6 * 10^(7.6106 - 11382/T), no met or P dependence."""
    T = np.array([1000.0])
    P = np.array([1e6])
    expected = 1e6 * 10 ** (7.6106 - TCOEFF_KCL / 1000.0)
    result = saturation_vapor_pressure(P, T, 0.0, "KCl")
    np.testing.assert_allclose(result, expected, rtol=1e-10)


@pytest.mark.unit
@pytest.mark.parametrize("species", ["KCl", "ZnS", "TiO2", "Mg2SiO4", "Fe", "Cr"])
def test_svp_positive(species):
    """SVP must be positive at physically reasonable temperatures."""
    T = np.array([1000.0, 1500.0, 2000.0])
    P = np.ones(3) * 1e6
    result = saturation_vapor_pressure(P, T, 0.0, species)
    assert np.all(result > 0), f"{species} SVP not positive"


@pytest.mark.unit
@pytest.mark.parametrize("species", ["KCl", "ZnS", "TiO2", "Mg2SiO4", "Fe", "Cr"])
def test_svp_increases_with_temperature(species):
    """SVP must be monotonically increasing with temperature."""
    T = np.array([500.0, 1000.0, 1500.0, 2000.0])
    P = np.ones(4) * 1e6
    result = saturation_vapor_pressure(P, T, 0.0, species)
    assert np.all(np.diff(result) > 0), f"{species} SVP not monotone increasing"


@pytest.mark.unit
def test_svp_zns_decreases_with_metallicity():
    """ZnS has vp_metcoeff > 0, so higher log_met → lower SVP."""
    T = np.array([1200.0])
    P = np.array([1e6])
    svp_solar = saturation_vapor_pressure(P, T, 0.0, "ZnS")
    svp_metal = saturation_vapor_pressure(P, T, 1.0, "ZnS")
    assert svp_metal < svp_solar


@pytest.mark.unit
def test_svp_kcl_no_metallicity_dependence():
    """KCl has vp_metcoeff=0, so SVP should not change with metallicity."""
    T = np.array([1000.0])
    P = np.array([1e6])
    svp_solar = saturation_vapor_pressure(P, T, 0.0, "KCl")
    svp_metal = saturation_vapor_pressure(P, T, 2.0, "KCl")
    np.testing.assert_allclose(svp_solar, svp_metal, rtol=1e-10)


@pytest.mark.unit
def test_svp_gas_object_matches_string():
    """Passing a Gas object should give the same result as the species name string."""
    from carmapy.carmapy import Gas
    T = np.array([1200.0])
    P = np.array([1e6])
    gas_obj = Gas("KCl", 1)
    svp_str = saturation_vapor_pressure(P, T, 0.0, "KCl")
    svp_obj = saturation_vapor_pressure(P, T, 0.0, gas_obj)
    np.testing.assert_allclose(svp_str, svp_obj, rtol=1e-10)
