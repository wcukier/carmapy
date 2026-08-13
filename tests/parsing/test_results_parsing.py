"""Tests for Results() parsing using the committed mini_output fixture.

All tests in this file use the `mini_carma` and `mini_results` session fixtures
from conftest.py, which skip automatically if the fixture files don't exist.
Run `python tests/generate_fixture.py` once to generate them.
"""
import pytest
import numpy as np


@pytest.mark.parsing
def test_numden_ndim(mini_results):
    assert mini_results.numden.ndim == 4


@pytest.mark.parsing
def test_numden_shape_nz(mini_results, mini_carma):
    assert mini_results.numden.shape[0] == mini_carma.NZ


@pytest.mark.parsing
def test_numden_shape_ngroup(mini_results, mini_carma):
    assert mini_results.numden.shape[1] == len(mini_carma.groups)


@pytest.mark.parsing
def test_numden_shape_nbin(mini_results, mini_carma):
    assert mini_results.numden.shape[2] == mini_carma.NBIN


@pytest.mark.parsing
def test_numden_nonnegative(mini_results):
    assert np.all(mini_results.numden >= 0)


@pytest.mark.parsing
def test_gas_abund_shape(mini_results, mini_carma):
    assert mini_results.gas_abund.shape[0] == mini_carma.NZ
    assert mini_results.gas_abund.shape[1] == len(mini_carma.gases)


@pytest.mark.parsing
def test_sat_vp_shape(mini_results, mini_carma):
    assert mini_results.sat_vp.shape[0] == mini_carma.NZ
    assert mini_results.sat_vp.shape[1] == len(mini_carma.gases)


@pytest.mark.parsing
def test_sat_vp_nonnegative(mini_results):
    assert np.all(mini_results.sat_vp >= 0)


@pytest.mark.parsing
def test_pressure_shape(mini_results, mini_carma):
    assert mini_results.P.shape == (mini_carma.NZ,)


@pytest.mark.parsing
def test_temperature_shape(mini_results, mini_carma):
    assert mini_results.T.shape == (mini_carma.NZ,)


@pytest.mark.parsing
def test_altitude_shape(mini_results, mini_carma):
    assert mini_results.Z.shape == (mini_carma.NZ,)


@pytest.mark.parsing
def test_pressure_monotonic(mini_results):
    """Pressure must be monotone (either all increasing or all decreasing)."""
    P = mini_results.P
    assert np.all(np.diff(P) < 0) or np.all(np.diff(P) > 0)


@pytest.mark.parsing
def test_ts_monotonically_increasing(mini_results):
    assert np.all(np.diff(mini_results.ts) > 0)


@pytest.mark.parsing
def test_dt_timestep(mini_results, mini_carma):
    assert mini_results.dt_timestep == mini_carma.dt * mini_carma.output_gap


@pytest.mark.parsing
def test_group_names_match(mini_results, mini_carma):
    assert mini_results.group_names == list(mini_carma.groups.keys())


@pytest.mark.parsing
def test_gas_names_match(mini_results, mini_carma):
    assert mini_results.gas_names == list(mini_carma.gases.keys())


@pytest.mark.parsing
def test_clouds_has_all_groups(mini_results, mini_carma):
    for name in mini_carma.groups:
        assert name in mini_results.clouds


@pytest.mark.parsing
def test_clouds_numden_shape(mini_results, mini_carma):
    for name in mini_carma.groups:
        nd = mini_results.clouds[name]["numden"]
        assert nd.ndim == 3
        assert nd.shape[0] == mini_carma.NZ
        assert nd.shape[1] == mini_carma.NBIN


@pytest.mark.parsing
def test_clouds_r_shape(mini_results, mini_carma):
    for name in mini_carma.groups:
        r = mini_results.clouds[name]["r"]
        assert r.shape == (mini_carma.NBIN,)


@pytest.mark.parsing
def test_clouds_r_positive(mini_results, mini_carma):
    for name in mini_carma.groups:
        assert np.all(mini_results.clouds[name]["r"] > 0)


@pytest.mark.parsing
def test_clouds_r_monotonically_increasing(mini_results, mini_carma):
    """Particle size bins must be ordered smallest to largest."""
    for name in mini_carma.groups:
        r = mini_results.clouds[name]["r"]
        assert np.all(np.diff(r) > 0), f"Radii not monotone for {name}"


@pytest.mark.parsing
def test_gases_dict_excludes_h2o(mini_results):
    """H2O (index 0) is excluded from results.gases by convention."""
    assert "H2O" not in mini_results.gases


@pytest.mark.parsing
def test_gases_dict_has_tio2(mini_results):
    assert "TiO2" in mini_results.gases
