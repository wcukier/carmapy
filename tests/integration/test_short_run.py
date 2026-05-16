"""Short end-to-end integration tests.

These tests actually invoke the Fortran binary via Carma.run(). They are
skipped by default; set CARMAPY_RUN_INTEGRATION=1 to enable them.
"""
import os
import pytest
import numpy as np


@pytest.fixture
def short_carma(tmp_path, example_levels):
    """A minimal TiO2 simulation configured for 50 timesteps."""
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base

    P, T, kzz, mu = example_levels
    name = str(tmp_path / "short_run")
    c = Carma(name)
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.set_stepping(dt=250, output_gap=5, n_tstep=50)
    c.add_gas("TiO2")
    c.add_hom_group("TiO2", 1e-8)
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)
    populate_abundances_at_cloud_base(c)
    return c


@pytest.mark.integration
def test_short_run_completes(short_carma):
    short_carma.run(suppress_output=True)


@pytest.mark.integration
def test_short_run_creates_output_file(short_carma):
    short_carma.run(suppress_output=True)
    name = os.path.basename(short_carma.name)
    out = os.path.join(short_carma.name, f"{name}.txt")
    assert os.path.exists(out)
    assert os.path.getsize(out) > 100


@pytest.mark.integration
def test_short_run_creates_flux_file(short_carma):
    short_carma.run(suppress_output=True)
    name = os.path.basename(short_carma.name)
    flux = os.path.join(short_carma.name, f"flux_{name}.txt")
    assert os.path.exists(flux)


@pytest.mark.integration
def test_short_run_output_parseable(short_carma):
    from carmapy.results import Results
    short_carma.run(suppress_output=True)
    results = Results(short_carma, read_diag=False)
    assert results.numden.shape[0] == short_carma.NZ
    assert results.numden.shape[1] == len(short_carma.groups)
    assert results.numden.shape[2] == short_carma.NBIN


@pytest.mark.integration
def test_short_run_numden_nonnegative(short_carma):
    from carmapy.results import Results
    short_carma.run(suppress_output=True)
    results = Results(short_carma, read_diag=False)
    assert np.all(results.numden >= 0)


@pytest.mark.integration
def test_short_run_gas_abund_finite(short_carma):
    from carmapy.results import Results
    short_carma.run(suppress_output=True)
    results = Results(short_carma, read_diag=False)
    assert np.all(np.isfinite(results.gas_abund))


@pytest.mark.integration
def test_short_run_ts_monotone(short_carma):
    from carmapy.results import Results
    short_carma.run(suppress_output=True)
    results = Results(short_carma, read_diag=False)
    assert np.all(np.diff(results.ts) > 0)
