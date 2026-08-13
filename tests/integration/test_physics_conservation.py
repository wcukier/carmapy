"""Physics invariant integration tests.

Each test asserts a universally required property of any valid CARMA run —
no NaN/Inf, non-negative number density, supersaturated initial conditions
must produce particles, and growth/evaporation must not fire simultaneously
(they are opposing processes). Comparative or configuration-dependent
"physics" tests were removed because they don't cleanly isolate a single
process.

All tests require CARMAPY_RUN_INTEGRATION=1.
"""
import pytest
import numpy as np


def _build_tio2(path, example_levels, *, n_tstep=50, output_gap=5, dt=250):
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base
    P, T, kzz, mu = example_levels
    c = Carma(path)
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.set_stepping(dt=dt, output_gap=output_gap, n_tstep=n_tstep)
    c.add_gas("TiO2")
    c.add_hom_group("TiO2", 1e-8)
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)
    populate_abundances_at_cloud_base(c)
    return c


@pytest.mark.integration
def test_no_nans_or_infs(tmp_path, example_levels):
    """Numerical sanity: no NaN or Inf in any field."""
    from carmapy.results import Results
    c = _build_tio2(str(tmp_path / "nan_check"), example_levels)
    c.run(suppress_output=True)
    r = Results(c, read_diag=False)
    assert np.all(np.isfinite(r.numden)), "NaN or Inf in numden"
    assert np.all(np.isfinite(r.gas_abund)), "NaN or Inf in gas_abund"


@pytest.mark.integration
def test_numden_nonneg_all_timesteps(tmp_path, example_levels):
    """Physical invariant: number density cannot be negative at any timestep."""
    from carmapy.results import Results
    c = _build_tio2(str(tmp_path / "nonneg"), example_levels)
    c.run(suppress_output=True)
    r = Results(c, read_diag=False)
    assert np.all(r.numden >= 0), "Negative numden detected"


@pytest.mark.integration
def test_nucleation_creates_particles(tmp_path, example_levels):
    """Given a supersaturated initial condition (populate_abundances_at_cloud_base
    fills gas to the saturation limit), at least one timestep must show particles
    forming. If total numden is identical at all times, nucleation is dead."""
    from carmapy.results import Results
    c = _build_tio2(str(tmp_path / "nuc_parts"), example_levels)
    c.run(suppress_output=True)
    r = Results(c, read_diag=False)
    total_initial = np.sum(r.numden[:, :, :, 0])
    total_final = np.sum(r.numden[:, :, :, -1])
    assert total_final > total_initial, (
        "Total numden did not increase from initial state — nucleation is not firing"
    )


@pytest.mark.integration
def test_grow_and_evap_not_simultaneously_large(tmp_path, example_levels):
    """Growth and evaporation are opposing processes at any single point:
    a level/bin/timestep is either net-condensing or net-evaporating, not
    both at large magnitudes. Both rates being large simultaneously would
    indicate a logic bug in CARMA's grow/evap dispatcher."""
    from carmapy.results import Results
    c = _build_tio2(str(tmp_path / "grow_evap"), example_levels)
    c.run(suppress_output=True)
    r = Results(c, read_diag=True)
    grow = r.clouds["Pure TiO2"]["grow_gain_rate"]  # (NZ, NBIN, NT)
    evap = r.clouds["Pure TiO2"]["evap_gain_rate"]
    threshold = 0.1 * max(grow.max(), evap.max())
    both_large = np.sum((grow > threshold) & (evap > threshold))
    assert both_large == 0, (
        f"{both_large} level/bin/time points have both large growth and evaporation rates"
    )
