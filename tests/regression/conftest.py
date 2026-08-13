"""Shared fixture for regression tests.

The regression tests (tutorial-config numden/gas_abund check and brightness-T
spectrum check) both require the same expensive simulation: a 24000-step
warm-up at output_gap=800 followed by a 10000-step fine restart at
output_gap=10. This module provides a single session-scoped fixture so the
warm-up+restart runs once per pytest invocation regardless of how many
regression tests consume it.

Configuration mirrors docs/source/notebooks/1_my_first_carma.ipynb exactly
(dt=100, wt_mol=mu[0], Mg2SiO4 het nucleation on TiO2 seed,
populate_abundances_at_cloud_base, extended atmosphere).
"""
import pytest


WARMUP_NSTEP = 24000
WARMUP_GAP = 800
FINE_NSTEP = 10000
FINE_GAP = 10
DT = 100


def _build_canonical(name: str):
    from carmapy.carmapy import Carma
    from carmapy.example import example_levels
    from carmapy.chemistry import populate_abundances_at_cloud_base

    P, T, kzz, mu = example_levels()
    c = Carma(name)
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    # The condensate (group) methods create their reservoir gases internally.
    c.add_hom_group("TiO2", 1e-8)
    c.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2 ** (1 / 3))
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.set_physical_params(surface_grav=31600, wt_mol=mu[0])
    c.calculate_z()
    c.calculate_z(mu)
    c.extend_atmosphere(1e10)
    populate_abundances_at_cloud_base(c)
    return c


@pytest.fixture(scope="session")
def regression_run(tmp_path_factory):
    """Run warm-up + fine restart once per session; reused by all regression tests.

    Returns the Carma object with `results` populated from the fine-restart
    output (~1000 snapshots at gap=10 covering the post-saturation steady
    state). The simulation directory path is `carma.name`.
    """
    from carmapy.results import Results

    sim_path = tmp_path_factory.mktemp("regression") / "sim"
    carma = _build_canonical(str(sim_path))
    carma.set_stepping(dt=DT, output_gap=WARMUP_GAP, n_tstep=WARMUP_NSTEP)
    carma.run(suppress_output=True)
    carma.restart = 1
    carma.set_stepping(dt=DT, output_gap=FINE_GAP, n_tstep=FINE_NSTEP)
    carma.run(suppress_output=True)
    carma.results = Results(carma, read_diag=False)
    return carma
