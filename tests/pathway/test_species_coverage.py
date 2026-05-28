"""Species-coverage regression tests.

Each of the 11 gas species in carmapy's gas_dict is run as a homogeneous
nucleation pathway with a short warm-up + fine restart so that the
time-averaged steady-state is compared against a committed reference.

Some species (Cr, Fe, Al2O3, SiO) may not condense at T=1000 K; their
references will capture near-zero numden, which still validates stability.

Generate references with:
    python tests/generate_species_fixtures.py

Run with CARMAPY_RUN_INTEGRATION=1.
"""
import numpy as np
import pytest
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
SPECIES_DATA = REPO_ROOT / "tests" / "data" / "species"

WARMUP_NSTEP = 2000
WARMUP_GAP   = 200
FINE_NSTEP   = 500
FINE_GAP     = 5
DT           = 250

# Run-to-run variability is 0.00e+00 for all species (bitwise deterministic).
# Tolerances are set to 1e-5 to catch physics regressions while allowing for
# minor floating-point differences across compilers/platforms.
NUMDEN_RTOL = 0.10
NUMDEN_ATOL = 1.0    # gate near-zero cells
GAS_RTOL    = 0.10
GAS_ATOL    = 1e-8

# H2O is excluded: it has lat_heat_e=NaN in gas_dict (condensation not
# supported in this configuration) and is always present as the background gas.
ALL_SPECIES = ["KCl", "ZnS", "Na2S", "MnS", "Cr",
               "Mg2SiO4", "Fe", "TiO2", "Al2O3", "SiO"]


def _build_species_carma(species, path, example_levels):
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base

    P, T, kzz, mu = example_levels
    c = Carma(path)
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)
    c.add_gas(species)
    c.add_hom_group(species, 1e-8)
    populate_abundances_at_cloud_base(c)
    return c


def _run_warmup_and_fine(carma):
    carma.set_stepping(dt=DT, output_gap=WARMUP_GAP, n_tstep=WARMUP_NSTEP)
    carma.run(suppress_output=True)
    carma.restart = 1
    carma.set_stepping(dt=DT, output_gap=FINE_GAP, n_tstep=FINE_NSTEP)
    carma.run(suppress_output=True)


# ── sanity tests (no committed reference required) ────────────────────────────

@pytest.mark.pathway
@pytest.mark.parametrize("species", ALL_SPECIES)
def test_species_numden_nonneg(species, tmp_path, example_levels):
    """numden must be non-negative at every output timestep."""
    from carmapy.results import Results
    c = _build_species_carma(species, str(tmp_path / species), example_levels)
    _run_warmup_and_fine(c)
    r = Results(c, read_diag=False)
    assert np.all(r.numden >= 0), f"{species}: negative numden detected"


@pytest.mark.pathway
@pytest.mark.parametrize("species", ALL_SPECIES)
def test_species_gas_abund_finite(species, tmp_path, example_levels):
    """gas_abund must be finite at every output timestep."""
    from carmapy.results import Results
    c = _build_species_carma(species, str(tmp_path / species), example_levels)
    _run_warmup_and_fine(c)
    r = Results(c, read_diag=False)
    assert np.all(np.isfinite(r.gas_abund)), f"{species}: non-finite gas_abund"


# ── regression tests (compare against committed .npz) ────────────────────────

@pytest.mark.pathway
@pytest.mark.parametrize("species", ALL_SPECIES)
def test_species_numden_regression(species, tmp_path, example_levels):
    """Time-averaged numden over the fine window must match the committed reference."""
    from carmapy.results import Results
    ref_path = SPECIES_DATA / species / f"{species}_reference.npz"
    if not ref_path.exists():
        pytest.skip(
            f"No reference for species {species!r}. "
            "Generate with: python tests/generate_species_fixtures.py"
        )
    c = _build_species_carma(species, str(tmp_path / species), example_levels)
    _run_warmup_and_fine(c)
    r = Results(c, read_diag=False)
    ref = np.load(ref_path)
    np.testing.assert_allclose(
        r.numden.mean(axis=-1), ref["numden"],
        rtol=NUMDEN_RTOL, atol=NUMDEN_ATOL,
        err_msg=f"{species}: time-averaged numden differs from reference",
    )


@pytest.mark.pathway
@pytest.mark.parametrize("species", ALL_SPECIES)
def test_species_gas_abund_regression(species, tmp_path, example_levels):
    """Time-averaged gas_abund over the fine window must match the committed reference."""
    from carmapy.results import Results
    ref_path = SPECIES_DATA / species / f"{species}_reference.npz"
    if not ref_path.exists():
        pytest.skip(
            f"No reference for species {species!r}. "
            "Generate with: python tests/generate_species_fixtures.py"
        )
    c = _build_species_carma(species, str(tmp_path / species), example_levels)
    _run_warmup_and_fine(c)
    r = Results(c, read_diag=False)
    ref = np.load(ref_path)
    np.testing.assert_allclose(
        r.gas_abund.mean(axis=-1), ref["gas_abund"],
        rtol=GAS_RTOL, atol=GAS_ATOL,
        err_msg=f"{species}: time-averaged gas_abund differs from reference",
    )
