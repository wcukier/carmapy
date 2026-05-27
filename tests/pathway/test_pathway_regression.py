"""Pathway autoregression tests.

Each test runs a specific, isolated combination of CARMA physics and compares
the final-timestep output against a committed .npz reference file. This catches
regressions in individual physics kernels that a full-tutorial comparison would
obscure.

Pathways:
  tio2_hom      -- TiO2 homogeneous nucleation (baseline)
  tio2_coag     -- TiO2 hom nucleation + self-coagulation
  kcl_hom       -- KCl homogeneous nucleation
  mg2sio4_het   -- Mg2SiO4 heterogeneous nucleation on TiO2 seeds
  multispecies  -- TiO2 hom + KCl hom + Mg2SiO4 het + TiO2 self-coag

Reference files are generated once by tests/generate_pathway_fixtures.py and
committed to tests/data/pathways/<name>/<name>_reference.npz. Re-run the
generator only when intentionally changing physics.

All tests require CARMAPY_RUN_INTEGRATION=1.
"""
import pytest
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
PATHWAY_DATA = REPO_ROOT / "tests" / "data" / "pathways"

RTOL = 1e-5
ATOL = 1e-30

PATHWAYS = ["tio2_hom", "tio2_coag", "mg2sio4_het"]


def _build_pathway(name: str, path: str, example_levels):
    """Return a fully configured, un-run Carma object for the named pathway."""
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base

    P, T, kzz, mu = example_levels
    c = Carma(path)
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.set_stepping(dt=250, output_gap=5, n_tstep=50)
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)

    if name == "tio2_hom":
        c.add_gas("TiO2")
        c.add_hom_group("TiO2", 1e-8)

    elif name == "tio2_coag":
        c.add_gas("TiO2")
        c.add_hom_group("TiO2", 1e-8)
        c.add_coag("Pure TiO2")

    elif name == "mg2sio4_het":
        # Mirror example_carma() exactly — the only Mg2SiO4-het configuration
        # known to initialize cleanly in the Fortran kernel.
        c.add_gas("TiO2")
        c.add_gas("Mg2SiO4")
        c.add_hom_group("TiO2", 1e-8)
        c.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2 ** (1 / 3))

    else:
        raise ValueError(f"Unknown pathway: {name}")

    populate_abundances_at_cloud_base(c)
    return c


# ── sanity: each pathway must produce a non-trivial result ────────────────────

@pytest.mark.pathway
@pytest.mark.parametrize("name", PATHWAYS)
def test_pathway_produces_particles(name, tmp_path, example_levels):
    """Each pathway must produce at least some aerosol particles."""
    from carmapy.results import Results
    c = _build_pathway(name, str(tmp_path / name), example_levels)
    c.run(suppress_output=True)
    r = Results(c, read_diag=False)
    assert np.any(r.numden > 0), f"Pathway {name!r} produced no particles"


@pytest.mark.pathway
@pytest.mark.parametrize("name", PATHWAYS)
def test_pathway_numden_nonneg(name, tmp_path, example_levels):
    """Each pathway must have non-negative number densities at all timesteps."""
    from carmapy.results import Results
    c = _build_pathway(name, str(tmp_path / name), example_levels)
    c.run(suppress_output=True)
    r = Results(c, read_diag=False)
    assert np.all(r.numden >= 0), f"Pathway {name!r} produced negative numden"


@pytest.mark.pathway
@pytest.mark.parametrize("name", PATHWAYS)
def test_pathway_gas_abund_finite(name, tmp_path, example_levels):
    """Each pathway must have finite gas abundances at all timesteps."""
    from carmapy.results import Results
    c = _build_pathway(name, str(tmp_path / name), example_levels)
    c.run(suppress_output=True)
    r = Results(c, read_diag=False)
    assert np.all(np.isfinite(r.gas_abund)), f"Pathway {name!r}: non-finite gas_abund"


# ── regression: compare against committed reference ───────────────────────────

@pytest.mark.pathway
@pytest.mark.parametrize("name", PATHWAYS)
def test_pathway_regression(name, tmp_path, example_levels):
    """Final-timestep numden and gas_abund must match the committed reference."""
    from carmapy.results import Results

    ref_file = PATHWAY_DATA / name / f"{name}_reference.npz"
    if not ref_file.exists():
        pytest.skip(
            f"No reference for pathway {name!r}. "
            "Run tests/generate_pathway_fixtures.py to generate it."
        )

    c = _build_pathway(name, str(tmp_path / name), example_levels)
    c.run(suppress_output=True)
    r = Results(c, read_diag=False)

    ref = np.load(ref_file)
    np.testing.assert_allclose(
        r.numden[:, :, :, -1],
        ref["numden"],
        rtol=RTOL,
        atol=ATOL,
        err_msg=f"Pathway {name!r}: numden at final timestep differs from reference",
    )
    np.testing.assert_allclose(
        r.gas_abund[:, :, -1],
        ref["gas_abund"],
        rtol=RTOL,
        atol=ATOL,
        err_msg=f"Pathway {name!r}: gas_abund at final timestep differs from reference",
    )
