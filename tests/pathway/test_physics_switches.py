"""Physics-switch regression tests.

Each variant modifies a single physics knob relative to the TiO2 hom
baseline (or Mg2SiO4 het baseline) and is compared against its own
committed time-averaged reference. This catches regressions in individual
physics paths that the tutorial comparison would obscure.

Variants:
  coag_on       -- TiO2 hom + self-coagulation
  coag_off      -- TiO2 hom, no coagulation
  het_on        -- Mg2SiO4 het nucleation on TiO2 seed
  het_off       -- Mg2SiO4 hom nucleation (no seed)
  nbin_20       -- TiO2 hom, NBIN=20
  nbin_40       -- TiO2 hom, NBIN=40
  high_lat_heat -- TiO2 hom, latent heat ×100 (suppresses condensation)
  high_viscosity-- TiO2 hom, rmu_1 ×1000 (suppresses sedimentation)
  wind_up       -- TiO2 hom, uniform +100 cm/s upward wind
  wind_down     -- TiO2 hom, uniform −100 cm/s downward wind

Generate references with:
    python tests/generate_physics_fixtures.py

Run with CARMAPY_RUN_INTEGRATION=1.
"""
import math
import numpy as np
import pytest
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
PHYSICS_DATA = REPO_ROOT / "tests" / "data" / "physics"

WARMUP_NSTEP = 2000
WARMUP_GAP   = 200
FINE_NSTEP   = 500
FINE_GAP     = 5
DT           = 250

# Tolerances: set after running generate_physics_fixtures.py --measure-error.
# Current values are intentionally generous until empirically tightened.
NUMDEN_RTOL = 0.10
NUMDEN_ATOL = 1.0    # gate near-zero cells
GAS_RTOL    = 0.10
GAS_ATOL    = 1e-8

ALL_VARIANTS = [
    "coag_on", "coag_off",
    "het_on", "het_off",
    "nbin_20", "nbin_40",
    "high_lat_heat", "high_viscosity",
    "wind_up", "wind_down",
]

# Pure H2 viscosity parameters (Gao et al. 2023)
_RMU_1 = 1.7970e-6
_RMU_2 = 0.685
_RMU_3 = -0.59
_RMU_4 = 140
_THCOND_0 = 7992.77
_THCOND_1 = 38.08
_THCOND_2 = -1.2585e-4
_CP = 1.3e8


def _build_switch_carma(variant, path, example_levels):
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base

    P, T, kzz, mu = example_levels
    c = Carma(path)
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)

    if variant in ("nbin_20", "nbin_40"):
        c.NBIN = int(variant.split("_")[1])

    c.calculate_z(mu)

    if variant in ("coag_on", "coag_off",
                   "nbin_20", "nbin_40",
                   "high_lat_heat", "high_viscosity",
                   "wind_up", "wind_down"):
        c.add_gas("TiO2")
        group = c.add_hom_group("TiO2", 1e-8)
        if variant == "high_lat_heat":
            # Condensate (latent heat, vapor pressure) properties are owned by
            # the group. Natural latent heat: L = vp_tcoeff * ln(10) * R_gas / wtmol_dif
            natural_lhe = group.vp_tcoeff * math.log(10) * 8.314e7 / group.gas.wtmol_dif
            group.lat_heat_e = natural_lhe * 100
        if variant == "coag_on":
            c.add_coag("Pure TiO2")

    elif variant == "het_on":
        c.add_gas("TiO2")
        c.add_gas("Mg2SiO4")
        c.add_hom_group("TiO2", 1e-8)
        c.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2 ** (1 / 3))

    elif variant == "het_off":
        c.add_gas("Mg2SiO4")
        c.add_hom_group("Mg2SiO4", 1e-8)

    if variant == "high_viscosity":
        c.set_atmospheric_parameters(
            rmu_1=_RMU_1 * 1000,
            rmu_2=_RMU_2,
            rmu_3=_RMU_3,
            rmu_4=_RMU_4,
            thcond_0=_THCOND_0,
            thcond_1=_THCOND_1,
            thcond_2=_THCOND_2,
            cp=_CP,
        )

    if variant == "wind_up":
        c.add_vertical_winds(np.full(c.NZ, +100.0))
    elif variant == "wind_down":
        c.add_vertical_winds(np.full(c.NZ, -100.0))

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
@pytest.mark.parametrize("variant", ALL_VARIANTS)
def test_switch_numden_nonneg(variant, tmp_path, example_levels):
    """numden must be non-negative at every output timestep."""
    from carmapy.results import Results
    c = _build_switch_carma(variant, str(tmp_path / variant), example_levels)
    _run_warmup_and_fine(c)
    r = Results(c, read_diag=False)
    assert np.all(r.numden >= 0), f"{variant}: negative numden detected"


@pytest.mark.pathway
@pytest.mark.parametrize("variant", ALL_VARIANTS)
def test_switch_gas_abund_finite(variant, tmp_path, example_levels):
    """gas_abund must be finite at every output timestep."""
    from carmapy.results import Results
    c = _build_switch_carma(variant, str(tmp_path / variant), example_levels)
    _run_warmup_and_fine(c)
    r = Results(c, read_diag=False)
    assert np.all(np.isfinite(r.gas_abund)), f"{variant}: non-finite gas_abund"


# ── regression tests (compare against committed .npz) ────────────────────────

@pytest.mark.pathway
@pytest.mark.parametrize("variant", ALL_VARIANTS)
def test_switch_numden_regression(variant, tmp_path, example_levels):
    """Time-averaged numden over the fine window must match the committed reference."""
    from carmapy.results import Results
    ref_path = PHYSICS_DATA / variant / f"{variant}_reference.npz"
    if not ref_path.exists():
        pytest.skip(
            f"No reference for variant {variant!r}. "
            "Generate with: python tests/generate_physics_fixtures.py"
        )
    c = _build_switch_carma(variant, str(tmp_path / variant), example_levels)
    _run_warmup_and_fine(c)
    r = Results(c, read_diag=False)
    ref = np.load(ref_path)
    np.testing.assert_allclose(
        r.numden.mean(axis=-1), ref["numden"],
        rtol=NUMDEN_RTOL, atol=NUMDEN_ATOL,
        err_msg=f"{variant}: time-averaged numden differs from reference",
    )


@pytest.mark.pathway
@pytest.mark.parametrize("variant", ALL_VARIANTS)
def test_switch_gas_abund_regression(variant, tmp_path, example_levels):
    """Time-averaged gas_abund over the fine window must match the committed reference."""
    from carmapy.results import Results
    ref_path = PHYSICS_DATA / variant / f"{variant}_reference.npz"
    if not ref_path.exists():
        pytest.skip(
            f"No reference for variant {variant!r}. "
            "Generate with: python tests/generate_physics_fixtures.py"
        )
    c = _build_switch_carma(variant, str(tmp_path / variant), example_levels)
    _run_warmup_and_fine(c)
    r = Results(c, read_diag=False)
    ref = np.load(ref_path)
    np.testing.assert_allclose(
        r.gas_abund.mean(axis=-1), ref["gas_abund"],
        rtol=GAS_RTOL, atol=GAS_ATOL,
        err_msg=f"{variant}: time-averaged gas_abund differs from reference",
    )
