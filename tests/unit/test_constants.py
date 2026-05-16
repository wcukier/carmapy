import pytest
import numpy as np
from carmapy.constants import (
    gas_dict, JUPITER_RADIUS, k_B, PROTON_MASS, BAR_TO_BARYE,
    TCOEFF_TIO2_HELLING,
)

KNOWN_SPECIES = [
    "H2O", "KCl", "ZnS", "Na2S", "MnS", "Cr",
    "Mg2SiO4", "Fe", "TiO2", "Al2O3", "SiO",
]

REQUIRED_FIELDS = [
    "rho_cond", "wtmol", "wtmol_dif", "coldia",
    "vp_offset", "vp_tcoeff", "vp_metcoeff", "vp_logpcoeff",
    "surften_0", "surften_slope", "is_typeIII", "stofact",
    "hill_formula", "mucos_dict",
]


@pytest.mark.unit
@pytest.mark.parametrize("species", KNOWN_SPECIES)
def test_species_in_gas_dict(species):
    assert species in gas_dict, f"{species} missing from gas_dict"


@pytest.mark.unit
@pytest.mark.parametrize("species", KNOWN_SPECIES)
def test_required_fields_present(species):
    for field in REQUIRED_FIELDS:
        assert field in gas_dict[species], f"{species} missing field '{field}'"


@pytest.mark.unit
@pytest.mark.parametrize("species", KNOWN_SPECIES)
def test_rho_cond_positive(species):
    assert gas_dict[species]["rho_cond"] > 0


@pytest.mark.unit
@pytest.mark.parametrize("species", KNOWN_SPECIES)
def test_wtmol_positive(species):
    assert gas_dict[species]["wtmol"] > 0
    assert gas_dict[species]["wtmol_dif"] > 0


@pytest.mark.unit
@pytest.mark.parametrize("species", KNOWN_SPECIES)
def test_coldia_reasonable(species):
    """Collision diameters should be ~1–10 Ångström (1e-8 to 1e-7 cm)."""
    cd = gas_dict[species]["coldia"]
    assert 5e-9 < cd < 5e-7, f"{species} coldia={cd} looks unreasonable"


@pytest.mark.unit
def test_tio2_vp_tcoeff_helling():
    """TiO2 T-coefficient from Helling et al. 2001 — guard against accidental edits."""
    assert gas_dict["TiO2"]["vp_tcoeff"] == pytest.approx(TCOEFF_TIO2_HELLING, rel=1e-6)


@pytest.mark.unit
def test_zns_has_metcoeff():
    """ZnS has a non-zero metallicity coefficient (positive)."""
    assert gas_dict["ZnS"]["vp_metcoeff"] > 0


@pytest.mark.unit
def test_mucos_dict_is_dict():
    for species in KNOWN_SPECIES:
        assert isinstance(gas_dict[species]["mucos_dict"], dict)


@pytest.mark.unit
def test_jupiter_radius_plausible():
    assert 6e9 < JUPITER_RADIUS < 8e9  # cm


@pytest.mark.unit
def test_kb_and_proton_mass_scale_height():
    """Scale height at T=1000K, g=980, mu=2 should be physically plausible (~1e7 cm)."""
    H = k_B * 1000.0 / (2.0 * PROTON_MASS * 980.0)
    assert 1e5 < H < 1e9


@pytest.mark.unit
def test_bar_to_barye():
    """1 bar = 1e6 barye (CGS)."""
    assert BAR_TO_BARYE == pytest.approx(1e6)
