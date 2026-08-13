import pytest
from carmapy.constants import gas_dict, group_dict, TCOEFF_TIO2_HELLING

KNOWN_SPECIES = [
    "H2O", "KCl", "ZnS", "Na2S", "MnS", "Cr",
    "Mg2SiO4", "Fe", "TiO2", "Al2O3", "SiO",
]

# Condensate properties now live on the group (group_dict).
GROUP_FIELDS = [
    "rho_cond", "wtmol", "coldia",
    "vp_offset", "vp_tcoeff", "vp_metcoeff", "vp_logpcoeff",
    "surften_0", "surften_slope", "is_typeIII", "stofact",
    "mucos_dict", "gas_phase",
]

# Gas-phase properties stay on the gas (gas_dict), keyed by the gas-phase name.
GAS_FIELDS = ["wtmol_dif", "hill_formula"]


def _gas_of(species):
    """Gas-phase reservoir name for a condensate group."""
    return group_dict[species]["gas_phase"]


@pytest.mark.unit
@pytest.mark.parametrize("species", KNOWN_SPECIES)
def test_species_in_group_dict(species):
    assert species in group_dict, f"{species} missing from group_dict"


@pytest.mark.unit
def test_required_fields_present():
    for species in KNOWN_SPECIES:
        for field in GROUP_FIELDS:
            assert field in group_dict[species], \
                f"{species} missing group field '{field}'"
        gas = _gas_of(species)
        for field in GAS_FIELDS:
            assert field in gas_dict[gas], \
                f"gas '{gas}' (for {species}) missing field '{field}'"


@pytest.mark.unit
@pytest.mark.parametrize("species", KNOWN_SPECIES)
def test_rho_cond_positive(species):
    assert group_dict[species]["rho_cond"] > 0


@pytest.mark.unit
@pytest.mark.parametrize("species", KNOWN_SPECIES)
def test_wtmol_positive(species):
    assert group_dict[species]["wtmol"] > 0
    assert gas_dict[_gas_of(species)]["wtmol_dif"] > 0


@pytest.mark.unit
@pytest.mark.parametrize("species", KNOWN_SPECIES)
def test_coldia_reasonable(species):
    """Collision diameters should be ~1–10 Ångström (1e-8 to 1e-7 cm)."""
    cd = group_dict[species]["coldia"]
    assert 5e-9 < cd < 5e-7, f"{species} coldia={cd} looks unreasonable"


@pytest.mark.unit
def test_tio2_vp_tcoeff_helling():
    """TiO2 T-coefficient from Helling et al. 2001 — guard against accidental edits."""
    assert group_dict["TiO2"]["vp_tcoeff"] == pytest.approx(TCOEFF_TIO2_HELLING, rel=1e-6)


@pytest.mark.unit
def test_zns_has_metcoeff():
    """ZnS has a non-zero metallicity coefficient (positive)."""
    assert group_dict["ZnS"]["vp_metcoeff"] > 0


@pytest.mark.unit
def test_mucos_dict_is_dict():
    for species in KNOWN_SPECIES:
        assert isinstance(group_dict[species]["mucos_dict"], dict)
