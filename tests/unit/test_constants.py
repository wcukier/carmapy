import pytest
from carmapy.constants import gas_dict, group_dict, TCOEFF_TIO2_HELLING

KNOWN_SPECIES = [
    "H2O", "KCl", "ZnS", "Na2S", "MnS", "Cr",
    "Mg2SiO4", "Fe", "TiO2", "Al2O3", "SiO", "MgSiO3",
]

# Condensate properties now live on the group (group_dict).
GROUP_FIELDS = [
    "rho_cond", "wtmol",
    "vp_offset", "vp_tcoeff", "vp_metcoeff", "vp_logpcoeff",
    "surften_0", "surften_slope", "is_typeIII", "stofact",
    "mucos_dict", "gas_phase", "stoich",
]

# Gas-phase properties stay on the gas (gas_dict), keyed by the gas-phase name.
GAS_FIELDS = ["wtmol_dif", "coldia", "hill_formula"]


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
@pytest.mark.parametrize("gas", sorted(gas_dict))
def test_coldia_reasonable(gas):
    """Collision diameters should be ~1–10 Ångström (1e-8 to 1e-7 cm)."""
    cd = gas_dict[gas]["coldia"]
    assert cd is not None, f"gas '{gas}' has no collision diameter (TODO in constants.py)"
    assert 5e-9 < cd < 5e-7, f"gas '{gas}' coldia={cd} looks unreasonable"


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


@pytest.mark.unit
@pytest.mark.parametrize("species", KNOWN_SPECIES)
def test_stoich_schema(species):
    """stoich is the full depletion vector; its key-gas entry equals stofact
    and every gas it names is a defined gas reservoir."""
    entry = group_dict[species]
    stoich = entry["stoich"]
    assert isinstance(stoich, dict) and len(stoich) >= 1
    assert entry["gas_phase"] in stoich, \
        f"{species}: key gas '{entry['gas_phase']}' missing from stoich"
    assert stoich[entry["gas_phase"]] == entry["stofact"], \
        f"{species}: stoich[key gas] != stofact"
    for gas, nu in stoich.items():
        assert gas in gas_dict, f"{species}: stoich names unknown gas '{gas}'"
        assert isinstance(nu, int) and nu > 0


@pytest.mark.unit
def test_gcomp_unique():
    """Each gas reservoir has a distinct composition flag."""
    gcomps = [v["gcomp"] for v in gas_dict.values()]
    assert len(gcomps) == len(set(gcomps)), "duplicate gcomp values in gas_dict"


@pytest.mark.unit
def test_silicates_share_sio_reservoir():
    """MgSiO3 draws from the SiO reservoir; Mg2SiO4 is keyed to Mg but also
    depletes SiO."""
    assert group_dict["MgSiO3"]["gas_phase"] == "SiO"
    assert group_dict["Mg2SiO4"]["gas_phase"] == "Mg"
    assert group_dict["Mg2SiO4"]["stoich"] == {"Mg": 2, "SiO": 1}
    assert group_dict["MgSiO3"]["stoich"] == {"SiO": 1, "Mg": 1}
