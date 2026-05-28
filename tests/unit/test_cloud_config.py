import pytest
from carmapy.carmapy import Carma
from carmapy.constants import RHO_TIO2


@pytest.mark.unit
def test_add_gas_creates_entry():
    c = Carma("_t")
    c.add_gas("TiO2")
    assert "TiO2" in c.gases


@pytest.mark.unit
def test_add_gas_sets_rho_cond():
    c = Carma("_t")
    c.add_gas("TiO2")
    assert c.gases["TiO2"].rho_cond == pytest.approx(RHO_TIO2)


@pytest.mark.unit
def test_add_gas_idempotent():
    c = Carma("_t")
    c.add_gas("TiO2")
    c.add_gas("TiO2")
    assert len([k for k in c.gases if k == "TiO2"]) == 1


@pytest.mark.unit
def test_add_gas_igas_sequential():
    c = Carma("_t")
    # H2O is added at index 1 in __init__
    c.add_gas("TiO2")
    assert c.gases["TiO2"].igas == 2


@pytest.mark.unit
def test_add_hom_group_creates_group():
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    assert "Pure TiO2" in c.groups


@pytest.mark.unit
def test_add_hom_group_returns_group():
    c = Carma("_t")
    grp = c.add_hom_group("TiO2", 1e-8)
    assert grp.name == "Pure TiO2"


@pytest.mark.unit
def test_add_hom_group_adds_gas():
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    assert "TiO2" in c.gases


@pytest.mark.unit
def test_add_hom_group_creates_hom_nuc():
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    assert len(c.nucs) == 1
    assert c.nucs[0].is_het is False


@pytest.mark.unit
def test_add_het_group_creates_group():
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    c.add_het_group("Mg2SiO4", "TiO2", 1e-8)
    assert "Mg2SiO4 on TiO2" in c.groups


@pytest.mark.unit
def test_add_het_group_creates_het_nuc():
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    c.add_het_group("Mg2SiO4", "TiO2", 1e-8)
    het_nucs = [n for n in c.nucs if n.is_het]
    assert len(het_nucs) == 1


@pytest.mark.unit
def test_add_coag_unknown_group_raises():
    c = Carma("_t")
    with pytest.raises(ValueError):
        c.add_coag("NonExistentGroup")


@pytest.mark.unit
def test_add_coag_appends_coag():
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    c.add_coag("Pure TiO2")
    assert len(c.coags) == 1


@pytest.mark.unit
def test_igroup_indices_sequential():
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    c.add_het_group("Mg2SiO4", "TiO2", 1e-8)
    igroups = [g.igroup for g in c.groups.values()]
    assert igroups == list(range(1, len(igroups) + 1))


@pytest.mark.unit
def test_ielem_indices_sequential():
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    c.add_het_group("Mg2SiO4", "TiO2", 1e-8)
    ielems = [e.ielem for e in c.elems.values()]
    assert ielems == list(range(1, len(ielems) + 1))
