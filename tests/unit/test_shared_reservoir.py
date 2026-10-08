"""Unit tests for condensates sharing a gas reservoir.

The competitive-condensation configuration is Mg2SiO4 (keyed to Mg), MgSiO3
and SiO2 (both keyed to SiO): one gas reservoir must be shared by several
condensates, and growth.txt must carry a per-pathway stoichiometric
coefficient (v2 format: ielem igas nu is_key).  No Fortran binary is executed
-- subprocess is mocked.
"""
import os
import filecmp
import numpy as np
import pytest
from unittest.mock import patch, MagicMock

from carmapy.carmapy import Carma


def _mock_popen():
    mp = MagicMock()
    mp.return_value.poll.return_value = 0
    mp.return_value.stdout.readline.return_value = b""
    mp.return_value.stdout.read.return_value = b""
    mp.return_value.__enter__ = MagicMock(return_value=mp.return_value)
    mp.return_value.__exit__ = MagicMock(return_value=False)
    return mp


def _silicate_carma(name):
    from carmapy.example import example_levels
    P, T, kzz, mu = example_levels()
    c = Carma(name)
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.add_P(P); c.add_T(T); c.add_kzz(kzz); c.calculate_z(mu)
    c.add_hom_group("TiO2", 1e-8)
    r_het = 1e-8 * 2 ** (1 / 3)
    c.add_het_group("Mg2SiO4", "TiO2", r_het)
    c.add_het_group("MgSiO3", "TiO2", r_het)
    c.add_het_group("SiO2", "TiO2", r_het)
    c.set_nmr({"TiO2": 1e-7, "Mg": 1e-5, "SiO": 1e-5})
    c.set_stepping(dt=100, output_gap=5, n_tstep=10)
    return c


@pytest.mark.unit
def test_shared_gas_deduplicated():
    """MgSiO3 and SiO2 reuse one SiO reservoir; Mg2SiO4 gets its own Mg."""
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    r_het = 1e-8 * 2 ** (1 / 3)
    c.add_het_group("Mg2SiO4", "TiO2", r_het)
    c.add_het_group("MgSiO3", "TiO2", r_het)
    c.add_het_group("SiO2", "TiO2", r_het)

    # H2O is the always-present background gas.
    assert sorted(c.gases) == ["H2O", "Mg", "SiO", "TiO2"]
    assert c.groups["MgSiO3 on TiO2"].gas is c.groups["SiO2 on TiO2"].gas
    assert c.groups["Mg2SiO4 on TiO2"].gas.name == "Mg"


@pytest.mark.unit
def test_growth_pathways_share_igas():
    """Each condensate gets one key growth pathway; the SiO-keyed pathways
    point at the same Gas object and carry the group's stofact as nu."""
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    r_het = 1e-8 * 2 ** (1 / 3)
    c.add_het_group("Mg2SiO4", "TiO2", r_het)
    c.add_het_group("MgSiO3", "TiO2", r_het)
    c.add_het_group("SiO2", "TiO2", r_het)

    rows = {(g.elem.group.material, g.gas.name): (g.nu, g.is_key)
            for g in c.growth}
    # The FINESST worked example: SiO2 limited by SiO; Mg2SiO4 limited by Mg
    # but also depleting SiO; MgSiO3 limited by SiO but also depleting Mg.
    assert rows[("Mg2SiO4", "Mg")] == (2, True)
    assert rows[("Mg2SiO4", "SiO")] == (1, False)
    assert rows[("MgSiO3", "SiO")] == (1, True)
    assert rows[("MgSiO3", "Mg")] == (1, False)
    assert rows[("SiO2", "SiO")] == (1, True)
    key_rows = [g for g in c.growth if g.is_key]
    assert len(key_rows) == 4  # one per condensate incl. TiO2


@pytest.mark.unit
def test_growth_file_v3_format(tmp_path):
    """growth.txt is written in v3 format with nu, is_key and xref columns."""
    c = _silicate_carma(str(tmp_path / "silicates"))
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c.run(suppress_output=True)

    with open(os.path.join(c.name, "inputs", "growth.txt")) as f:
        header = f.readline().split()
        rows = [line.split() for line in f if line.strip()]

    assert header == ["ielem", "igas", "nu", "is_key", "xref"]
    assert all(len(r) == 5 for r in rows)
    # One key row per volatile element (TiO2, Mg2SiO4, MgSiO3, SiO2) plus the
    # two secondary depletion rows (Mg2SiO4->SiO, MgSiO3->Mg).
    assert len(rows) == 6
    key_rows = [r for r in rows if int(r[3]) == 1]
    sec_rows = [r for r in rows if int(r[3]) == 0]
    assert len(key_rows) == 4
    assert len(sec_rows) == 2
    # Each element has exactly one key pathway.
    assert len({int(r[0]) for r in key_rows}) == 4
    # Forsterite's key pathway carries nu = 2.
    assert sorted(int(r[2]) for r in key_rows) == [1, 1, 1, 2]
    assert sorted(int(r[2]) for r in sec_rows) == [1, 1]
    # Key rows carry no reference abundance; each secondary row carries its
    # gas's bottom abundance as a mass mixing ratio (nmr 1e-5 in this deck).
    assert all(float(r[4]) == 0.0 for r in key_rows)
    gas_by_igas = {g.igas: g for g in c.gases.values()}
    for r in sec_rows:
        g = gas_by_igas[int(r[1])]
        assert float(r[4]) == pytest.approx(1e-5 * g.wtmol_dif / c.wt_mol, rel=1e-8)


@pytest.mark.unit
def test_shared_reservoir_roundtrip(tmp_path):
    """run() -> load_carma() -> run() reproduces identical input files for the
    three-silicate configuration."""
    from carmapy.base import load_carma

    dir1 = str(tmp_path / "first")
    c1 = _silicate_carma(dir1)
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c1.run(suppress_output=True)

    c2 = load_carma(dir1)
    dir2 = str(tmp_path / "second")
    c2.name = dir2
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c2.run(suppress_output=True)

    for fname in ["groups.txt", "gases.txt", "elements.txt",
                  "nucleation.txt", "growth.txt", "coagulation.txt",
                  "gas_input.txt"]:
        f1 = os.path.join(dir1, "inputs", fname)
        f2 = os.path.join(dir2, "inputs", fname)
        assert filecmp.cmp(f1, f2, shallow=False), \
            f"{fname} differs after load_carma roundtrip"


@pytest.mark.unit
def test_cloud_base_uses_deepest_condensate(tmp_path):
    """With MgSiO3 and SiO2 sharing the SiO reservoir, the reservoir abundance
    is evaluated at the deeper (higher-pressure) of the two cloud bases."""
    import carmapy.chemistry as chemistry

    c = _silicate_carma(str(tmp_path / "base"))

    bases = {
        "Pure TiO2":       (2e6, 1900.0),
        "Mg2SiO4 on TiO2": (6e6, 1750.0),
        "MgSiO3 on TiO2":  (5e6, 1700.0),
        "SiO2 on TiO2":    (9e6, 1800.0),   # deepest SiO consumer
    }

    def fake_find_cloud_base(carma, group):
        return bases[group.name]

    calls = {}

    def fake_abundances(T, P, species, metallicity):
        calls[species[0]] = (float(T[0]), float(P[0]))
        return np.array([1e-6])

    with patch.object(chemistry, "find_cloud_base", fake_find_cloud_base), \
         patch.object(chemistry, "get_fastchem_abundances", fake_abundances):
        chemistry.populate_abundances_at_cloud_base(c)

    sio_hill = c.gases["SiO"].hill_formula
    assert calls[sio_hill] == (1800.0, 9e6)


@pytest.mark.unit
def test_legacy_growth_file_still_loads(tmp_path):
    """A 2-column (pre-v2) growth.txt loads with nu = group stofact, is_key
    True."""
    dir1 = str(tmp_path / "legacy")
    c1 = _silicate_carma(dir1)
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c1.run(suppress_output=True)

    growth_path = os.path.join(dir1, "inputs", "growth.txt")
    with open(growth_path) as f:
        f.readline()
        rows = [line.split() for line in f if line.strip()]
    with open(growth_path, "w") as f:
        f.write("ielem\tigas\n")
        for r in rows:
            f.write(f"{r[0]}\t {r[1]}\n")

    from carmapy.base import load_carma
    c2 = load_carma(dir1)
    # load_carma does not reconstruct group.material; key by element name.
    by_name = {g.elem.name.split(" ")[0]: g for g in c2.growth}
    assert by_name["Mg2SiO4"].nu == 2
    assert by_name["SiO2"].nu == 1
    assert all(g.is_key for g in c2.growth)


@pytest.mark.unit
def test_gas_boundary_dicts_not_aliased():
    """Setting one gas's boundary condition must not change any other gas's
    (each Gas needs its own copy of the boundary defaults)."""
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    r_het = 1e-8 * 2 ** (1 / 3)
    c.add_het_group("Mg2SiO4", "TiO2", r_het)
    c.add_het_group("MgSiO3", "TiO2", r_het)

    gases = list(c.gases.values())
    assert len({id(g.boundary) for g in gases}) == len(gases)

    c.set_gas_boundary(c.gases["Mg"], bot_conc=6.8e-5)
    c.set_gas_boundary(c.gases["SiO"], bot_conc=5.1e-5)
    assert c.gases["Mg"].boundary["bot_conc"] == 6.8e-5
    assert c.gases["SiO"].boundary["bot_conc"] == 5.1e-5
    assert c.gases["TiO2"].boundary.get("bot_conc") in (None, -1)
