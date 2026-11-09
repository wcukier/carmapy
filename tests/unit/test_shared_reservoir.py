"""Condensates sharing gas reservoirs: Mg2SiO4 (keyed to Mg), MgSiO3 on TiO2 and
homogeneous SiO (both keyed to SiO). subprocess is mocked, so no Fortran binary
runs.
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
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)
    c.add_hom_group("TiO2", 1e-8)
    r_het = 1e-8 * 2 ** (1 / 3)
    c.add_het_group("Mg2SiO4", "TiO2", r_het)
    c.add_het_group("MgSiO3", "TiO2", r_het)
    c.add_hom_group("SiO", 1e-8)
    c.set_nmr({"TiO2": 1e-7, "Mg": 1e-5, "SiO": 1e-5})
    c.set_stepping(dt=100, output_gap=5, n_tstep=10)
    return c


@pytest.mark.unit
def test_shared_gas_deduplicated():
    """MgSiO3 and SiO reuse one SiO reservoir; Mg2SiO4 gets its own Mg."""
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    r_het = 1e-8 * 2 ** (1 / 3)
    c.add_het_group("Mg2SiO4", "TiO2", r_het)
    c.add_het_group("MgSiO3", "TiO2", r_het)
    c.add_hom_group("SiO", 1e-8)

    # H2O is the always-present background gas.
    assert sorted(c.gases) == ["H2O", "Mg", "SiO", "TiO2"]
    assert c.groups["MgSiO3 on TiO2"].gas is c.groups["Pure SiO"].gas
    assert c.groups["Mg2SiO4 on TiO2"].gas.name == "Mg"


@pytest.mark.unit
def test_growth_and_depletions_follow_stofact():
    """Each condensate gets one growth pathway from its key gas, and a
    depletion for every other gas in its stofact."""
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    r_het = 1e-8 * 2 ** (1 / 3)
    c.add_het_group("Mg2SiO4", "TiO2", r_het)
    c.add_het_group("MgSiO3", "TiO2", r_het)
    c.add_hom_group("SiO", 1e-8)

    growth = {(g.elem.group.material, g.gas.name) for g in c.growth}
    assert growth == {("TiO2", "TiO2"), ("Mg2SiO4", "Mg"),
                      ("MgSiO3", "SiO"), ("SiO", "SiO")}
    # Mg2SiO4 limited by Mg but also depleting SiO; MgSiO3 limited by SiO but
    # also depleting Mg.
    depletions = {(d.group.material, d.gas.name): d.nu for d in c.depletions}
    assert depletions == {("Mg2SiO4", "SiO"): 1, ("MgSiO3", "Mg"): 1}


@pytest.mark.unit
def test_input_files_split_growth_and_depletion(tmp_path):
    """growth.txt carries (ielem, igas) per key pathway and depletion.txt
    (igroup, igas, nu, xref) per secondary reactant."""
    import f90nml

    c = _silicate_carma(str(tmp_path / "silicates"))
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c.run(suppress_output=True)

    def rows(fname):
        with open(os.path.join(c.name, "inputs", fname)) as f:
            header = f.readline().split()
            return header, [line.split() for line in f if line.strip()]

    header, growth = rows("growth.txt")
    assert header == ["ielem", "igas"]
    assert len(growth) == 4
    assert len({int(r[0]) for r in growth}) == 4

    header, depletion = rows("depletion.txt")
    assert header == ["igroup", "igas", "nu", "xref"]
    assert sorted(int(r[2]) for r in depletion) == [1, 1]
    # Each secondary row carries its gas's bottom abundance as a mass mixing
    # ratio (nmr 1e-5 in this deck).
    gas_by_igas = {g.igas: g for g in c.gases.values()}
    for r in depletion:
        g = gas_by_igas[int(r[1])]
        assert float(r[3]) == pytest.approx(1e-5 * g.wtmol_dif / c.wt_mol, rel=1e-8)

    nml = f90nml.read(os.path.join(c.name, "inputs", "input.nml"))
    assert nml["input_params"]["NGROWTH"] == 4
    assert nml["input_params"]["NDEPLETION"] == 2


@pytest.mark.unit
def test_shared_reservoir_roundtrip(tmp_path):
    """run() -> load_carma() -> run() reproduces identical input files for the
    shared-reservoir configuration."""
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
                  "nucleation.txt", "growth.txt", "depletion.txt",
                  "coagulation.txt", "gas_input.txt"]:
        f1 = os.path.join(dir1, "inputs", fname)
        f2 = os.path.join(dir2, "inputs", fname)
        assert filecmp.cmp(f1, f2, shallow=False), \
            f"{fname} differs after load_carma roundtrip"


@pytest.mark.unit
def test_cloud_base_uses_deepest_condensate(tmp_path):
    """SiO gas is the key gas of MgSiO3 and of the SiO condensate, and a
    secondary gas of Mg2SiO4; its reservoir abundance is evaluated at the deepest
    (highest-pressure) of the three cloud bases, here Mg2SiO4's."""
    import carmapy.chemistry as chemistry

    c = _silicate_carma(str(tmp_path / "base"))

    bases = {
        "Pure TiO2":       (9e6, 1900.0),
        "Mg2SiO4 on TiO2": (6e6, 1750.0),   # deepest SiO consumer
        "MgSiO3 on TiO2":  (5e6, 1700.0),
        "Pure SiO":        (1e6, 1500.0),
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
    assert calls[sio_hill] == (1750.0, 6e6)


@pytest.mark.unit
def test_run_without_depletion_file_loads(tmp_path):
    """A run written on main has a growth.txt but no depletion.txt or
    depletion_file entry; it loads with no secondary reactants."""
    import f90nml
    from carmapy.base import load_carma

    dir1 = str(tmp_path / "main")
    c1 = _silicate_carma(dir1)
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c1.run(suppress_output=True)

    os.remove(os.path.join(dir1, "inputs", "depletion.txt"))
    nml_path = os.path.join(dir1, "inputs", "input.nml")
    nml = f90nml.read(nml_path)
    del nml["io_files"]["depletion_file"]
    del nml["input_params"]["NDEPLETION"]
    nml.write(nml_path, force=True)

    c2 = load_carma(dir1)
    assert len(c2.growth) == len(c1.growth)
    assert c2.depletions == []


@pytest.mark.unit
def test_add_gas_depletion():
    """add_gas_depletion creates a missing reservoir, and refuses the group's
    key gas and a non-positive nu."""
    c = Carma("_t")
    c.do_secondary_depletion = False
    c.add_hom_group("TiO2", 1e-8)
    c.add_het_group("MgSiO3", "TiO2", 1e-8 * 2 ** (1 / 3))
    assert c.depletions == []

    d = c.add_gas_depletion("MgSiO3 on TiO2", "Mg", 1, xref=3e-5)
    assert d.gas is c.gases["Mg"]
    assert (d.nu, d.xref) == (1, 3e-5)
    assert c.depletions == [d]

    with pytest.raises(ValueError):
        c.add_gas_depletion("MgSiO3 on TiO2", "SiO", 1)
    with pytest.raises(ValueError):
        c.add_gas_depletion("MgSiO3 on TiO2", "Mg", 0)
    with pytest.raises(ValueError):
        c.add_gas_depletion("MgSiO3 on TiO2", "Mg", 1, xref=-1.0)


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


@pytest.mark.unit
def test_stofact_int_or_dict():
    """An int stofact is the key gas's factor and adds no depletions; a dict
    adds one for each non-key gas, and must include the key gas."""
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    r_het = 1e-8 * 2 ** (1 / 3)
    g = c.add_het_group("MgSiO3", "TiO2", r_het, stofact=1)
    assert g.nu_key == 1
    assert c.depletions == []

    g = c.add_het_group("Mg2SiO4", "TiO2", r_het)
    assert g.nu_key == 2
    assert [(d.gas.name, d.nu) for d in c.depletions] == [("SiO", 1)]

    with pytest.raises(ValueError):
        c.add_hom_group("Fe", 1e-8, stofact={"Mg": 1})
    with pytest.raises(ValueError):
        c.add_hom_group("Fe", 1e-8, stofact=0)
