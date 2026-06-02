"""Unit tests for the user-defined (custom) condensate workflow.

Custom condensates are defined by:
  1. add_gas(name, wtmol_dif=..., hill_formula=..., gcomp=...)  -- gas-phase props
  2. add_het_group(name, seed, rmin, **condensate_props)        -- condensate props
     (or add_hom_group(name, rmin, **condensate_props))

The gas-phase properties live on the Gas; the condensate properties live on the
Group. These tests cover construction, the file-writing layer, and a load/run
roundtrip (no Fortran binary is executed — subprocess is mocked).
"""
import os
import filecmp
import numpy as np
import pytest
from unittest.mock import patch, MagicMock

from carmapy.carmapy import Carma, Gas, Group


# Properties of a made-up condensate "Al_oxide" (gas-phase reservoir "Al_oxide").
CUSTOM_GAS = dict(wtmol_dif=26.98, hill_formula="Al", gcomp=0)
CUSTOM_COND = dict(
    rho_cond=3.99, wtmol=101.926, coldia=3.825e-8,
    vp_offset=17.7, vp_tcoeff=45892.6, vp_metcoeff=1.66, vp_logpcoeff=0,
    surften_0=690, surften_slope=0, stofact=2, is_typeIII=True,
)


def _mock_popen():
    mp = MagicMock()
    mp.return_value.poll.return_value = 0
    mp.return_value.stdout.readline.return_value = b""
    mp.return_value.stdout.read.return_value = b""
    mp.return_value.__enter__ = MagicMock(return_value=mp.return_value)
    mp.return_value.__exit__ = MagicMock(return_value=False)
    return mp


@pytest.mark.unit
def test_add_gas_custom_gas_phase():
    """A custom gas takes its gas-phase properties from add_gas kwargs."""
    c = Carma("_t")
    gas = c.add_gas("Al_oxide", **CUSTOM_GAS)
    assert "Al_oxide" in c.gases
    assert gas.wtmol_dif == pytest.approx(26.98)
    assert gas.hill_formula == "Al"


@pytest.mark.unit
def test_add_het_group_custom_condensate_props_on_group():
    """add_het_group stores the condensate properties on the Group and reuses
    the pre-added custom gas."""
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    c.add_gas("Al_oxide", **CUSTOM_GAS)
    group = c.add_het_group("Al_oxide", "TiO2", 1e-8 * 2 ** (1 / 3),
                            mucos=0.724172, **CUSTOM_COND)

    assert group.material == "Al_oxide"
    assert group.name == "Al_oxide on TiO2"
    assert group.rho_cond == pytest.approx(3.99)
    assert group.wtmol == pytest.approx(101.926)
    assert group.stofact == 2
    assert group.is_typeIII is True
    # The gas was reused (not duplicated) and keeps its gas-phase properties.
    assert group.gas.name == "Al_oxide"
    assert group.gas.wtmol_dif == pytest.approx(26.98)
    assert sum(k == "Al_oxide" for k in c.gases) == 1


@pytest.mark.unit
def test_add_het_group_custom_gas_phase_kwarg():
    """gas_phase kwarg lets the condensate reuse a differently-named gas."""
    c = Carma("_t")
    c.add_hom_group("TiO2", 1e-8)
    c.add_gas("Al", wtmol_dif=26.98, hill_formula="Al")
    group = c.add_het_group("Al_oxide", "TiO2", 1e-8 * 2 ** (1 / 3),
                            mucos=0.724172, gas_phase="Al", **CUSTOM_COND)
    assert group.gas.name == "Al"
    assert "Al" in c.gases and "Al_oxide" not in c.gases


@pytest.mark.unit
def test_custom_condensate_writes_input_files(tmp_path):
    """run() should write the custom condensate's properties into groups.txt and
    its gas into gases.txt."""
    from carmapy.example import example_levels
    P, T, kzz, mu = example_levels()
    c = Carma(str(tmp_path / "custom"))
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.add_P(P); c.add_T(T); c.add_kzz(kzz); c.calculate_z(mu)
    c.add_hom_group("TiO2", 1e-8)
    c.add_gas("Al_oxide", **CUSTOM_GAS)
    c.add_het_group("Al_oxide", "TiO2", 1e-8 * 2 ** (1 / 3),
                    mucos=0.724172, **CUSTOM_COND)
    c.set_nmr({"Al_oxide": 1e-7, "TiO2": 1e-7})
    c.set_stepping(dt=100, output_gap=5, n_tstep=10)

    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c.run(suppress_output=True)

    with open(os.path.join(c.name, "inputs", "groups.txt")) as f:
        groups = f.read()
    assert "Al_oxide on TiO2" in groups
    # condensate density 3.99 appears in the row
    assert "3.99" in groups

    with open(os.path.join(c.name, "inputs", "gases.txt")) as f:
        gases = f.read()
    assert "Al_oxide Vapor" in gases


@pytest.mark.unit
def test_custom_condensate_roundtrip(tmp_path):
    """run() -> load_carma() -> run() reproduces identical input files for a
    simulation that uses a custom condensate."""
    from carmapy.base import load_carma
    from carmapy.example import example_levels

    P, T, kzz, mu = example_levels()

    dir1 = str(tmp_path / "first")
    c1 = Carma(dir1)
    c1.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c1.set_atmospheric_parameters_from_defaults("Pure H2")
    c1.add_P(P); c1.add_T(T); c1.add_kzz(kzz); c1.calculate_z(mu)
    c1.add_hom_group("TiO2", 1e-8)
    c1.add_gas("Al_oxide", **CUSTOM_GAS)
    c1.add_het_group("Al_oxide", "TiO2", 1e-8 * 2 ** (1 / 3),
                     mucos=0.724172, **CUSTOM_COND)
    c1.set_nmr({"Al_oxide": 1e-7, "TiO2": 1e-7})
    c1.set_stepping(dt=100, output_gap=5, n_tstep=10)
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c1.run(suppress_output=True)

    c2 = load_carma(dir1)
    dir2 = str(tmp_path / "second")
    c2.name = dir2
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c2.run(suppress_output=True)

    for fname in ["groups.txt", "gases.txt", "elements.txt",
                  "nucleation.txt", "growth.txt", "coagulation.txt"]:
        f1 = os.path.join(dir1, "inputs", fname)
        f2 = os.path.join(dir2, "inputs", fname)
        assert filecmp.cmp(f1, f2, shallow=False), \
            f"{fname} differs after load_carma roundtrip"
