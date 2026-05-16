"""Tests for run() input-file generation.

subprocess.Popen and shutil.copy are mocked so the Fortran binary is never
executed or required. These tests only verify that run() writes the expected
files with the expected content.
"""
import os
import copy
import pytest
import numpy as np
from unittest.mock import patch, MagicMock
import f90nml


def _mock_popen():
    mp = MagicMock()
    mp.return_value.poll.return_value = 0
    mp.return_value.stdout.readline.return_value = b""
    mp.return_value.stdout.read.return_value = b""
    mp.return_value.__enter__ = MagicMock(return_value=mp.return_value)
    mp.return_value.__exit__ = MagicMock(return_value=False)
    return mp


@pytest.fixture
def ready_carma(tmp_path, example_levels):
    """A fully-configured Carma object ready to call run()."""
    from carmapy.carmapy import Carma
    P, T, kzz, mu = example_levels
    c = Carma(str(tmp_path / "run_test"))
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.set_stepping(dt=250, output_gap=10, n_tstep=100)
    c.add_gas("TiO2")
    c.add_hom_group("TiO2", 1e-8)
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)
    return c


@pytest.mark.unit
def test_run_creates_input_nml(ready_carma):
    c = ready_carma
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c.run(suppress_output=True)
    nml_path = os.path.join(c.name, "inputs", "input.nml")
    assert os.path.exists(nml_path)


@pytest.mark.unit
def test_run_nml_has_correct_dimensions(ready_carma):
    c = ready_carma
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c.run(suppress_output=True)
    nml = f90nml.read(os.path.join(c.name, "inputs", "input.nml"))
    params = nml["input_params"]
    assert params["NZ"] == c.NZ
    assert params["NELEM"] == len(c.elems)
    assert params["NGROUP"] == len(c.groups)
    assert params["NGAS"] == len(c.gases)
    assert params["NBIN"] == c.NBIN


@pytest.mark.unit
def test_run_nml_has_correct_stepping(ready_carma):
    c = ready_carma
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c.run(suppress_output=True)
    nml = f90nml.read(os.path.join(c.name, "inputs", "input.nml"))
    params = nml["input_params"]
    assert params["iskip"] == c.output_gap
    assert params["nstep"] == c.n_tstep
    assert params["dtime"] == c.dt


@pytest.mark.unit
@pytest.mark.parametrize("filename", [
    "inputs/input.nml",
    "inputs/groups.txt",
    "inputs/gases.txt",
    "inputs/elements.txt",
    "inputs/nucleation.txt",
    "inputs/growth.txt",
    "inputs/coagulation.txt",
    "inputs/centers.txt",
    "inputs/levels.txt",
    "inputs/temps.txt",
    "inputs/winds.txt",
    "inputs/gas_input.txt",
    "inputs/pbound.txt",
    "inputs/gbound.txt",
])
def test_run_creates_expected_files(ready_carma, filename):
    c = ready_carma
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c.run(suppress_output=True)
    assert os.path.exists(os.path.join(c.name, filename)), f"Missing {filename}"


@pytest.mark.unit
def test_run_gases_file_has_header(ready_carma):
    c = ready_carma
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c.run(suppress_output=True)
    with open(os.path.join(c.name, "inputs", "gases.txt")) as f:
        header = f.readline()
    for col in ["name", "wtmol", "rho_cond", "coldia", "vp_offset", "vp_tcoeff"]:
        assert col in header, f"Column '{col}' missing from gases.txt header"


@pytest.mark.unit
def test_run_centers_file_row_count(ready_carma):
    c = ready_carma
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c.run(suppress_output=True)
    with open(os.path.join(c.name, "inputs", "centers.txt")) as f:
        rows = f.readlines()
    assert len(rows) == c.NZ + 1  # 1 header + NZ data rows


@pytest.mark.unit
def test_run_levels_file_row_count(ready_carma):
    c = ready_carma
    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c.run(suppress_output=True)
    with open(os.path.join(c.name, "inputs", "levels.txt")) as f:
        rows = f.readlines()
    assert len(rows) == c.NZ + 2  # 1 header + (NZ+1) level rows


@pytest.mark.unit
def test_run_raises_without_wt_mol(tmp_path):
    from carmapy.carmapy import Carma
    c = Carma(str(tmp_path / "run_err"))
    c.surface_grav = 980.0
    # wt_mol not set
    with pytest.raises(RuntimeError):
        c.run()


@pytest.mark.unit
def test_run_raises_without_surface_grav(tmp_path):
    from carmapy.carmapy import Carma
    c = Carma(str(tmp_path / "run_err2"))
    c.wt_mol = 2.0
    # surface_grav not set
    with pytest.raises(RuntimeError):
        c.run()
