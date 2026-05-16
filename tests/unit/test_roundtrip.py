"""Roundtrip test: run() → load_carma() → run() should produce identical input files.

This test does not execute the Fortran binary (subprocess.Popen and shutil.copy
are mocked). It only verifies the Python file-writing layer is deterministic
across a load.
"""
import os
import filecmp
import pytest
import numpy as np
from unittest.mock import patch, MagicMock


DETERMINISTIC_FILES = [
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
    "inputs/input.nml",
]


def _mock_popen():
    mp = MagicMock()
    mp.return_value.poll.return_value = 0
    mp.return_value.stdout.readline.return_value = b""
    mp.return_value.stdout.read.return_value = b""
    mp.return_value.__enter__ = MagicMock(return_value=mp.return_value)
    mp.return_value.__exit__ = MagicMock(return_value=False)
    return mp


@pytest.mark.unit
@pytest.mark.xfail(
    strict=True,
    reason="load_carma has bugs: (1) igridv subtraction L61, (2) gases.txt unpack mismatch L99",
)
def test_load_then_run_reproduces_input_files(tmp_path, example_levels):
    """After load_carma(dir1), re-running into dir2 should produce identical inputs."""
    from carmapy.carmapy import Carma
    from carmapy.base import load_carma
    from carmapy.chemistry import populate_abundances_at_cloud_base

    P, T, kzz, mu = example_levels

    # ── First run: build Carma from scratch, write input files to dir1 ────────
    dir1 = str(tmp_path / "first")
    c1 = Carma(dir1)
    c1.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c1.set_atmospheric_parameters_from_defaults("Pure H2")
    c1.set_stepping(dt=250, output_gap=5, n_tstep=50)
    c1.add_gas("TiO2")
    c1.add_hom_group("TiO2", 1e-8)
    c1.add_P(P)
    c1.add_T(T)
    c1.add_kzz(kzz)
    c1.calculate_z(mu)
    populate_abundances_at_cloud_base(c1)

    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c1.run(suppress_output=True)

    # ── Roundtrip: load_carma, then re-run into a new directory ──────────────
    c2 = load_carma(dir1)
    dir2 = str(tmp_path / "second")
    c2.name = dir2

    with patch("subprocess.Popen", _mock_popen()), patch("shutil.copy"):
        c2.run(suppress_output=True)

    # ── Assertion: every input file should be byte-identical ─────────────────
    mismatched = []
    for rel in DETERMINISTIC_FILES:
        f1 = os.path.join(dir1, rel)
        f2 = os.path.join(dir2, rel)
        if not filecmp.cmp(f1, f2, shallow=False):
            mismatched.append(rel)
    assert not mismatched, f"Files differ after roundtrip: {mismatched}"
