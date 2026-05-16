import os
import shutil
import numpy as np
import pytest
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
MINI_DIR = Path(__file__).parent / "data" / "mini_output" / "mini"
REF_DIR = REPO_ROOT / "docs" / "source" / "notebooks" / "my_first_carma"


def _binary_available():
    pkg_exe = REPO_ROOT / "src" / "carmapy" / "carmapy.exe"
    return pkg_exe.exists() or shutil.which("carmapy.exe") is not None


def pytest_configure(config):
    config.addinivalue_line("markers", "unit: fast, no Fortran binary required")
    config.addinivalue_line("markers", "parsing: uses committed fixture output files")
    config.addinivalue_line("markers", "integration: requires carmapy.exe, ~seconds to 2 min")
    config.addinivalue_line("markers", "regression: full tutorial run, ~25 min")
    config.addinivalue_line("markers", "slow: alias for regression")


def pytest_collection_modifyitems(config, items):
    run_integration = os.environ.get("CARMAPY_RUN_INTEGRATION", "0") == "1"
    run_slow = os.environ.get("CARMAPY_RUN_SLOW", "0") == "1"
    binary_ok = _binary_available()

    skip_no_binary = pytest.mark.skip(reason="carmapy.exe not found")
    skip_integration = pytest.mark.skip(
        reason="Set CARMAPY_RUN_INTEGRATION=1 to run integration tests"
    )
    skip_slow = pytest.mark.skip(
        reason="Set CARMAPY_RUN_SLOW=1 to run slow regression tests"
    )

    for item in items:
        if "integration" in item.keywords:
            if not binary_ok:
                item.add_marker(skip_no_binary)
            elif not run_integration:
                item.add_marker(skip_integration)
        if "regression" in item.keywords or "slow" in item.keywords:
            if not binary_ok:
                item.add_marker(skip_no_binary)
            elif not run_slow:
                item.add_marker(skip_slow)


@pytest.fixture(scope="session")
def example_levels():
    from carmapy.example import example_levels as _load
    return _load(t=1000)


@pytest.fixture(scope="session")
def mini_carma():
    """Reconstruct the Carma object matching tests/data/mini_output/mini."""
    from carmapy.carmapy import Carma
    from carmapy.example import example_levels as _load

    if not (MINI_DIR / "mini.txt").exists():
        pytest.skip("Mini fixture not generated; run tests/generate_fixture.py first")

    P, T, kzz, mu = _load(t=1000)
    carma = Carma(str(MINI_DIR))
    carma.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    carma.set_atmospheric_parameters_from_defaults("Pure H2")
    carma.set_stepping(dt=250, output_gap=5, n_tstep=50)
    carma.add_gas("TiO2")
    carma.add_hom_group("TiO2", 1e-8)
    carma.add_P(P)
    carma.add_T(T)
    carma.add_kzz(kzz)
    carma.calculate_z(mu)
    return carma


@pytest.fixture(scope="session")
def mini_results(mini_carma):
    from carmapy.results import Results
    return Results(mini_carma, read_diag=False)
