import os
import shutil
import subprocess
import numpy as np
import pytest
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
MINI_DIR = Path(__file__).parent / "data" / "mini_output" / "mini"
REF_DIR = REPO_ROOT / "docs" / "source" / "notebooks" / "my_first_carma"


def _pkg_dir() -> Path:
    """Return the directory of the installed carmapy package."""
    import carmapy
    return Path(carmapy.__file__).parent


def _pkg_exe() -> Path:
    return _pkg_dir() / "carmapy.exe"


def _carma_src_root() -> Path:
    """Return src/CARMA/ relative to the installed package (works for editable installs)."""
    return _pkg_dir().parent / "CARMA"


def _test_build_dir() -> Path:
    return _carma_src_root() / "build" / "carma_test"


def _test_exe() -> Path:
    return _test_build_dir() / "carmapy.exe"


def _binary_available():
    return _pkg_exe().exists() or shutil.which("carmapy.exe") is not None


def _build_test_binary():
    """Build carmapy.exe with -DCARMA_TEST_CHECKS into the carma_test build dir.

    Requires an editable install so that src/CARMA/ is accessible.
    Returns the path to the built binary, or None if the build failed.
    Only rebuilds if the binary is missing or older than carma_carmapy.F90.
    """
    carma_src = _carma_src_root()
    src_f90 = carma_src / "tests" / "carma_carmapy.F90"
    test_exe = _test_exe()
    test_build = _test_build_dir()

    if not carma_src.exists():
        print(
            "\nWARNING: CARMA source not found at expected location for editable install. "
            "Cannot build test binary. Integration tests will use the production binary."
        )
        return None

    if (
        test_exe.exists()
        and src_f90.exists()
        and test_exe.stat().st_mtime >= src_f90.stat().st_mtime
    ):
        return str(test_exe)

    test_build.mkdir(parents=True, exist_ok=True)
    shutil.copy(carma_src / "Makefile", test_build / "Makefile")

    result = None
    for compiler in ("ifort", "gfortran"):
        try:
            result = subprocess.run(
                ["make", "all", f"FORTRAN={compiler}",
                 "EXTRA_FFLAGS=-DCARMA_TEST_CHECKS"],
                cwd=str(test_build),
                capture_output=True,
                text=True,
            )
            if result.returncode == 0 and test_exe.exists():
                return str(test_exe)
        except FileNotFoundError:
            continue

    print(
        "\nWARNING: Could not build test binary with -DCARMA_TEST_CHECKS. "
        "Integration tests will use the production binary (invariant checks disabled).\n"
        f"Build output:\n{result.stderr[-2000:] if result else '(no output)'}"
    )
    return None


def pytest_sessionstart(session):
    """Before any tests run, build the test binary if integration tests are enabled.

    The test binary (compiled with -DCARMA_TEST_CHECKS) is only injected for
    `integration`-marked tests. Pathway and regression tests deliberately use
    the production binary so that committed references stay stable across
    compiler/flag changes that don't affect physics.
    """
    if os.environ.get("CARMAPY_RUN_INTEGRATION", "0") != "1":
        return
    if "CARMAPY_EXE" in os.environ:
        return  # caller already pinned a binary; respect it
    exe = _build_test_binary()
    # Store the test binary path in a private variable; the autouse fixture in
    # tests/integration/conftest.py will inject it only for integration tests.
    if exe:
        os.environ["_CARMAPY_TEST_EXE"] = exe


def pytest_configure(config):
    config.addinivalue_line("markers", "unit: fast, no Fortran binary required")
    config.addinivalue_line("markers", "parsing: uses committed fixture output files")
    config.addinivalue_line("markers", "integration: requires carmapy.exe, ~seconds to 2 min")
    config.addinivalue_line("markers", "pathway: per-physics-pathway autoregression vs committed .npz references, requires carmapy.exe")
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

    skip_pathway = pytest.mark.skip(
        reason="Set CARMAPY_RUN_INTEGRATION=1 to run pathway autoregression tests"
    )

    for item in items:
        if "integration" in item.keywords:
            if not binary_ok:
                item.add_marker(skip_no_binary)
            elif not run_integration:
                item.add_marker(skip_integration)
        if "pathway" in item.keywords:
            if not binary_ok:
                item.add_marker(skip_no_binary)
            elif not run_integration:
                item.add_marker(skip_pathway)
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
