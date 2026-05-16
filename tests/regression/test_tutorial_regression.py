"""Full tutorial regression tests.

Runs the standard my_first_carma tutorial and compares numerical output against
the committed reference files in docs/source/notebooks/my_first_carma/.

These tests are very slow (~25 min). Enable with CARMAPY_RUN_SLOW=1.
"""
import os
import filecmp
import pytest
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent.parent
REF_DIR = REPO_ROOT / "docs" / "source" / "notebooks" / "my_first_carma"

RTOL = 1e-5
ATOL = 1e-30

DETERMINISTIC_INPUT_FILES = [
    "groups.txt",
    "gases.txt",
    "elements.txt",
    "nucleation.txt",
    "growth.txt",
    "coagulation.txt",
]


@pytest.fixture(scope="module")
def tutorial_results(tmp_path_factory):
    """Run the full my_first_carma tutorial in a temp directory."""
    from carmapy.example import example_carma
    from carmapy.results import Results

    if not (REF_DIR / "my_first_carma.txt").exists():
        pytest.skip("Reference tutorial output not found in docs/")

    tmp = tmp_path_factory.mktemp("regression")
    carma = example_carma(str(tmp / "my_first_carma"))
    carma.run(suppress_output=True)
    results = Results(carma, read_diag=False)
    return carma, results


@pytest.fixture(scope="module")
def ref_results(tutorial_results):
    """Load the committed reference results for comparison."""
    from carmapy.results import Results
    carma, _ = tutorial_results
    # Reconstruct a carma object pointing at the reference directory
    import copy
    ref_carma = copy.deepcopy(carma)
    ref_carma.name = str(REF_DIR)
    return Results(ref_carma, read_diag=False)


@pytest.mark.regression
@pytest.mark.slow
def test_tutorial_creates_output_file(tutorial_results):
    carma, _ = tutorial_results
    name = os.path.basename(carma.name)
    assert os.path.exists(os.path.join(carma.name, f"{name}.txt"))


@pytest.mark.regression
@pytest.mark.slow
def test_tutorial_creates_flux_file(tutorial_results):
    carma, _ = tutorial_results
    name = os.path.basename(carma.name)
    assert os.path.exists(os.path.join(carma.name, f"flux_{name}.txt"))


@pytest.mark.regression
@pytest.mark.slow
def test_tutorial_numden_regression(tutorial_results, ref_results):
    _, new = tutorial_results
    np.testing.assert_allclose(
        new.numden[:, :, :, -1],
        ref_results.numden[:, :, :, -1],
        rtol=RTOL,
        atol=ATOL,
        err_msg="numden at final timestep differs from reference",
    )


@pytest.mark.regression
@pytest.mark.slow
def test_tutorial_gas_abund_regression(tutorial_results, ref_results):
    _, new = tutorial_results
    np.testing.assert_allclose(
        new.gas_abund[:, :, -1],
        ref_results.gas_abund[:, :, -1],
        rtol=RTOL,
        atol=ATOL,
        err_msg="gas_abund at final timestep differs from reference",
    )


@pytest.mark.regression
@pytest.mark.slow
@pytest.mark.parametrize("filename", DETERMINISTIC_INPUT_FILES)
def test_tutorial_input_files_match_reference(tutorial_results, filename):
    """Python-generated input files should be byte-identical to reference."""
    carma, _ = tutorial_results
    new_file = os.path.join(carma.name, "inputs", filename)
    ref_file = REF_DIR / "inputs" / filename
    if not ref_file.exists():
        pytest.skip(f"Reference input file {filename} not found")
    assert filecmp.cmp(new_file, str(ref_file), shallow=False), (
        f"{filename} differs from reference"
    )
