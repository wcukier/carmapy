"""Integration-test fixtures.

Sets CARMAPY_EXE to the test binary (compiled with -DCARMA_TEST_CHECKS) for
the duration of each integration test, then restores the previous value.
Pathway and regression tests are not affected — they use the production binary.
"""
import os
import pytest


@pytest.fixture(autouse=True)
def _use_test_binary():
    """Temporarily point CARMAPY_EXE at the test binary for integration tests."""
    test_exe = os.environ.get("_CARMAPY_TEST_EXE")
    if not test_exe:
        yield
        return
    prev = os.environ.get("CARMAPY_EXE")
    os.environ["CARMAPY_EXE"] = test_exe
    try:
        yield
    finally:
        if prev is None:
            os.environ.pop("CARMAPY_EXE", None)
        else:
            os.environ["CARMAPY_EXE"] = prev
