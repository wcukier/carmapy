"""``carma_linalg``'s dense solver, against numpy.

This build links no LAPACK, so the LU used by the implicit temperature step is
written out by hand. It is checked here on its own rather than only through the
physics that consumes it: a solver that is subtly wrong shows up in a coupled
run as a profile that drifts, which is indistinguishable from a dozen other
things. Here it is either right or it is not.
"""
import shutil
import subprocess
import textwrap
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).parent.parent.parent
RT_SRC = REPO_ROOT / "src" / "CARMA" / "source" / "base"
MODULES = ["carma_precision_mod", "carma_linalg"]

PROBE = textwrap.dedent("""\
    program linalg_probe
      use carma_precision_mod
      use carma_linalg
      implicit none
      integer :: n, nrhs, i, j, k
      logical :: ok
      real(kind=f), allocatable :: a(:,:), lu(:,:), b(:), x(:)
      integer, allocatable :: piv(:)

      open(unit=31, file='linalg_case.txt', status='old', action='read')
      open(unit=32, file='linalg_out.txt', status='replace')

      read(31,*) n, nrhs
      allocate(a(n,n), lu(n,n), b(n), x(n), piv(n))

      read(31,*) ((a(i,j), j = 1, n), i = 1, n)

      lu(:,:) = a(:,:)
      call lu_factor(n, lu, piv, ok)

      if (ok) then
        write(32,'(i2)') 1
      else
        write(32,'(i2)') 0
      end if

      do k = 1, nrhs
        read(31,*) (b(i), i = 1, n)
        if (ok) then
          call lu_solve(n, lu, piv, b, x)
        else
          x(:) = 0._f
        end if
        write(32,'(1e26.17)') (x(i), i = 1, n)
      end do

      close(31)
      close(32)
    end program linalg_probe
""")


@pytest.fixture(scope="module")
def probe(tmp_path_factory):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")

    build = tmp_path_factory.mktemp("linalg")
    for mod in MODULES:
        src = RT_SRC / f"{mod}.F90"
        assert src.exists(), f"missing {src}"
        subprocess.run(["gfortran", "-ffree-line-length-512", "-O2", "-c", str(src)],
                       cwd=build, check=True)

    (build / "linalg_probe.F90").write_text(PROBE)
    subprocess.run(["gfortran", "-ffree-line-length-512", "-O2", "-o", "linalg_probe",
                    "linalg_probe.F90"] + [f"{m}.o" for m in MODULES],
                   cwd=build, check=True)
    return build


def _solve(probe, A, rhs):
    """Factor A once and solve for every right-hand side. Returns (ok, X)."""
    A = np.asarray(A, dtype=float)
    rhs = np.atleast_2d(np.asarray(rhs, dtype=float))
    n = A.shape[0]

    with open(probe / "linalg_case.txt", "w") as f:
        f.write(f"{n} {rhs.shape[0]}\n")
        for row in A:
            f.write(" ".join(repr(float(v)) for v in row) + "\n")
        for b in rhs:
            f.write(" ".join(repr(float(v)) for v in b) + "\n")

    subprocess.run(["./linalg_probe"], cwd=probe, check=True)

    tok = (probe / "linalg_out.txt").read_text().split()
    ok = bool(int(tok[0]))
    x = np.array([float(v) for v in tok[1:]], dtype=float).reshape(rhs.shape[0], n)
    return ok, x


def test_random_matrices_match_numpy(probe):
    """The bread-and-butter case, over a range of sizes."""
    rng = np.random.default_rng(20260813)

    for n in (1, 2, 5, 17, 64, 105):        # 105 is the production column
        A = rng.normal(size=(n, n))
        b = rng.normal(size=n)

        ok, x = _solve(probe, A, b)
        assert ok, f"n={n} reported singular"
        np.testing.assert_allclose(x[0], np.linalg.solve(A, b), rtol=1e-10, atol=1e-12)


def test_pivoting_is_actually_needed(probe):
    """A zero on the diagonal is fatal without pivoting and harmless with it."""
    A = np.array([[0.0, 1.0], [1.0, 0.0]])
    b = np.array([3.0, 7.0])

    ok, x = _solve(probe, A, b)
    assert ok, "a permutation matrix is not singular"
    np.testing.assert_allclose(x[0], np.linalg.solve(A, b), rtol=1e-12, atol=0)


def test_multiple_right_hand_sides_reuse_one_factorisation(probe):
    """lu_factor / lu_solve are split so the expensive half happens once."""
    rng = np.random.default_rng(7)
    n = 24
    A = rng.normal(size=(n, n))
    rhs = rng.normal(size=(4, n))

    ok, x = _solve(probe, A, rhs)
    assert ok
    for k in range(rhs.shape[0]):
        np.testing.assert_allclose(x[k], np.linalg.solve(A, rhs[k]),
                                   rtol=1e-10, atol=1e-12)


def test_singular_matrix_is_reported_not_silently_wrong(probe):
    """The implicit step falls back to the diagonal on `ok = false`, so a
    singular system has to be detected rather than returning noise."""
    A = np.array([[1.0, 2.0, 3.0],
                  [2.0, 4.0, 6.0],      # exactly twice row 1
                  [1.0, 0.0, 1.0]])
    ok, _ = _solve(probe, A, np.ones(3))
    assert not ok


def test_matrix_shaped_like_the_implicit_step(probe):
    """(I/dt - J) with a dense J: diagonally dominant, non-symmetric, and with
    the mixed signs the measured Jacobian actually has."""
    rng = np.random.default_rng(99)
    n = 105

    J = rng.normal(scale=1e-4, size=(n, n))
    J[np.diag_indices(n)] = -np.abs(rng.normal(scale=1e-2, size=n))
    J[0, 0] = +4.8e-2          # a positive diagonal entry, as measured up top

    A = np.eye(n) / 5000.0 - J
    b = rng.normal(scale=1e-3, size=n)

    ok, x = _solve(probe, A, b)
    assert ok
    np.testing.assert_allclose(x[0], np.linalg.solve(A, b), rtol=1e-9, atol=1e-14)
