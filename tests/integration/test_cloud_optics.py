"""Verify the Fortran cloud-optics sum and the optics file it reads.

``carma_cloudopt`` reproduces ``results.py::_get_cloud_opacities``, so it is
checked against a numpy transcription of that function's sums rather than
against a tolerance pulled from thin air -- both compute the same quantity from
the same inputs, so agreement should be at round-off.

Unlike ``test_radiation_port.py`` these need neither pyharp nor torch: the whole
point of Phase 2 is that the run path has no such dependency.  ``gfortran`` is
still required for the solver test; the file-format test is pure Python.
"""
import subprocess
import textwrap
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).parent.parent.parent
RT_SRC = REPO_ROOT / "src" / "CARMA" / "source" / "base"
MODULES = ["carma_precision_mod", "carma_cloudopt"]

PROBE_SRC = textwrap.dedent("""\
    program cloudopt_probe
      use carma_precision_mod
      use carma_cloudopt
      implicit none
      integer :: nz, nbin, nwave, ngroup, iz, ibin, igroup, iw
      real(kind=f), allocatable :: numden(:,:,:), radius(:,:), dz(:)
      real(kind=f), allocatable :: qext(:,:,:), ssa(:,:,:), asym(:,:,:)
      real(kind=f), allocatable :: tau(:,:), w0(:,:), gasym(:,:)

      open(unit=31, file='cloudopt_case.txt', status='old', action='read')
      open(unit=32, file='cloudopt_out.txt', status='replace')
      read(31,*) nz, nbin, nwave, ngroup

      allocate(numden(nz, nbin, ngroup), radius(nbin, ngroup), dz(nz))
      allocate(qext(nwave, nbin, ngroup), ssa(nwave, nbin, ngroup), &
               asym(nwave, nbin, ngroup))
      allocate(tau(nwave, nz), w0(nwave, nz), gasym(nwave, nz))

      read(31,*) (((numden(iz, ibin, igroup), iz = 1, nz), ibin = 1, nbin), igroup = 1, ngroup)
      read(31,*) ((radius(ibin, igroup), ibin = 1, nbin), igroup = 1, ngroup)
      read(31,*) (dz(iz), iz = 1, nz)
      read(31,*) (((qext(iw, ibin, igroup), iw = 1, nwave), ibin = 1, nbin), igroup = 1, ngroup)
      read(31,*) (((ssa(iw, ibin, igroup),  iw = 1, nwave), ibin = 1, nbin), igroup = 1, ngroup)
      read(31,*) (((asym(iw, ibin, igroup), iw = 1, nwave), ibin = 1, nbin), igroup = 1, ngroup)

      call cloud_optics_column(nz, nbin, nwave, ngroup, numden, radius, dz, &
                               qext, ssa, asym, tau, w0, gasym)

      write(32,'(1e26.17)') ((tau(iw, iz),   iw = 1, nwave), iz = 1, nz)
      write(32,'(1e26.17)') ((w0(iw, iz),    iw = 1, nwave), iz = 1, nz)
      write(32,'(1e26.17)') ((gasym(iw, iz), iw = 1, nwave), iz = 1, nz)

      close(31)
      close(32)
    end program cloudopt_probe
""")


@pytest.fixture(scope="module")
def probe(tmp_path_factory):
    """Compile carma_cloudopt and a driver that exercises the bin sum."""
    shutil = pytest.importorskip("shutil")
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")

    build = tmp_path_factory.mktemp("cloudopt")
    for mod in MODULES:
        src = RT_SRC / f"{mod}.F90"
        assert src.exists(), f"missing {src}"
        subprocess.run(
            ["gfortran", "-ffree-line-length-512", "-O2", "-c", str(src)],
            cwd=build, check=True)

    (build / "cloudopt_probe.F90").write_text(PROBE_SRC)
    subprocess.run(
        ["gfortran", "-ffree-line-length-512", "-O2", "-o", "cloudopt_probe",
         "cloudopt_probe.F90"] + [f"{m}.o" for m in MODULES],
        cwd=build, check=True)

    return build


def test_cloud_optics_matches_reference(probe):
    """carma_cloudopt's bin sum vs a numpy transcription of results.py."""
    rng = np.random.default_rng(7)
    nz, nbin, nwave, ngroup = 40, 30, 24, 2

    # Number densities spanning many orders of magnitude with most bins empty,
    # which is what a real cloud field looks like and what exercises the
    # zero-cross-section skip.
    numden = 10.0 ** rng.uniform(-20, 4, (nz, nbin, ngroup))
    numden[rng.random((nz, nbin, ngroup)) < 0.6] = 0.0
    numden[7, :, :] = 0.0            # one wholly cloud-free layer

    radius = np.empty((nbin, ngroup))
    for g in range(ngroup):
        radius[:, g] = 1e-8 * 2.0 ** (np.arange(nbin) / 3.0) * (1.0 + 0.3 * g)

    dz = 10.0 ** rng.uniform(3, 6, nz)
    qext = rng.uniform(1e-8, 4.0, (nwave, nbin, ngroup))
    ssa = rng.uniform(0.0, 1.0, (nwave, nbin, ngroup))
    asym = rng.uniform(-0.4, 0.95, (nwave, nbin, ngroup))

    with open(probe / "cloudopt_case.txt", "w") as f:
        f.write(f"{nz} {nbin} {nwave} {ngroup}\n")
        for arr in (numden, radius, dz, qext, ssa, asym):
            f.write(" ".join(repr(float(v))
                             for v in arr.ravel(order="F")) + "\n")

    subprocess.run(["./cloudopt_probe"], cwd=probe, check=True)

    vals = np.array([float(x) for x in
                     (probe / "cloudopt_out.txt").read_text().split()])
    n = nwave * nz
    got = {"tau": vals[:n].reshape(nwave, nz, order="F"),
           "w0": vals[n:2 * n].reshape(nwave, nz, order="F"),
           "g": vals[2 * n:].reshape(nwave, nz, order="F")}

    # results.py's weight_term (pi r^2 numden) and its three weighted sums.
    weight = np.pi * radius[None, :, :] ** 2 * numden
    bext = np.einsum("zbg,wbg->zw", weight, qext)
    bsca = np.einsum("zbg,wbg->zw", weight, qext * ssa)
    bg = np.einsum("zbg,wbg->zw", weight, qext * ssa * asym)

    ref = {
        "tau": (bext * dz[:, None]).T,
        "w0": np.divide(bsca, bext, out=np.zeros_like(bsca),
                        where=bext > 0).T,
        "g": np.divide(bg, bsca, out=np.zeros_like(bg), where=bsca > 0).T,
    }

    for name in ("tau", "w0", "g"):
        scale = max(np.abs(ref[name]).max(), 1e-300)
        rel = np.abs(got[name] - ref[name]).max() / scale
        assert rel < 1e-13, f"{name} mismatch ({rel:.2e})"

    # A cloud-free layer must contribute nothing at all, not merely something
    # small -- the solver divides by these.
    assert np.all(got["tau"][:, 7] == 0.0)
    assert np.all(got["w0"][:, 7] == 0.0)
    assert np.all(got["g"][:, 7] == 0.0)


def test_optics_file_ordering(tmp_path):
    """The optics file must be written in the order the Fortran reads it.

    carma_carmapy.F90 loops group, then wavelength, then bin, and stores into
    ``(NWAVE, NBIN, NGROUP)``, whereas the Mie tables are ``(nbin, nband, 3)``.
    A transpose here would silently put one band's optics on every bin, which
    is exactly the kind of error that produces plausible-looking output, so the
    layout is pinned rather than assumed.
    """
    from carmapy import radiation

    nbin, nband = 4, 3
    # Values encode their own indices, so a transposed write is unambiguous.
    tables = {}
    for igroup, name in enumerate(("Pure TiO2", "Mg2SiO4 on TiO2")):
        t = np.zeros((nbin, nband, 3))
        for ibin in range(nbin):
            for iband in range(nband):
                t[ibin, iband, 0] = 100 * igroup + 10 * iband + ibin + 1
                t[ibin, iband, 1] = t[ibin, iband, 0] / 2.0   # Q_sca
                t[ibin, iband, 2] = 0.01 * (ibin + 1)         # g
        tables[name] = t

    path = radiation.write_optics_file(str(tmp_path / "optics.txt"), tables)

    lines = Path(path).read_text().splitlines()
    assert lines[0].split() == ["2", str(nband), str(nbin)]

    rows = [ln.split("\t") for ln in lines[2:]]
    assert len(rows) == 2 * nband * nbin

    k = 0
    for igroup in range(2):
        for iband in range(nband):
            for ibin in range(nbin):
                ig, iw, ib, qext, ssa, asym = rows[k]
                assert (int(ig), int(iw), int(ib)) == (igroup + 1, iband + 1,
                                                       ibin + 1)
                expect = 100 * igroup + 10 * iband + ibin + 1
                assert float(qext) == pytest.approx(expect)
                # ssa is Q_sca/Q_ext, which is 1/2 by construction above.
                assert float(ssa) == pytest.approx(0.5)
                assert float(asym) == pytest.approx(0.01 * (ibin + 1))
                k += 1


def test_optics_file_zero_qext_gives_zero_ssa(tmp_path):
    """Bins where the Mie call underflowed must not produce a 0/0 nan."""
    from carmapy import radiation

    tables = {"Pure TiO2": np.zeros((2, 2, 3))}
    path = radiation.write_optics_file(str(tmp_path / "optics.txt"), tables)

    vals = np.genfromtxt(path, delimiter="\t", skip_header=2)
    assert np.all(np.isfinite(vals))
    assert np.all(vals[:, 3:] == 0.0)
