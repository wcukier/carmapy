"""``carma_rayleigh``'s cross-section, against published values.

Rayleigh scattering is the one piece of the shortwave that is new physics
rather than a port, so it has no reference implementation inside this project
to check against. It is checked here against measured/tabulated H2 and He
cross-sections instead, plus the two structural properties the rest of the
shortwave depends on: the lambda^-4 scaling, and that a mixture is the
number-fraction weighted sum.

The mixture case matters more than it looks. The only composition shipped is
solar H/He, so a single-species bug would be invisible in every production run
while still being wrong the moment anyone supplies a different mixture -- which
is exactly what Phase 8 is for.
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
MODULES = ["carma_precision_mod", "carma_rayleigh"]

PROBE = textwrap.dedent("""\
    program rayleigh_probe
      use carma_precision_mod
      use carma_rayleigh
      implicit none
      integer :: ncase, i
      real(kind=f) :: lam, xh2, xhe
      real(kind=f) :: x(RAY_NSPEC)

      open(unit=31, file='ray_case.txt', status='old', action='read')
      open(unit=32, file='ray_out.txt', status='replace')

      read(31,*) ncase
      do i = 1, ncase
        read(31,*) lam, xh2, xhe
        x(RAY_H2) = xh2
        x(RAY_HE) = xhe
        write(32,'(3e26.17)') ray_sigma(RAY_H2, lam), &
                              ray_sigma(RAY_HE, lam), &
                              ray_sigma_mix(x, lam)
      end do

      close(31)
      close(32)
    end program rayleigh_probe
""")


@pytest.fixture(scope="module")
def probe(tmp_path_factory):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")

    build = tmp_path_factory.mktemp("rayleigh")
    for mod in MODULES:
        src = RT_SRC / f"{mod}.F90"
        assert src.exists(), f"missing {src}"
        subprocess.run(["gfortran", "-ffree-line-length-512", "-O2", "-c", str(src)],
                       cwd=build, check=True)

    (build / "rayleigh_probe.F90").write_text(PROBE)
    subprocess.run(["gfortran", "-ffree-line-length-512", "-O2", "-o", "rayleigh_probe",
                    "rayleigh_probe.F90"] + [f"{m}.o" for m in MODULES],
                   cwd=build, check=True)
    return build


def _sigma(probe, cases):
    """cases: list of (wavelength_cm, x_H2, x_He). Returns (n, 3) of
    [sigma_H2, sigma_He, sigma_mix] in cm^2/molecule."""
    with open(probe / "ray_case.txt", "w") as f:
        f.write(f"{len(cases)}\n")
        for lam, xh2, xhe in cases:
            f.write(f"{lam!r} {xh2!r} {xhe!r}\n")

    subprocess.run(["./rayleigh_probe"], cwd=probe, check=True)

    tok = (probe / "ray_out.txt").read_text().split()
    return np.array([float(v) for v in tok], dtype=float).reshape(len(cases), 3)


def test_h2_against_published_cross_sections(probe):
    """H2 Rayleigh cross-section at optical wavelengths.

    Reference values from the Dalgarno & Williams (1962) series, which is the
    standard tabulation. The module deliberately uses the polarizability form
    instead -- it generalises to other species -- so a few percent of
    disagreement is expected and is the documented cost of that choice.
    """
    def dalgarno(lam_cm):
        a = lam_cm * 1e8                      # cm -> Angstrom
        return 8.14e-13 / a**4 + 1.28e-6 / a**6 + 1.61 / a**8

    lams = [0.4e-4, 0.5e-4, 0.7e-4, 1.0e-4]   # 0.4 - 1.0 micron
    got = _sigma(probe, [(l, 1.0, 0.0) for l in lams])[:, 0]

    for lam, s in zip(lams, got):
        ref = dalgarno(lam)
        assert s == pytest.approx(ref, rel=0.10), (
            f"H2 at {lam*1e4:.2f} um: got {s:.3e}, Dalgarno & Williams {ref:.3e}")


def test_helium_is_far_weaker_than_hydrogen(probe):
    """He's polarizability is ~1/4 of H2's, so its cross-section is ~1/16."""
    got = _sigma(probe, [(0.5e-4, 1.0, 0.0)])[0]
    ratio = got[1] / got[0]
    assert ratio == pytest.approx((0.205 / 0.802) ** 2, rel=1e-6)
    # Sanity against the absolute scale, ~8.8e-29 cm^2 at 500 nm.
    assert got[1] == pytest.approx(8.8e-29, rel=0.10)


def test_lambda_to_the_minus_fourth(probe):
    """Doubling the wavelength must cut the cross-section by exactly 16."""
    got = _sigma(probe, [(0.5e-4, 1.0, 0.0), (1.0e-4, 1.0, 0.0)])[:, 0]
    assert got[0] / got[1] == pytest.approx(16.0, rel=1e-12)


def test_mixture_is_the_weighted_sum(probe):
    """The composition sum, which the shipped H/He default alone would not
    distinguish from a hard-coded H2 cross-section."""
    x_h2, x_he = 0.86, 0.14
    got = _sigma(probe, [(0.5e-4, x_h2, x_he)])[0]
    s_h2, s_he, s_mix = got

    assert s_mix == pytest.approx(x_h2 * s_h2 + x_he * s_he, rel=1e-12)
    # And it must actually differ from pure H2, or the test proves nothing.
    assert abs(s_mix - s_h2) / s_h2 > 0.10


def test_pure_helium_is_not_hydrogen(probe):
    """A composition with no H2 at all must return helium's cross-section --
    the case a hard-coded H2 assumption would get badly wrong."""
    got = _sigma(probe, [(0.5e-4, 0.0, 1.0)])[0]
    assert got[2] == pytest.approx(got[1], rel=1e-12)


def test_degenerate_inputs_return_zero(probe):
    """A non-positive wavelength returns zero rather than dividing by it."""
    got = _sigma(probe, [(0.0, 1.0, 0.0)])[0]
    assert np.all(got == 0.0)
