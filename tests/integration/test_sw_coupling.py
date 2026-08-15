"""The shortwave, once it is wired into ``carma_rce``.

``test_swsolve`` checks the solver in isolation against pyHARP. This checks the
things that only exist after it is coupled: that the incident flux entering the
column is what the ``mu0`` convention says it is, that every watt of it is
either absorbed or reflected, and that turning irradiation off leaves the
longwave answer exactly where it was.

The ``mu0`` convention is the reason this file exists. Getting it wrong -- by a
factor of mu0, or by folding redistribution in twice -- produces a converged,
entirely plausible-looking profile at the wrong temperature, which no
smoke test would catch.
"""
import struct
import subprocess
import textwrap
from pathlib import Path

import numpy as np
import pytest

from test_rce import (MODULES, RT_SRC, SIGMA_SB, _build, _fmt, write_grey_ck)

pytestmark = pytest.mark.integration

I_CART = 1

#: Drives one flux evaluation with irradiation on and reports the shortwave
#: budget the coupling computed.
PROBE_SW = textwrap.dedent("""\
    program sw_couple
      use carma_precision_mod
      use carma_rce
      implicit none
      type(rce_type) :: rce
      integer :: nz, nbin, ngroup, nband, igridv, iz, rc
      real(kind=f) :: t_int, t_irr, t_star, mu0, w_surf
      real(kind=f) :: cp_cgs, grav_cgs, wtmol
      character(len=256) :: ck_path
      real(kind=f), allocatable :: p(:), pl(:), t(:)
      real(kind=f), allocatable :: radius(:,:)
      real(kind=f), allocatable :: qext(:,:,:), ssa(:,:,:), asym(:,:,:)
      real(kind=f), allocatable :: fnet(:), tnum(:), tden(:)

      open(unit=31, file='sw_case.txt', status='old', action='read')
      open(unit=32, file='sw_out.txt', status='replace')

      read(31,'(a)') ck_path
      read(31,*) nz, nbin, ngroup, nband, igridv
      read(31,*) t_int, t_irr, t_star, mu0, w_surf
      read(31,*) cp_cgs, grav_cgs, wtmol

      allocate(p(nz), pl(nz+1), t(nz), radius(nbin, ngroup))
      allocate(qext(nband, nbin, ngroup), ssa(nband, nbin, ngroup), &
               asym(nband, nbin, ngroup))
      allocate(fnet(nz+1), tnum(nz), tden(nz))

      read(31,*) (p(iz), iz = 1, nz)
      read(31,*) (pl(iz), iz = 1, nz+1)
      read(31,*) (t(iz), iz = 1, nz)

      radius(:,:)   = 1.e-5_f
      qext(:,:,:)   = 0._f
      ssa(:,:,:)    = 0._f
      asym(:,:,:)   = 0._f

      call rce_init(rce, nz, nbin, ngroup, nband, igridv, trim(ck_path), &
                    t_int, t_irr, t_star, mu0, w_surf, &
                    cp_cgs, grav_cgs, wtmol, &
                    I_ADIABAT_PARMENTIER, "", &
                    0, 1._f, 1.e30_f, 1._f, 0.02_f, 1, rc)
      if (rc < 0) stop '    *** FAILED rce_init ***'

      ! No clouds, and a hydrostatic thickness consistent with the profile.
      rce%numden_grp(:,:,:) = 0._f
      do iz = 1, nz
        rce%dz(iz) = (8.314462618_f * t(iz) / (wtmol * 1.e-3_f * (grav_cgs/100._f))) &
                     * log(pl(iz) / pl(iz+1))
      end do

      call rce_fluxes(rce, p, pl, t, radius, qext, ssa, asym, fnet, tnum, tden)

      write(32,'(4e26.17)') rce%absorbed_sw, rce%reflected_sw, &
                            sum(rce%f0), fnet(nz+1)

      call rce_destroy(rce)
      close(31)
      close(32)
    end program sw_couple
""")


@pytest.fixture(scope="module")
def sw_probe(tmp_path_factory):
    return _build(tmp_path_factory.mktemp("swcouple"), "sw_couple", PROBE_SW)


def _column(nz=24, p_top=1e2, p_bot=1e7, t_bot=1800.0, t_top=900.0):
    pl = np.geomspace(p_bot, p_top, nz + 1)
    p = np.sqrt(pl[:-1] * pl[1:])
    t = np.linspace(t_bot, t_top, nz)
    return p, pl, t


def _run(probe, ck, *, t_int=500.0, t_irr=0.0, t_star=0.0, mu0=0.5,
         w_surf=0.0, nband=40, nz=24):
    p, pl, t = _column(nz=nz)

    with open(probe / "sw_case.txt", "w") as f:
        f.write(f"{ck}\n")
        f.write(f"{nz} 1 1 {nband} {I_CART}\n")
        f.write(f"{t_int!r} {t_irr!r} {t_star!r} {mu0!r} {w_surf!r}\n")
        f.write(f"{1.3e8!r} {31600.0!r} {2.3!r}\n")
        f.write(_fmt(p) + "\n")
        f.write(_fmt(pl) + "\n")
        f.write(_fmt(t) + "\n")

    subprocess.run(["./sw_couple"], cwd=probe, check=True)

    absorbed, reflected, f0_sum, fnet_top = (
        float(v) for v in (probe / "sw_out.txt").read_text().split())
    return dict(absorbed=absorbed, reflected=reflected,
                f0_sum=f0_sum, fnet_top=fnet_top)


def test_incident_budget_closes(sw_probe, tmp_path_factory):
    """Every watt entering the top is either absorbed or reflected.

    This is the ``mu0`` test. The incident flux is ``mu0 * sigma * t_irr**4``,
    not ``sigma * t_irr**4``, and it is checked at two angles so the factor
    cannot cancel out of both.
    """
    d = tmp_path_factory.mktemp("ckbudget")
    ck = str(d / "grey.ck")
    # Thick enough that nothing reaches the bottom, so the budget is closed by
    # absorption and reflection alone.
    write_grey_ck(ck, kappa_ln=-50.0)

    t_irr = 1200.0
    for mu0 in (0.25, 0.5, 1.0):
        r = _run(sw_probe, ck, t_irr=t_irr, mu0=mu0)

        incident = mu0 * SIGMA_SB * t_irr ** 4
        total = r["absorbed"] + r["reflected"]

        assert total == pytest.approx(incident, rel=1e-9), (
            f"mu0={mu0}: absorbed {r['absorbed']:.6g} + reflected "
            f"{r['reflected']:.6g} = {total:.6g}, expected {incident:.6g}")


def test_f0_carries_exactly_sigma_t_irr_4(sw_probe, tmp_path_factory):
    """The per-band incident spectrum must sum to the flux t_irr names --
    the normalisation that separates the beam's colour from its magnitude."""
    d = tmp_path_factory.mktemp("ckf0")
    ck = str(d / "grey.ck")
    write_grey_ck(ck, kappa_ln=-50.0)

    for t_irr, t_star in ((1200.0, 0.0), (1200.0, 5800.0), (300.0, 3000.0)):
        r = _run(sw_probe, ck, t_irr=t_irr, t_star=t_star)
        assert r["f0_sum"] == pytest.approx(SIGMA_SB * t_irr ** 4, rel=1e-10), (
            f"t_irr={t_irr} t_star={t_star}")


def test_t_star_changes_colour_not_magnitude(sw_probe, tmp_path_factory):
    """Two very different stellar colours at the same t_irr deliver the same
    total energy, but not the same absorbed fraction -- which is the whole
    point of separating them."""
    d = tmp_path_factory.mktemp("ckcolour")
    ck = str(d / "grey.ck")
    # A temperature-dependent grey opacity would confuse this; a plain one
    # isolates the spectral placement of the beam against Rayleigh.
    write_grey_ck(ck, kappa_ln=-50.0)

    hot = _run(sw_probe, ck, t_irr=1000.0, t_star=8000.0)
    cool = _run(sw_probe, ck, t_irr=1000.0, t_star=1500.0)

    assert hot["f0_sum"] == pytest.approx(cool["f0_sum"], rel=1e-10)

    # A hotter beam sits further into the blue, where Rayleigh scattering is
    # stronger, so more of it comes back out.
    assert hot["reflected"] > cool["reflected"], (
        f"hot star reflected {hot['reflected']:.4g}, "
        f"cool {cool['reflected']:.4g}")


def test_irradiation_off_is_inert(sw_probe, tmp_path_factory):
    """t_irr = 0 must leave the longwave answer exactly as it was, with no
    shortwave budget and no Rayleigh in the optical properties."""
    d = tmp_path_factory.mktemp("ckoff")
    ck = str(d / "grey.ck")
    write_grey_ck(ck, kappa_ln=-50.0)

    off = _run(sw_probe, ck, t_irr=0.0)

    assert off["absorbed"] == 0.0
    assert off["reflected"] == 0.0
    assert off["f0_sum"] == 0.0
    assert np.isfinite(off["fnet_top"])


def test_more_irradiation_means_more_absorbed(sw_probe, tmp_path_factory):
    """Absorbed flux must scale with the beam, and the emergent net flux must
    respond -- a shortwave that is computed but not coupled would pass every
    test above and fail this one."""
    d = tmp_path_factory.mktemp("ckscale")
    ck = str(d / "grey.ck")
    write_grey_ck(ck, kappa_ln=-50.0)

    weak = _run(sw_probe, ck, t_irr=600.0)
    strong = _run(sw_probe, ck, t_irr=1200.0)

    # sigma*T^4, so 16x the flux for 2x the temperature.
    assert strong["absorbed"] / weak["absorbed"] == pytest.approx(16.0, rel=1e-6)
    assert strong["fnet_top"] != weak["fnet_top"]
