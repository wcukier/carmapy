"""Verify the radiative-convective temperature update in ``carma_rce``.

Three of these are port checks: the Parmentier adiabat, the level-temperature
reconstruction and the hydrostatic regrid all have Python counterparts, so they
are compared against those rather than against a hand-picked tolerance.

The energy-conservation checks are the load-bearing ones.  The internal flux
enters only as the net flux imposed at the base of the column, and nothing
carries it to the top except the radiative solve and the convective
adjustment, so relaxing a grey column and confirming the emergent flux reaches
``sigma*Teff^4`` exercises both end to end.  It is a sharp test: the adjustment
conserves enthalpy exactly, so any leak in the handoff between the two shows up
directly in the emergent flux.

Needs ``gfortran`` and numpy only -- no pyharp, no torch.  Run with
``CARMAPY_RUN_INTEGRATION=1 pytest tests/integration/test_rce.py``.
"""
import subprocess
import textwrap
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).parent.parent.parent
RT_SRC = REPO_ROOT / "src" / "CARMA" / "source" / "base"
MODULES = ["carma_precision_mod", "carma_enums_mod", "carma_planck",
           "carma_ckopacity", "carma_cloudopt", "carma_rtsolve",
           "carma_swsolve", "carma_linalg", "carma_rayleigh", "carma_rce"]

# Must match carmapy.radiation._CK_MAGIC / _CK_VERSION.
CK_MAGIC = 0x43524D41434B3031
CK_VERSION = 1

SIGMA_SB = 5.670374419e-8
K_B = 1.380649e-23
M_PROT = 1.67262192369e-27
R_GAS = 8.31446261815324

# The Fortran holds these as parameters in carma_rce.F90; the Python side of
# the same pair lives in carmapy.constants. Importing rather than restating
# them keeps the port check honest.
from carmapy.constants import (PARMENTIER_A_COEFF as PARMENTIER_A,
                               PARMENTIER_B_COEFF as PARMENTIER_B)

I_CART = 1
I_LOGP = 8


# --------------------------------------------------------------------------
# probes
# --------------------------------------------------------------------------

PROBE_UNITS = textwrap.dedent("""\
    program rce_units
      use carma_precision_mod
      use carma_rce
      implicit none
      integer :: nz, iz, igridv
      real(kind=f) :: t0, p0, wtmol, grav
      real(kind=f), allocatable :: p(:), pl(:), t(:), tl(:), zc(:), zl(:)
      real(kind=f), allocatable :: padia(:), tadia(:)

      open(unit=31, file='units_case.txt', status='old', action='read')
      open(unit=32, file='units_out.txt', status='replace')

      ! --- adiabat ---
      read(31,*) nz, t0, p0
      allocate(padia(nz), tadia(nz))
      read(31,*) (padia(iz), iz = 1, nz)
      do iz = 1, nz
        tadia(iz) = rce_adiabat(padia(iz), t0, p0)
      end do
      write(32,'(1e26.17)') (tadia(iz), iz = 1, nz)
      deallocate(padia, tadia)

      ! --- level temperatures and regrid ---
      read(31,*) nz, igridv, wtmol, grav
      allocate(p(nz), pl(nz+1), t(nz), tl(nz+1), zc(nz), zl(nz+1))
      read(31,*) (p(iz),  iz = 1, nz)
      read(31,*) (pl(iz), iz = 1, nz+1)
      read(31,*) (t(iz),  iz = 1, nz)

      zc(:) = 0._f
      zl(:) = 0._f

      call rce_level_temps(nz, p, pl, t, tl)
      call rce_regrid_z(nz, igridv, pl, tl, wtmol, grav, zc, zl)

      write(32,'(1e26.17)') (tl(iz), iz = 1, nz+1)
      write(32,'(1e26.17)') (zl(iz), iz = 1, nz+1)
      write(32,'(1e26.17)') (zc(iz), iz = 1, nz)

      close(31)
      close(32)
    end program rce_units
""")

PROBE_ADIABAT = textwrap.dedent("""\
    program rce_adiabat_probe
      use carma_precision_mod
      use carma_rce
      implicit none
      type(rce_type) :: rce
      integer :: n, i, rc
      real(kind=f) :: t0, p0, x, y
      real(kind=f), allocatable :: pres(:), temp(:)
      character(len=256) :: table_path

      open(unit=31, file='adiabat_case.txt', status='old', action='read')
      open(unit=32, file='adiabat_out.txt', status='replace')

      read(31,'(A)') table_path
      read(31,*) n, t0, p0
      allocate(pres(n), temp(n))
      read(31,*) (pres(i), i = 1, n)

      rc = 0
      call rce_load_adiabat(trim(table_path), rce, rc)
      if (rc < 0) stop "    *** FAILED rce_load_adiabat ***"

      ! The gradient itself, so an interpolation difference is distinguishable
      ! from an integration difference.
      do i = 1, n
        temp(i) = rce_grad_ad(rce, t0, pres(i))
      end do
      write(32,'(1e26.17)') (temp(i), i = 1, n)

      ! The integrated adiabat, swept outward from the anchor the way
      ! rce_update rebuilds the convective zone.
      x = log(p0)
      y = log(t0)
      do i = 1, n
        y = rce_integrate_adiabat(rce, y, x, log(pres(i)))
        x = log(pres(i))
        temp(i) = exp(y)
      end do
      write(32,'(1e26.17)') (temp(i), i = 1, n)

      close(31)
      close(32)
      call rce_destroy(rce)
    end program rce_adiabat_probe
""")

PROBE_KZZ = textwrap.dedent("""\
    program rce_kzz_probe
      use carma_precision_mod
      use carma_rce
      implicit none
      type(rce_type) :: rce
      integer :: nz, i, rc, adiabat
      real(kind=f) :: grav_cgs, rescale
      real(kind=f), allocatable :: p(:), pl(:), t(:), fnet(:), kzz(:)
      character(len=256) :: table_path

      open(unit=31, file='kzz_case.txt', status='old', action='read')
      open(unit=32, file='kzz_out.txt', status='replace')

      read(31,'(A)') table_path
      read(31,*) nz, adiabat, rce%t_int, rce%wtmol, grav_cgs, rce%kzz_mixl_scale

      ! rce%grav is SI; the namelist and the Python reference are both cgs.
      rce%grav    = grav_cgs * 1.e-2_f
      rce%nz      = nz
      rce%adiabat = adiabat
      rce%sw_on   = .false.

      if (adiabat == I_ADIABAT_TABLE) then
        rc = 0
        call rce_load_adiabat(trim(table_path), rce, rc)
        if (rc < 0) stop "    *** FAILED rce_load_adiabat ***"
      end if

      allocate(p(nz), pl(nz+1), t(nz), fnet(nz+1), kzz(nz+1))
      read(31,*) (p(i), i = 1, nz)
      read(31,*) (pl(i), i = 1, nz+1)
      read(31,*) (t(i), i = 1, nz)
      read(31,*) (fnet(i), i = 1, nz+1)

      call rce_kzz(rce, p, pl, t, fnet, kzz, rescale)

      write(32,'(1e26.17)') (kzz(i), i = 1, nz+1)
      write(32,'(1e26.17)') rescale

      close(31)
      close(32)
      call rce_destroy(rce)
    end program rce_kzz_probe
""")

PROBE_RUN = textwrap.dedent("""\
    program rce_run
      use carma_precision_mod
      use carma_rce
      implicit none
      type(rce_type) :: rce
      integer :: nz, nbin, ngroup, nband, nelem, igridv, mode, gap_max
      integer :: nstep, istep, iz, rc
      real(kind=f) :: t_int, cp_cgs, grav_cgs, wtmol
      real(kind=f) :: accel, dt_max, dt_tol, dtau_tol, dtime
      character(len=256) :: ck_path
      real(kind=f), allocatable :: p(:), pl(:), t(:), zc(:), zl(:)
      real(kind=f), allocatable :: numden(:,:,:), radius(:,:)
      real(kind=f), allocatable :: qext(:,:,:), ssa(:,:,:), asym(:,:,:)
      integer, allocatable :: elem2group(:)
      logical, allocatable :: elem_is_number(:)

      open(unit=31, file='run_case.txt', status='old', action='read')
      open(unit=32, file='run_out.txt', status='replace')

      read(31,'(a)') ck_path
      read(31,*) nz, nbin, ngroup, nband, nelem, igridv, mode, gap_max, nstep
      read(31,*) t_int, cp_cgs, grav_cgs, wtmol
      read(31,*) accel, dt_max, dt_tol, dtau_tol, dtime

      allocate(p(nz), pl(nz+1), t(nz), zc(nz), zl(nz+1))
      allocate(numden(nz, nelem, nbin), radius(nbin, ngroup))
      allocate(qext(nband, nbin, ngroup), ssa(nband, nbin, ngroup), &
               asym(nband, nbin, ngroup))
      allocate(elem2group(nelem), elem_is_number(nelem))

      read(31,*) (p(iz),  iz = 1, nz)
      read(31,*) (pl(iz), iz = 1, nz+1)
      read(31,*) (t(iz),  iz = 1, nz)
      read(31,*) (zl(iz), iz = 1, nz+1)
      read(31,*) (zc(iz), iz = 1, nz)

      numden(:,:,:) = 0._f
      radius(:,:)   = 1.e-5_f
      qext(:,:,:)   = 0._f
      ssa(:,:,:)    = 0._f
      asym(:,:,:)   = 0._f
      elem2group(:)     = 1
      elem_is_number(:) = .true.

      rc = 0
      call rce_init(rce, nz, nbin, ngroup, nband, igridv, trim(ck_path), &
                    t_int, 0._f, 0._f, 0.5_f, 0._f, &
                    cp_cgs, grav_cgs, wtmol, &
                    I_ADIABAT_PARMENTIER, "", I_KZZ_STATIC, 1._f, &
                    mode, accel, dt_max, dt_tol, dtau_tol, gap_max, rc)
      if (rc < 0) stop "    *** FAILED rce_init ***"

      do istep = 1, nstep
        call rce_update(rce, dtime, p, pl, t, zc, zl, nelem, numden, &
                        elem2group, elem_is_number, radius, qext, ssa, asym)
      end do

      write(32,'(1e26.17)') (t(iz), iz = 1, nz)
      write(32,'(1e26.17)') (rce%fnet(iz), iz = 1, nz+1)
      write(32,'(1e26.17)') (rce%dtdt(iz), iz = 1, nz)
      write(32,'(1e26.17)') (rce%tau_rad(iz), iz = 1, nz)
      write(32,'(1e26.17)') (merge(1._f, 0._f, rce%is_conv(iz)), iz = 1, nz)
      write(32,'(1e26.17)') rce%conv_resid
      write(32,'(4i10)') rce%nsolve, rce%nclamp, rce%nz_rcb, rce%nzone

      call rce_destroy(rce)
      close(31)
      close(32)
    end program rce_run
""")


PROBE_JAC = textwrap.dedent("""\
    program rce_jac
      use carma_precision_mod
      use carma_rce
      implicit none
      type(rce_type) :: rce
      integer :: nz, nbin, ngroup, nband, igridv, iz, k, rc
      real(kind=f) :: t_int, cp_cgs, grav_cgs, wtmol, dmass
      character(len=256) :: ck_path
      real(kind=f), allocatable :: p(:), pl(:), t(:), zl(:), dt_in(:), tp(:)
      real(kind=f), allocatable :: radius(:,:)
      real(kind=f), allocatable :: qext(:,:,:), ssa(:,:,:), asym(:,:,:)
      real(kind=f), allocatable :: f1(:), tn(:), td(:), h0(:), h1(:), pred(:)

      open(unit=31, file='jac_case.txt', status='old', action='read')
      open(unit=32, file='jac_out.txt', status='replace')

      read(31,'(a)') ck_path
      read(31,*) nz, nbin, ngroup, nband, igridv
      read(31,*) t_int, cp_cgs, grav_cgs, wtmol

      allocate(p(nz), pl(nz+1), t(nz), zl(nz+1), dt_in(nz), tp(nz))
      allocate(radius(nbin, ngroup))
      allocate(qext(nband, nbin, ngroup), ssa(nband, nbin, ngroup), &
               asym(nband, nbin, ngroup))
      allocate(f1(nz+1), tn(nz), td(nz), h0(nz), h1(nz), pred(nz))

      read(31,*) (p(iz),  iz = 1, nz)
      read(31,*) (pl(iz), iz = 1, nz+1)
      read(31,*) (t(iz),  iz = 1, nz)
      read(31,*) (zl(iz), iz = 1, nz+1)
      read(31,*) (dt_in(iz), iz = 1, nz)

      radius(:,:) = 1.e-5_f
      qext(:,:,:) = 0._f
      ssa(:,:,:)  = 0._f
      asym(:,:,:) = 0._f

      rc = 0
      call rce_init(rce, nz, nbin, ngroup, nband, igridv, trim(ck_path), &
                    t_int, 0._f, 0._f, 0.5_f, 0._f, &
                    cp_cgs, grav_cgs, wtmol, &
                    I_ADIABAT_PARMENTIER, "", I_KZZ_STATIC, 1._f, &
                    0, 1._f, 1.e30_f, 1._f, 0.02_f, 1, rc)
      if (rc < 0) stop 'rce_init failed'

      do iz = 1, nz
        rce%dz(iz) = zl(iz+1) - zl(iz)
      end do
      rce%numden_grp(:,:,:) = 0._f

      call rce_optics_prep(rce, radius, qext, ssa, asym)

      call rce_solve(rce, p, pl, t)
      h0(:) = rce%dtdt(:)

      call rce_jacobian(rce, p, pl, t)

      ! The Jacobian's own prediction for a finite step, J*dT ...
      do iz = 1, nz
        pred(iz) = 0._f
        do k = 1, nz
          pred(iz) = pred(iz) + rce%jac(iz, k) * dt_in(k)
        end do
      end do

      ! ... against what the solver actually does for that same step.
      tp(:) = t(:) + dt_in(:)
      call rce_fluxes(rce, p, pl, tp, f1, tn, td)
      f1(1) = rce%fnet(1)
      do iz = 1, nz
        dmass  = (pl(iz) - pl(iz+1)) / (grav_cgs * 1.e-2_f)
        h1(iz) = -(f1(iz+1) - f1(iz)) / (dmass * cp_cgs * 1.e-4_f)
      end do

      write(32,'(1e26.17)') ((rce%jac(iz,k), iz = 1, nz), k = 1, nz)
      write(32,'(1e26.17)') (pred(iz), iz = 1, nz)
      write(32,'(1e26.17)') (h1(iz) - h0(iz), iz = 1, nz)

      call rce_destroy(rce)
      close(31)
      close(32)
    end program rce_jac
""")


def _build(tmp, probe_name, probe_src):
    """Compile the RT modules plus one probe program into ``tmp``."""
    for mod in MODULES:
        src = RT_SRC / f"{mod}.F90"
        assert src.exists(), f"missing {src}"
        subprocess.run(
            ["gfortran", "-ffree-line-length-512", "-O2", "-c", str(src)],
            cwd=tmp, check=True)

    (tmp / f"{probe_name}.F90").write_text(probe_src)
    subprocess.run(
        ["gfortran", "-ffree-line-length-512", "-O2", "-o", probe_name,
         f"{probe_name}.F90"] + [f"{m}.o" for m in MODULES],
        cwd=tmp, check=True)
    return tmp


@pytest.fixture(scope="module")
def units_probe(tmp_path_factory):
    import shutil
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")
    return _build(tmp_path_factory.mktemp("rce_units"), "rce_units",
                  PROBE_UNITS)


@pytest.fixture(scope="module")
def adiabat_probe(tmp_path_factory):
    import shutil
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")
    return _build(tmp_path_factory.mktemp("rce_adiabat"), "rce_adiabat_probe",
                  PROBE_ADIABAT)


@pytest.fixture(scope="module")
def kzz_probe(tmp_path_factory):
    import shutil
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")
    return _build(tmp_path_factory.mktemp("rce_kzz"), "rce_kzz_probe",
                  PROBE_KZZ)


@pytest.fixture(scope="module")
def run_probe(tmp_path_factory):
    import shutil
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")
    return _build(tmp_path_factory.mktemp("rce_run"), "rce_run", PROBE_RUN)


def _fmt(arr):
    return " ".join(repr(float(v)) for v in np.ravel(arr)) + "\n"


# --------------------------------------------------------------------------
# a synthetic grey opacity table
# --------------------------------------------------------------------------

def write_grey_ck(path, kappa_ln, nband=40, ng=4,
                  wn_lo=1.0, wn_hi=30000.0, dlnk_dlnt=0.0):
    """Write a ck table with a single wavelength-independent absorption cross
    section, in the byte layout ``carmapy.radiation.export_ck_table`` uses.

    Grey opacity makes the correlated-k machinery exact rather than
    approximate, so any discrepancy in the flux is the solver's, not the
    opacity's.  The bands span 1 to 30000 cm^-1, which holds better than
    99.99% of the Planck function between 300 and 3000 K -- the emergent flux
    can therefore be compared against the full ``sigma*T^4`` without a
    band-coverage correction.

    ``dlnk_dlnt`` makes the cross section vary as ``T**dlnk_dlnt``, and
    defaults to 0 -- constant in temperature, which is what the energy and
    cadence tests want.

    **Leaving it at 0 is what made a whole class of bug invisible.** With the
    opacity independent of temperature, perturbing one layer changes no other
    layer's absorption, so every off-diagonal term of the heating Jacobian is
    zero and the column decouples. Three successive stabilisers passed this
    suite and then failed, oscillated or produced NaN on a real 196-band
    column for reasons the suite structurally could not see. Any test about
    how layers couple has to set this non-zero; 2 is a reasonable stand-in for
    a molecular band.
    """
    edges = np.geomspace(wn_lo, wn_hi, nband + 1)
    wmin, wmax = edges[:-1], edges[1:]

    gauss_pts = (np.arange(ng) + 0.5) / ng
    gauss_wts = np.full(ng, 1.0 / ng)

    wavenumber = np.repeat(0.5 * (wmin + wmax), ng)
    weights = np.tile(gauss_wts, nband)

    pres = np.geomspace(1.0, 1e8, 12)
    temp = np.linspace(200.0, 4000.0, 10)

    # kappa is stored as ln(cross section), so a linear term in ln T here is a
    # power law in T.
    kappa = np.full((nband * ng, pres.size, temp.size), kappa_ln)
    if dlnk_dlnt != 0.0:
        kappa = kappa + dlnk_dlnt * np.log(temp / temp[0])[None, None, :]

    with open(path, "wb") as f:
        np.array([CK_MAGIC], dtype=np.int64).tofile(f)
        np.array([CK_VERSION, nband, ng, pres.size, temp.size],
                 dtype=np.int64).tofile(f)
        for arr in (wmin, wmax, gauss_pts, gauss_wts,
                    wavenumber, weights, pres, temp):
            arr.astype(np.float64).tofile(f)
        kappa.astype(np.float64).ravel(order="F").tofile(f)

    return nband


# --------------------------------------------------------------------------
# tests
# --------------------------------------------------------------------------

def test_adiabat_matches_python(units_probe, tmp_path_factory):
    """``rce_adiabat`` vs ``extend_atmosphere``'s Parmentier expressions."""
    p_adia = np.geomspace(1e3, 1e8, 25)
    t0, p0 = 1900.0, 1e7

    # rce_regrid_z takes gravity in m/s^2; only rce_init does the cgs
    # conversion. 316 m/s^2 is log g = 4.5 cgs, a typical brown dwarf.
    nz_g, igridv, wtmol, grav = 6, I_CART, 2.3, 316.0
    p = np.geomspace(1e7, 1e2, nz_g)
    pl = np.geomspace(2e7, 5e1, nz_g + 1)
    t = np.linspace(1900.0, 700.0, nz_g)

    with open(units_probe / "units_case.txt", "w") as f:
        f.write(f"{len(p_adia)} {t0!r} {p0!r}\n")
        f.write(_fmt(p_adia))
        f.write(f"{nz_g} {igridv} {wtmol!r} {grav!r}\n")
        f.write(_fmt(p))
        f.write(_fmt(pl))
        f.write(_fmt(t))

    subprocess.run(["./rce_units"], cwd=units_probe, check=True)
    tokens = (units_probe / "units_out.txt").read_text().split()
    got = np.array([float(x) for x in tokens[:len(p_adia)]])

    k = t0 / (PARMENTIER_A - PARMENTIER_B * t0) * (p_adia / p0) ** PARMENTIER_A
    ref = PARMENTIER_A * k / (1 + PARMENTIER_B * k)

    rel = np.abs(got - ref).max() / np.abs(ref).max()
    assert rel < 1e-14, f"adiabat mismatch, worst relative difference {rel:.2e}"

    # The anchor has to be an exact fixed point, or the convective zone would
    # step discontinuously across the boundary every time it is rebuilt.
    at_p0 = PARMENTIER_A * (t0 / (PARMENTIER_A - PARMENTIER_B * t0)) / (
        1 + PARMENTIER_B * t0 / (PARMENTIER_A - PARMENTIER_B * t0))
    assert at_p0 == pytest.approx(t0, rel=1e-14)


def test_tabulated_adiabat_matches_python(adiabat_probe):
    """``rce_grad_ad`` / ``rce_integrate_adiabat`` vs ``carmapy.adiabat``.

    Both sides read the same shipped table, so the gradient must agree to
    round-off. The integrated profile is compared more loosely because the two
    use different schemes -- RK4 with a capped step in Fortran, adaptive RK45
    in scipy -- and agreement there is the statement that both solve the same
    ODE, not that they take the same steps.
    """
    from carmapy.adiabat import TABLE_FILE, grad_ad, tabulated_adiabat

    # Pa in the Fortran, barye in Python: a factor of 10.
    t0, p0_pa = 1800.0, 1e6
    pres_pa = np.geomspace(1.2e6, 2.5e7, 20)   # stays inside the table

    with open(adiabat_probe / "adiabat_case.txt", "w") as f:
        f.write(f"{TABLE_FILE}\n")
        f.write(f"{len(pres_pa)} {t0!r} {p0_pa!r}\n")
        f.write(_fmt(pres_pa))

    subprocess.run(["./rce_adiabat_probe"], cwd=adiabat_probe, check=True)
    tokens = (adiabat_probe / "adiabat_out.txt").read_text().split()

    n = len(pres_pa)
    got_grad = np.array([float(x) for x in tokens[:n]])
    got_temp = np.array([float(x) for x in tokens[n:2 * n]])

    ref_grad = grad_ad(t0, pres_pa * 10.0)
    np.testing.assert_allclose(got_grad, ref_grad, rtol=1e-12)

    # The two agree to about 1e-6, which is the residual integration error of
    # the schemes rather than a difference in what they integrate: the
    # gradient is only C0 (bilinear, so it has kinks at cell edges), which
    # caps how far either scheme's order of accuracy carries.
    ref_temp = tabulated_adiabat(pres_pa * 10.0, t0, p0_pa * 10.0)
    np.testing.assert_allclose(got_temp, ref_temp, rtol=1e-5)


# --------------------------------------------------------------------------
# eddy diffusion
# --------------------------------------------------------------------------

def _kzz_case(nz=40, t_int=1200.0, wt_mol=2.2, grav=1e5, f_frac=None):
    """A column for the Kzz probe: pressures [Pa], temperatures, net flux.

    The profile is a Parmentier adiabat, and the net flux ramps from nothing at
    the base to the whole internal flux at the top, so the column has a
    convective interior under a radiative upper atmosphere -- both branches of
    the mixing length floor get exercised.
    """
    from carmapy.adiabat import parmentier_adiabat

    SIGMA = 0.56687e-4

    # Bottom to top, and stopping at 0.01 bar so the tabulated adiabat is
    # evaluated inside its own grid rather than against its clamped edge.
    pl_pa = np.geomspace(3e7, 1e3, nz + 1)[::-1]
    p_pa = np.sqrt(pl_pa[1:] * pl_pa[:-1])
    t = parmentier_adiabat(p_pa * 10.0, 2600.0, 3e8)

    if f_frac is None:
        # W/m^2: fnet is SI on both sides of the port.
        f_frac = np.linspace(0.0, 1.0, nz + 1) ** 2
    fnet = f_frac * SIGMA * t_int**4 * 1e-3

    return dict(nz=nz, t_int=t_int, wt_mol=wt_mol, grav=grav,
                pl_pa=pl_pa, p_pa=p_pa, t=t, fnet=fnet)


def _run_kzz_probe(probe, case, adiabat, mixl_scale=1.0):
    """Run the Fortran probe on ``case``, returning ``(kzz, rescale)``."""
    from carmapy.adiabat import ADIABAT_CODES, TABLE_FILE

    with open(probe / "kzz_case.txt", "w") as f:
        f.write(f"{TABLE_FILE}\n")
        f.write(f"{case['nz']} {ADIABAT_CODES[adiabat]} "
                f"{float(case['t_int'])!r} {float(case['wt_mol'])!r} "
                f"{float(case['grav'])!r} {float(mixl_scale)!r}\n")
        f.write(_fmt(case["p_pa"]))
        f.write(_fmt(case["pl_pa"]))
        f.write(_fmt(case["t"]))
        f.write(_fmt(case["fnet"]))

    subprocess.run(["./rce_kzz_probe"], cwd=probe, check=True)
    tokens = (probe / "kzz_out.txt").read_text().split()

    got = np.array([float(x) for x in tokens[:case["nz"] + 1]])
    return got, float(tokens[case["nz"] + 1])


def _kzz_reference(case, adiabat, mixl_scale=1.0):
    """``carmapy.kzz`` on the same column. Pa -> barye at the boundary."""
    from carmapy.kzz import mixing_length_kzz

    return mixing_length_kzz(case["p_pa"] * 10.0, case["t"],
                             case["pl_pa"] * 10.0, case["fnet"],
                             t_int=case["t_int"], wt_mol=case["wt_mol"],
                             surface_grav=case["grav"], adiabat=adiabat,
                             mixl_scale=mixl_scale, return_rescale=True)


@pytest.mark.parametrize("adiabat", ["parmentier", "table"])
@pytest.mark.parametrize("mixl_scale", [1.0, 0.4, 2.5])
def test_kzz_matches_python(kzz_probe, adiabat, mixl_scale):
    """``rce_kzz`` against ``carmapy.kzz.mixing_length_kzz``.

    A pure port check: both sides evaluate the same closed-form expression on
    the same column, so the only differences allowed are the ones arithmetic
    ordering makes. Run against both adiabats because the lapse ratio is where
    the two models enter, and a gradient wired to the wrong one would only
    misplace the mixing length by a factor of order unity -- big enough to
    matter, small enough to look plausible.

    The mixing length scale is swept alongside, above and below 1, because the
    two sides have to agree on *where* it is applied and not merely that it
    exists: multiplying before the 0.1*H floor rather than after it gives the
    same answer wherever the floor is slack and a different one everywhere it
    binds, which is exactly the stable upper atmosphere.
    """
    case = _kzz_case()

    got, got_rescale = _run_kzz_probe(kzz_probe, case, adiabat, mixl_scale)
    ref, ref_rescale = _kzz_reference(case, adiabat, mixl_scale)

    np.testing.assert_allclose(got, ref, rtol=1e-11)
    assert got_rescale == pytest.approx(ref_rescale, rel=1e-12)


def test_kzz_mixl_scale_is_a_four_thirds_power(kzz_probe):
    """Scaling the mixing length by ``f`` scales Kzz by ``f**(4/3)``.

    Kzz goes as ``H*(l/H)**(4/3)`` and nothing else in the expression touches
    ``l``, so the response is exact rather than approximate. Worth pinning
    because the 4/3 is the whole reason the knob needs documenting: a user
    reaching for "twice the mixing" wants 1.68, not 2.
    """
    case = _kzz_case()

    base, _ = _run_kzz_probe(kzz_probe, case, "parmentier", 1.0)
    for f in (0.25, 0.5, 2.0, 3.0):
        scaled, _ = _run_kzz_probe(kzz_probe, case, "parmentier", f)
        np.testing.assert_allclose(scaled / base, f ** (4 / 3), rtol=1e-9)


def test_kzz_mixl_scale_must_be_positive():
    """A non-positive scale would put zero or a negative into a 4/3 power."""
    from carmapy.kzz import mixing_length_kzz

    case = _kzz_case()
    for bad in (0.0, -1.0):
        with pytest.raises(ValueError, match="must be positive"):
            mixing_length_kzz(case["p_pa"] * 10.0, case["t"],
                              case["pl_pa"] * 10.0, case["fnet"],
                              t_int=case["t_int"], wt_mol=case["wt_mol"],
                              surface_grav=case["grav"], mixl_scale=bad)


def test_kzz_floors_a_radiative_column(kzz_probe):
    """A column carrying its whole flux radiatively sits on the flux floor.

    With ``fnet`` already equal to the internal flux everywhere there is no
    convective flux left to drive mixing, so ``chf`` collapses onto
    ``sigma*(0.05*Teff)^4``. What is being checked is that the result is a
    small positive number rather than a NaN: the cube root in the mixing
    length velocity is what makes the floor load-bearing rather than cosmetic.
    """
    case = _kzz_case(f_frac=np.ones(41))

    got, _ = _run_kzz_probe(kzz_probe, case, "parmentier")
    ref, _ = _kzz_reference(case, "parmentier")

    assert np.all(np.isfinite(got))
    assert np.all(got > 0.0)
    np.testing.assert_allclose(got, ref, rtol=1e-11)

    # Floored at 0.05^4 of the flux, and the velocity goes as its cube root, so
    # the column sits well below the convective case.
    #
    # The bottom level is excluded, and equal between the two by construction:
    # the deepest layer is *defined* to carry the whole flux convectively, so
    # no amount of radiative transport above it changes what it mixes with.
    convective, _ = _run_kzz_probe(kzz_probe, _kzz_case(), "parmentier")
    assert got[0] == pytest.approx(convective[0])
    assert np.all(got[1:-5] < convective[1:-5])


def test_kzz_rises_with_internal_flux(kzz_probe):
    """Kzz goes as the cube root of the convective flux.

    Doubling the internal flux at fixed profile multiplies ``chf`` by four and
    Kzz by ``4**(1/3)`` -- everything else in the expression depends only on
    the temperature profile, which is held. A sharper statement than a
    monotonicity check, and it pins the exponent.
    """
    lo = _kzz_case(t_int=1200.0)
    hi = _kzz_case(t_int=1200.0 * np.sqrt(2.0))   # sigma*T^4 doubles twice over

    got_lo, _ = _run_kzz_probe(kzz_probe, lo, "parmentier")
    got_hi, _ = _run_kzz_probe(kzz_probe, hi, "parmentier")

    np.testing.assert_allclose(got_hi / got_lo, 4.0 ** (1 / 3), rtol=1e-6)


def test_level_temps_and_regrid(units_probe):
    """Level reconstruction and the hydrostatic regrid vs numpy references.

    The regrid reference is ``Carma.calculate_z``'s ``I_CART`` branch,
    transcribed here so the test does not need the rest of the Python package
    (and so a change to either side shows up as a failure rather than as two
    matching bugs).
    """
    # rce_regrid_z takes gravity in m/s^2; only rce_init does the cgs
    # conversion. 316 m/s^2 is log g = 4.5 cgs, a typical brown dwarf.
    nz_g, igridv, wtmol, grav = 6, I_CART, 2.3, 316.0
    p = np.geomspace(1e7, 1e2, nz_g)
    pl = np.geomspace(2e7, 5e1, nz_g + 1)
    t = np.linspace(1900.0, 700.0, nz_g)

    tokens = (units_probe / "units_out.txt").read_text().split()
    off = 25                                   # past the adiabat block
    tl = np.array([float(x) for x in tokens[off:off + nz_g + 1]])
    off += nz_g + 1
    zl = np.array([float(x) for x in tokens[off:off + nz_g + 1]])
    off += nz_g + 1
    zc = np.array([float(x) for x in tokens[off:off + nz_g]])

    # --- level temperatures: linear in ln p, endpoints extrapolated ---
    ref_tl = np.zeros(nz_g + 1)
    for i in range(1, nz_g):
        w = (np.log(pl[i]) - np.log(p[i])) / (np.log(p[i - 1]) - np.log(p[i]))
        ref_tl[i] = t[i] + w * (t[i - 1] - t[i])
    w = (np.log(pl[0]) - np.log(p[0])) / (np.log(p[0]) - np.log(p[1]))
    ref_tl[0] = t[0] + w * (t[0] - t[1])
    w = (np.log(pl[nz_g]) - np.log(p[-1])) / (np.log(p[-1]) - np.log(p[-2]))
    ref_tl[nz_g] = t[-1] + w * (t[-1] - t[-2])

    assert np.abs(tl - ref_tl).max() / np.abs(ref_tl).max() < 1e-14

    # --- hydrostatic grid: Carma.calculate_z, I_CART branch, in SI ---
    scale_h = K_B * ref_tl / (wtmol * M_PROT * grav)
    ref_zl = np.zeros(nz_g + 1)
    for i in range(1, nz_g + 1):
        ref_zl[i] = ref_zl[i - 1] + scale_h[i] * np.log(pl[i - 1] / pl[i])
    ref_zc = 0.5 * (ref_zl[1:] + ref_zl[:-1])

    assert np.abs(zl - ref_zl).max() / np.abs(ref_zl).max() < 1e-13
    assert np.abs(zc - ref_zc).max() / np.abs(ref_zc).max() < 1e-13


def test_logp_grid_is_not_regridded(units_probe):
    """``I_LOGP`` altitudes do not depend on T, so the regrid must be a no-op.

    Silently moving them would put sedimentation and eddy diffusion on a grid
    inconsistent with the one CARMA was configured with.
    """
    nz_g, wtmol, grav = 6, 2.3, 316.0
    p = np.geomspace(1e7, 1e2, nz_g)
    pl = np.geomspace(2e7, 5e1, nz_g + 1)
    t = np.linspace(1900.0, 700.0, nz_g)
    p_adia = np.array([1e5])

    with open(units_probe / "units_case.txt", "w") as f:
        f.write("1 1900.0 1e7\n")
        f.write(_fmt(p_adia))
        f.write(f"{nz_g} {I_LOGP} {wtmol!r} {grav!r}\n")
        f.write(_fmt(p))
        f.write(_fmt(pl))
        f.write(_fmt(t))

    subprocess.run(["./rce_units"], cwd=units_probe, check=True)
    tokens = (units_probe / "units_out.txt").read_text().split()

    off = 1 + (nz_g + 1)
    zl = np.array([float(x) for x in tokens[off:off + nz_g + 1]])
    off += nz_g + 1
    zc = np.array([float(x) for x in tokens[off:off + nz_g]])

    # The probe zeroes both arrays before calling, so untouched means zero.
    assert np.all(zl == 0.0)
    assert np.all(zc == 0.0)


# --------------------------------------------------------------------------
# relaxation to radiative-convective equilibrium
# --------------------------------------------------------------------------

def _adiabat(p, t0, p0):
    k = t0 / (PARMENTIER_A - PARMENTIER_B * t0) * (p / p0) ** PARMENTIER_A
    return PARMENTIER_A * k / (1 + PARMENTIER_B * k)


def _relax(probe, kappa_ln=-57.0, t_int=1000.0, nstep=16000, accel=20.0,
           dt_max=5.0, gap_max=1, nz=40, grav_cgs=31600.0,
           wtmol=2.3, cp_cgs=1.3e8, dtime=250.0, isothermal=False):
    """Relax a clear grey column and return its final state.

    ``accel`` and ``dt_max`` are chosen so the step stays inside the local
    radiative time constant; see the note in ``test_energy_conservation``.
    """
    nband = write_grey_ck(probe / "grey.ck", kappa_ln)

    pl = np.geomspace(1e7, 1.0, nz + 1)
    p = np.sqrt(pl[:-1] * pl[1:])
    t = np.full(nz, 1000.0) if isothermal else _adiabat(p, 2000.0, 1e7)

    scale_h = K_B * t / (wtmol * M_PROT * grav_cgs / 100.0)
    zl = np.zeros(nz + 1)
    for i in range(1, nz + 1):
        zl[i] = zl[i - 1] + scale_h[i - 1] * np.log(pl[i - 1] / pl[i])
    zc = 0.5 * (zl[1:] + zl[:-1])

    with open(probe / "run_case.txt", "w") as f:
        f.write("grey.ck\n")
        f.write(f"{nz} 1 1 {nband} 1 {I_CART} 0 {gap_max} {nstep}\n")
        f.write(f"{t_int!r} {cp_cgs!r} {grav_cgs!r} {wtmol!r}\n")
        f.write(f"{accel!r} {dt_max!r} 1.0 0.02 {dtime!r}\n")
        for arr in (p, pl, t, zl, zc):
            f.write(_fmt(arr))

    subprocess.run(["./rce_run"], cwd=probe, check=True)
    tok = (probe / "run_out.txt").read_text().split()

    o = 0
    t_out = np.array([float(x) for x in tok[o:o + nz]]); o += nz
    fnet = np.array([float(x) for x in tok[o:o + nz + 1]]); o += nz + 1
    dtdt = np.array([float(x) for x in tok[o:o + nz]]); o += nz
    tau_rad = np.array([float(x) for x in tok[o:o + nz]]); o += nz
    is_conv = np.array([float(x) for x in tok[o:o + nz]]) > 0.5; o += nz
    resid = float(tok[o]); o += 1
    nsolve, nclamp, nz_rcb, nzone = (int(x) for x in tok[o:o + 4])

    return dict(t=t_out, fnet=fnet, dtdt=dtdt, tau_rad=tau_rad,
                is_conv=is_conv, resid=resid, nsolve=nsolve, nclamp=nclamp,
                nz_rcb=nz_rcb, nzone=nzone, p=p, pl=pl, nstep=nstep)


def _pfact(t_bar, p_lo, p_hi):
    """The adiabatic temperature ratio ``rce_pfact`` computes, in numpy."""
    return (p_lo / p_hi) ** (PARMENTIER_A - PARMENTIER_B * t_bar)


def test_energy_conservation(run_probe):
    """A clear grey column must relax so that F_TOA -> sigma*Teff^4.

    The internal flux is imposed as the net flux through the base of the
    column, and nothing carries it upward except the radiative solve and the
    convective adjustment, so reaching ``sigma*Teff^4`` at the top is an
    end-to-end check on both.

    ``accel`` is set so the effective step stays inside the local radiative
    time constant. Pushing it far beyond that leaves the update clamp-limited:
    it stays bounded, but it orbits equilibrium instead of reaching it, which
    is exactly what a rising ``nclamp`` is there to reveal.
    """
    r = _relax(run_probe)
    target = SIGMA_SB * 1000.0 ** 4

    ratio = r["fnet"][-1] / target
    assert abs(ratio - 1.0) < 1e-4, \
        f"F_TOA is {ratio:.8f} of sigma*Teff^4"

    # In radiative equilibrium the net flux is height-independent above the
    # convective boundary -- every layer passes on exactly what it receives.
    radiative = r["fnet"][r["nz_rcb"]:]
    spread = radiative.max() - radiative.min()
    assert spread / target < 1e-4, \
        f"net flux varies by {spread:.3f} W/m^2 across the radiative zone"

    assert r["nz_rcb"] > 0, "test case should have a convective interior"
    assert r["nzone"] >= 1


def test_convective_zone_is_neutral_and_the_rest_is_not(run_probe):
    """The adjustment must leave exactly the layers it marks on the adiabat.

    The mask is the run's only statement about where convection is, so it has
    to agree with the profile it hands back: neutral inside a zone, strictly
    subadiabatic outside one.
    """
    r = _relax(run_probe)
    p, t, mask = r["p"], r["t"], r["is_conv"]

    dp = r["pl"][:-1] - r["pl"][1:]
    t_bar = (dp[:-1] * t[:-1] + dp[1:] * t[1:]) / (dp[:-1] + dp[1:])
    ratio = t[:-1] / (_pfact(t_bar, p[:-1], p[1:]) * t[1:])

    # An interface inside a zone has both its layers marked.
    inside = mask[:-1] & mask[1:]
    assert inside.any(), "no convective interface to check"
    assert np.abs(ratio[inside] - 1.0).max() < 1e-4, \
        f"a marked interface is off the adiabat by {np.abs(ratio[inside]-1).max():.2e}"

    # Anything not marked must be genuinely stable, not merely unconverged.
    outside = ~(mask[:-1] | mask[1:])
    if outside.any():
        assert ratio[outside].max() < 1.0

    assert r["resid"] < 1e-3, \
        f"adjustment left {r['resid']:.2e} fractional superadiabaticity"


def test_isothermal_start_finds_its_own_interior(run_probe):
    """Relaxation from an isothermal start must find the interior itself.

    Starting far from equilibrium with no boundary supplied is the case the
    fixed-``p_rcb`` scheme could not handle: the convective zone has to be
    located, not assumed.
    """
    r = _relax(run_probe, isothermal=True, nstep=20000)
    target = SIGMA_SB * 1000.0 ** 4

    assert abs(r["fnet"][-1] / target - 1.0) < 1e-3, \
        f"F_TOA is {r['fnet'][-1]/target:.6f} of sigma*Teff^4"
    assert r["nz_rcb"] > 0, "an isothermal start must grow a convective interior"


def test_adaptive_cadence_matches_every_step(run_probe):
    """The adaptive solve cadence must not change the answer it converges to.

    Holding ``dTdt`` between solves is the whole reason a full radiative
    transfer solve is affordable, so the triggers have to be tight enough that
    skipping solves costs accuracy no larger than the ``rad_dT_tol`` budget
    they are set from.
    """
    every = _relax(run_probe, gap_max=1)
    adaptive = _relax(run_probe, gap_max=50)

    assert adaptive["nsolve"] < every["nsolve"], \
        "adaptive cadence did not skip any solves"

    # dt_tol is 1.0 K in _relax; allow a small multiple of it for drift that
    # accumulates across the run.
    dt = np.abs(adaptive["t"] - every["t"]).max()
    assert dt < 3.0, f"adaptive cadence moved the profile by {dt:.3f} K"

    target = SIGMA_SB * 1000.0 ** 4
    assert abs(adaptive["fnet"][-1] / target - 1.0) < 1e-4


# --------------------------------------------------------------------------
# the heating Jacobian
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def jac_probe(tmp_path_factory):
    import shutil
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")
    return _build(tmp_path_factory.mktemp("rce_jac"), "rce_jac", PROBE_JAC)


def _jacobian(probe, dt_in, dlnk_dlnt=2.0, kappa_ln=-57.0, nz=40,
              grav_cgs=31600.0, wtmol=2.3, cp_cgs=1.3e8):
    """Measure the Jacobian on a grey column, plus J*dT against the real change.

    ``dlnk_dlnt`` defaults non-zero here: with the opacity independent of
    temperature every off-diagonal term is identically zero and the test would
    be vacuous.
    """
    nband = write_grey_ck(probe / "grey.ck", kappa_ln, dlnk_dlnt=dlnk_dlnt)

    pl = np.geomspace(1e7, 1.0, nz + 1)
    p = np.sqrt(pl[:-1] * pl[1:])
    t = _adiabat(p, 2000.0, 1e7)

    scale_h = K_B * t / (wtmol * M_PROT * grav_cgs / 100.0)
    zl = np.zeros(nz + 1)
    for i in range(1, nz + 1):
        zl[i] = zl[i - 1] + scale_h[i - 1] * np.log(pl[i - 1] / pl[i])

    with open(probe / "jac_case.txt", "w") as f:
        f.write("grey.ck\n")
        f.write(f"{nz} 1 1 {nband} {I_CART}\n")
        f.write(f"1000.0 {cp_cgs!r} {grav_cgs!r} {wtmol!r}\n")
        for arr in (p, pl, t, zl, dt_in):
            f.write(_fmt(arr))

    subprocess.run(["./rce_jac"], cwd=probe, check=True)

    tok = (probe / "jac_out.txt").read_text().split()
    o = 0
    jac = np.array([float(x) for x in tok[o:o + nz * nz]]).reshape(nz, nz, order="F")
    o += nz * nz
    pred = np.array([float(x) for x in tok[o:o + nz]]); o += nz
    actual = np.array([float(x) for x in tok[o:o + nz]]); o += nz
    return jac, pred, actual


def test_jacobian_predicts_a_real_step(jac_probe):
    """``J*dT`` must match what the solver actually does for that ``dT``.

    This is the check that the Jacobian is *right*, as opposed to merely
    self-consistent: re-running the same finite differences would agree with
    itself no matter what it measured. Here the prediction is compared against
    an independent evaluation of the real flux solver on a perturbed profile,
    so a wrong sign, a transposed index or a missed temperature dependence all
    show up.

    The residual is checked for first-order convergence rather than against a
    fixed tolerance. What is left over after ``J*dT`` is the linearisation
    remainder, second order in the step, so measured against a first-order
    signal it must fall in proportion to the step. A Jacobian that is simply
    wrong leaves a residual that does not shrink at all.

    The steps are large on purpose. Measured on this column, the relative
    residual falls cleanly from 0.38 at 18 K to 0.038 at 2.3 K -- halving with
    each halving of the step -- and then stops falling below about 0.5 K, where
    the flux difference reaches the solver's own reproducibility floor and the
    signal is no longer above the noise. Testing in that regime would measure
    round-off. It is also why ``JAC_DT`` is 2 K rather than something smaller:
    the perturbation used to build the Jacobian has to sit above the same floor.
    """
    rng = np.random.default_rng(4)
    nz = 40
    step = rng.normal(scale=0.5, size=nz)

    ratios = []
    for scale_factor in (8.0, 4.0, 2.0):
        jac, pred, actual = _jacobian(jac_probe, step * scale_factor)
        assert np.all(np.isfinite(jac))
        ratios.append(np.abs(pred - actual).max() / np.abs(actual).max())

    for coarse, fine in zip(ratios, ratios[1:]):
        assert fine < 0.6 * coarse, \
            f"residual {ratios} does not fall with the step; J is not the derivative"


def test_jacobian_is_not_tridiagonal(jac_probe):
    """A tridiagonal truncation misses most of the operator, which is why the
    implicit step solves a dense system.

    A standing measurement rather than an assumption. The bandwidth is a
    property of the vertical discretisation -- of the layer-to-level
    reconstruction the Planck source is built from -- so a change there could
    make a banded solve legitimate and much cheaper. This is how you would find
    out.

    **This column understates the case.** A grey table gives every band the
    same optical depth, so a layer either exchanges with its neighbours or
    radiates to space, and the upper column here comes out almost perfectly
    tridiagonal. Long-range coupling needs spectral diversity -- opaque bands
    reaching distant layers while windows do not -- which only a real table
    has. On the 196-band Sonora column ``captured(1)`` is 0.47-0.58 at *every*
    depth and ``captured(16)`` still only 0.94 in the optically thin middle.
    Those numbers are what motivated the dense solve; what is asserted below is
    the part this table can actually show.
    """
    nz = 40
    jac, _, _ = _jacobian(jac_probe, np.zeros(nz))

    def captured(k, w):
        col = np.abs(jac[:, k])
        lo, hi = max(k - w, 0), min(k + w + 1, nz)
        return col[lo:hi].sum() / col.sum()

    # The optically thick lower column, where even grey opacity couples a layer
    # to more than its immediate neighbours.
    for k in (2, 5, 10, 15):
        c1 = captured(k, 1)
        assert c1 < 0.7, \
            f"layer {k}: a tridiagonal captures {c1:.2f}; a banded solve may now be viable"

    # Cooling to space has to dominate a layer's response to its own
    # temperature: warm it and it loses more than it gains back.
    for k in (2, 5, 10, 15, 25, 35):
        assert jac[k, k] < 0.0, f"layer {k} has a positive diagonal at depth"

    # Not asserted, but the reason a banded fit cannot be rescued: on the real
    # 196-band column the near-neighbour terms come out with the *same* sign as
    # the diagonal (measured around layer 51:
    # +1.1e-4 +2.3e-4 -2.4e-4 [-9.2e-4] -2.5e-4 +1.9e-4 +8.9e-5). That is the
    # shared-level reconstruction showing through -- perturbing a layer raises
    # the levels it shares with its neighbours, so their Planck source rises
    # and they cool, where genuine exchange would warm them. The pattern is not
    # asserted here because it is not robust across depth on a grey column.
