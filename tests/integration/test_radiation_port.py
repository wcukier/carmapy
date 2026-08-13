"""Verify the Fortran radiative-transfer port against pyHARP.

The RT modules added for radiative coupling (``carma_planck``,
``carma_ckopacity``, ``carma_rtsolve``) are transliterations of pyHARP's C++,
so they can be checked against pyHARP itself rather than against a tolerance
pulled from thin air. Agreement should be at round-off; anything larger is a
port bug.

Skipped unless both ``gfortran`` and ``pyharp`` are available -- ``pyharp`` is
an optional extra (``pip install carmapy[radiation]``).

Note on the gas-giant bottom boundary: pyHARP releases before the
``hard_surface=false`` branch existed treat that flag as a no-op and always use
the terrestrial surface. Where the installed pyHARP is such a release, these
tests compare on ``hard_surface=true``, which both implement identically. The
gas-giant branch is what a brown dwarf actually needs, and is covered instead
by the energy-conservation check (net upward flux at the base == sigma*Teff^4).
"""
import math
import shutil
import subprocess
import textwrap
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).parent.parent.parent
RT_SRC = REPO_ROOT / "src" / "CARMA" / "source" / "base"
RT_MODULES = ["carma_precision_mod", "carma_planck", "carma_ckopacity",
              "carma_rtsolve"]

@pytest.fixture(scope="module")
def torch():
    """Import torch lazily, and only for tests that actually need it.

    Importing torch at module scope would pull it in during *collection*, i.e.
    in every run of the suite including ones that skip these tests. That is not
    merely wasteful: torch and PyFastChem each load their own OpenMP runtime,
    and having both in one process aborts the interpreter partway through the
    chemistry tests.
    """
    torch = pytest.importorskip("torch", reason="pyharp extra not installed")
    pytest.importorskip("pyharp", reason="pyharp extra not installed")

    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")

    torch.set_default_dtype(torch.float64)
    return torch


# pyHARP's btop_factor is not exposed through its Python bindings, so it holds
# its default of exp(-1); the Fortran is handed the same value.
BTOP_FACTOR = math.exp(-1.0)

PROBE_SRC = textwrap.dedent("""\
    program toon_probe
      use carma_precision_mod
      use carma_rtsolve
      implicit none
      integer :: ncase, icase, nlay, nlev, k
      integer :: tef, ihard, idelta
      real(kind=f) :: a_surf, btop_factor
      real(kind=f), allocatable :: dtau(:), w0(:), g(:), be(:)
      real(kind=f), allocatable :: f_up(:), f_dn(:)

      open(unit=31, file='toon_cases.txt', status='old', action='read')
      open(unit=32, file='toon_out.txt', status='replace')
      read(31,*) ncase

      do icase = 1, ncase
        read(31,*) nlay, a_surf, tef, btop_factor, ihard, idelta
        nlev = nlay + 1
        allocate(dtau(nlay), w0(nlay), g(nlay), be(nlev), f_up(nlev), f_dn(nlev))
        read(31,*) (dtau(k), k = 1, nlay)
        read(31,*) (w0(k),   k = 1, nlay)
        read(31,*) (g(k),    k = 1, nlay)
        read(31,*) (be(k),   k = 1, nlev)

        call toon_lw_column(nlay, dtau, w0, g, be, a_surf, tef, btop_factor, &
                            ihard == 1, idelta == 1, f_up, f_dn)

        write(32,'(2i8)') icase, nlay
        write(32,'(1e26.17)') (f_up(k), k = 1, nlev)
        write(32,'(1e26.17)') (f_dn(k), k = 1, nlev)
        deallocate(dtau, w0, g, be, f_up, f_dn)
      end do

      close(31)
      close(32)
    end program toon_probe
""")


@pytest.fixture(scope="module")
def probe(tmp_path_factory, torch):
    """Compile the RT modules and a driver that exercises the Toon solver."""
    build = tmp_path_factory.mktemp("rt_port")

    for mod in RT_MODULES:
        src = RT_SRC / f"{mod}.F90"
        assert src.exists(), f"missing {src}"
        subprocess.run(
            ["gfortran", "-ffree-line-length-512", "-O2", "-c", str(src)],
            cwd=build, check=True)

    (build / "toon_probe.F90").write_text(PROBE_SRC)
    subprocess.run(
        ["gfortran", "-ffree-line-length-512", "-O2", "-o", "toon_probe",
         "toon_probe.F90"] + [f"{m}.o" for m in RT_MODULES],
        cwd=build, check=True)

    return build


def _pyharp_supports_gas_giant(torch):
    """True if the installed pyHARP implements the hard_surface=false branch."""
    from pyharp import ToonMcKay89, ToonMcKay89Options, bbflux_wavenumber

    nlay = 3
    temp = np.array([1800.0, 1400.0, 1000.0, 600.0])
    prop = torch.zeros((1, 1, nlay, 3))
    prop[0, 0, :, 0] = torch.tensor([1.0, 0.5, 0.1])

    out = []
    for flags in ("planck", "planck,hard_surface"):
        opts = ToonMcKay89Options()
        opts.flags(flags)
        opts.top_emission_flag(0)
        opts.wave_lower([500.0])
        opts.wave_upper([700.0])
        flx = ToonMcKay89(opts).forward(
            prop, "", torch.tensor(temp).reshape(1, nlay + 1),
            albedo=torch.zeros((1, 1)))
        out.append(flx[0, 0, 0, 0].item())

    # If the two surface types give the same base upward flux, the gas-giant
    # branch is not implemented in this build.
    return not math.isclose(out[0], out[1], rel_tol=1e-12)


def test_planck_matches_pyharp(probe, torch):
    """Band-integrated Planck vs pyHARP's bbflux_wavenumber."""
    from pyharp import bbflux_wavenumber

    src = textwrap.dedent("""\
        program p
          use carma_precision_mod
          use carma_planck
          implicit none
          integer :: i, j
          real(kind=f) :: wn1, wn2, t, tp(3)
          tp = (/ 500._f, 1500._f, 2500._f /)
          open(unit=21, file='planck_out.txt', status='replace')
          do i = 1, 12
            wn1 = 30.8_f * real(i, kind=f)
            wn2 = wn1 + 40._f * real(i, kind=f)
            do j = 1, 3
              t = tp(j)
              write(21,'(4e26.17)') wn1, wn2, t, bbflux_wavenumber(wn1, wn2, t)
            end do
          end do
          close(21)
        end program p
    """)
    (probe / "pp.F90").write_text(src)
    subprocess.run(["gfortran", "-ffree-line-length-512", "-O2", "-o", "pp",
                    "pp.F90", "carma_precision_mod.o", "carma_planck.o"],
                   cwd=probe, check=True)
    subprocess.run(["./pp"], cwd=probe, check=True)

    rows = np.loadtxt(probe / "planck_out.txt")
    assert len(rows) == 36

    for wn1, wn2, temp, got in rows:
        ref = float(bbflux_wavenumber(
            float(wn1), float(wn2), torch.tensor([temp]))[0])
        assert got == pytest.approx(ref, rel=1e-12), \
            f"Planck mismatch at wn=[{wn1},{wn2}] T={temp}"


def test_toon_longwave_matches_pyharp(probe, torch):
    """Toon-McKay89 level fluxes vs pyHARP's own ToonMcKay89."""
    from pyharp import ToonMcKay89, ToonMcKay89Options, bbflux_wavenumber

    hard_values = (0, 1) if _pyharp_supports_gas_giant(torch) else (1,)

    rng = np.random.default_rng(11)
    cases = []
    for nlay in (12, 40):
        p = np.logspace(7, 0, nlay + 1)                # Pa, bottom-to-top
        temp = 1800.0 * (p / p[0]) ** 0.18 + 120.0
        thick = 8.0 * (p[:-1] / p[0]) ** 0.9 + 1e-3

        for tef in (-1, 0, 1):
            for hard in hard_values:
                for delta in (0, 1):
                    for dtau, w0, g, band in (
                        (thick, np.zeros(nlay), np.zeros(nlay), (500., 700.)),
                        (thick, np.full(nlay, 0.9), np.full(nlay, 0.6),
                         (1500., 1900.)),
                        (10.0 ** rng.uniform(-5, 1.2, nlay),
                         rng.uniform(0.0, 0.999, nlay),
                         rng.uniform(-0.4, 0.9, nlay), (2000., 2400.)),
                    ):
                        be = bbflux_wavenumber(band[0], band[1],
                                               torch.tensor(temp)).numpy()
                        cases.append(dict(nlay=nlay, dtau=dtau, w0=w0, g=g,
                                          be=be, temp=temp, band=band,
                                          tef=tef, hard=hard, delta=delta))

    with open(probe / "toon_cases.txt", "w") as f:
        f.write(f"{len(cases)}\n")
        for c in cases:
            f.write(f"{c['nlay']} 0.0 {c['tef']} {BTOP_FACTOR!r} "
                    f"{c['hard']} {c['delta']}\n")
            for key in ("dtau", "w0", "g", "be"):
                f.write(" ".join(repr(float(v))
                                 for v in np.ravel(c[key])) + "\n")

    subprocess.run(["./toon_probe"], cwd=probe, check=True)

    tokens = (probe / "toon_out.txt").read_text().split()
    i = 0
    got = []
    while i < len(tokens):
        nlay = int(tokens[i + 1])
        i += 2
        up = np.array([float(x) for x in tokens[i:i + nlay + 1]])
        i += nlay + 1
        dn = np.array([float(x) for x in tokens[i:i + nlay + 1]])
        i += nlay + 1
        got.append((up, dn))

    assert len(got) == len(cases)

    for c, (f_up, f_dn) in zip(cases, got):
        opts = ToonMcKay89Options()
        flags = "planck"
        if c["hard"]:
            flags += ",hard_surface"
        if c["delta"]:
            flags += ",delta_eddington_lw"
        opts.flags(flags)
        opts.top_emission_flag(c["tef"])
        opts.wave_lower([c["band"][0]])
        opts.wave_upper([c["band"][1]])

        prop = torch.zeros((1, 1, c["nlay"], 3))
        prop[0, 0, :, 0] = torch.tensor(c["dtau"])
        prop[0, 0, :, 1] = torch.tensor(c["w0"])
        prop[0, 0, :, 2] = torch.tensor(c["g"])

        ref = ToonMcKay89(opts).forward(
            prop, "", torch.tensor(c["temp"]).reshape(1, c["nlay"] + 1),
            albedo=torch.zeros((1, 1)))[0, 0].numpy()

        where = f"nlay={c['nlay']} tef={c['tef']} hard={c['hard']} delta={c['delta']}"
        for name, mine, theirs in (("up", f_up, ref[:, 0]),
                                   ("dn", f_dn, ref[:, 1])):
            scale = max(np.abs(theirs).max(), 1e-300)
            rel = np.abs(mine - theirs).max() / scale
            assert rel < 1e-11, f"{name} flux mismatch ({rel:.2e}) at {where}"
