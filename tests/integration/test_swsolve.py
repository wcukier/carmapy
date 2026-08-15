"""``carma_swsolve``'s Toon-McKay89 shortwave, against pyHARP and against analysis.

The longwave port could only ever be checked on the terrestrial bottom
boundary, because pyHARP's dispatch has a bug that makes ``hard_surface``
a no-op (see ``claude/RADIATION_PORT_STATUS.md``). The shortwave dispatch reads
``data[1..4]`` against exactly four iterator inputs and is correct, so this port
*can* be checked end to end -- and is, here.

One real limitation. pyHARP expands ``umu0`` to ``{nwave, ncol, nlyr+1, 1}``,
the same value at every level, so it only ever exercises the uniform-mu branch.
The slant-path branch, which the reference implements and this port preserves,
has no reference to compare against. It is tested structurally instead: it must
converge on the uniform branch as the per-level angles converge.
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
MODULES = ["carma_precision_mod", "carma_rtsolve", "carma_swsolve"]

PROBE = textwrap.dedent("""\
    program sw_probe
      use carma_precision_mod
      use carma_swsolve
      implicit none
      integer :: ncase, ic, nlay, nlev, k
      real(kind=f) :: f0, w_surf
      real(kind=f), allocatable :: dtau(:), w0(:), g(:), mu(:)
      real(kind=f), allocatable :: fup(:), fdn(:)

      open(unit=31, file='sw_case.txt', status='old', action='read')
      open(unit=32, file='sw_out.txt', status='replace')

      read(31,*) ncase
      do ic = 1, ncase
        read(31,*) nlay, f0, w_surf
        nlev = nlay + 1
        allocate(dtau(nlay), w0(nlay), g(nlay), mu(nlev), fup(nlev), fdn(nlev))

        read(31,*) (dtau(k), k = 1, nlay)
        read(31,*) (w0(k),   k = 1, nlay)
        read(31,*) (g(k),    k = 1, nlay)
        read(31,*) (mu(k),   k = 1, nlev)

        call toon_sw_column(nlay, f0, mu, dtau, w0, g, w_surf, fup, fdn)

        write(32,'(i8)') nlev
        write(32,'(1e26.17)') (fup(k), k = 1, nlev)
        write(32,'(1e26.17)') (fdn(k), k = 1, nlev)

        deallocate(dtau, w0, g, mu, fup, fdn)
      end do

      close(31)
      close(32)
    end program sw_probe
""")


@pytest.fixture(scope="module")
def probe(tmp_path_factory):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran not available")

    build = tmp_path_factory.mktemp("swsolve")
    for mod in MODULES:
        src = RT_SRC / f"{mod}.F90"
        assert src.exists(), f"missing {src}"
        subprocess.run(["gfortran", "-ffree-line-length-512", "-O2", "-c", str(src)],
                       cwd=build, check=True)

    (build / "sw_probe.F90").write_text(PROBE)
    subprocess.run(["gfortran", "-ffree-line-length-512", "-O2", "-o", "sw_probe",
                    "sw_probe.F90"] + [f"{m}.o" for m in MODULES],
                   cwd=build, check=True)
    return build


@pytest.fixture(scope="module")
def torch():
    """Torch behind a fixture, never at collection time: its OpenMP runtime
    clashes with PyFastChem's and aborts the interpreter mid-suite."""
    return pytest.importorskip("torch")


def _run(probe, cases):
    """cases: list of dicts with nlay, f0, w_surf, dtau, w0, g, mu.
    Returns a list of (f_up, f_dn), bottom-to-top."""
    with open(probe / "sw_case.txt", "w") as f:
        f.write(f"{len(cases)}\n")
        for c in cases:
            f.write(f"{c['nlay']} {float(c['f0'])!r} {float(c['w_surf'])!r}\n")
            for key in ("dtau", "w0", "g", "mu"):
                f.write(" ".join(repr(float(v)) for v in c[key]) + "\n")

    subprocess.run(["./sw_probe"], cwd=probe, check=True)

    tok = (probe / "sw_out.txt").read_text().split()
    out, i = [], 0
    for _ in cases:
        nlev = int(tok[i]); i += 1
        up = np.array([float(v) for v in tok[i:i + nlev]]); i += nlev
        dn = np.array([float(v) for v in tok[i:i + nlev]]); i += nlev
        out.append((up, dn))
    return out


def _cases():
    """A spread over optical depth, scattering strength, asymmetry, incidence
    angle and surface albedo -- including the thick and thin limits where the
    exponentials in the solver are worst behaved."""
    rng = np.random.default_rng(20260813)
    cases = []

    for nlay, w0v, gv, mu0, alb in [
        (3,  0.9,  0.0,  1.0,  0.0),
        (3,  0.9,  0.8,  0.5,  0.0),
        (8,  0.5,  0.3,  0.7,  0.0),
        (8,  0.99, 0.85, 0.25, 0.0),
        (16, 0.6,  0.5,  1.0,  0.3),     # reflecting surface
        (16, 0.999, 0.9, 0.5,  0.8),     # near-conservative + bright surface
        (32, 0.3,  0.1,  0.1,  0.0),     # grazing incidence
    ]:
        dtau = rng.uniform(0.01, 2.0, nlay)
        cases.append(dict(nlay=nlay, f0=1000.0, w_surf=alb,
                          dtau=dtau,
                          w0=np.full(nlay, w0v),
                          g=np.full(nlay, gv),
                          mu=np.full(nlay + 1, mu0)))

    # Very thick and very thin columns, where exp() is capped or ~1.
    cases.append(dict(nlay=6, f0=500.0, w_surf=0.0,
                      dtau=np.full(6, 40.0), w0=np.full(6, 0.7),
                      g=np.full(6, 0.6), mu=np.full(7, 0.6)))
    cases.append(dict(nlay=6, f0=500.0, w_surf=0.0,
                      dtau=np.full(6, 1e-6), w0=np.full(6, 0.7),
                      g=np.full(6, 0.6), mu=np.full(7, 0.6)))

    # Vertically structured optical properties, the realistic case.
    nlay = 20
    cases.append(dict(nlay=nlay, f0=1361.0, w_surf=0.0,
                      dtau=np.logspace(-3, 0.5, nlay),
                      w0=np.linspace(0.2, 0.95, nlay),
                      g=np.linspace(0.0, 0.85, nlay),
                      mu=np.full(nlay + 1, 0.5)))
    return cases


def test_shortwave_matches_pyharp(probe, torch):
    """Level fluxes against pyHARP's own ToonMcKay89 in its shortwave branch."""
    from pyharp import ToonMcKay89, ToonMcKay89Options

    cases = _cases()
    got = _run(probe, cases)

    for c, (f_up, f_dn) in zip(cases, got):
        opts = ToonMcKay89Options()
        opts.flags("")                       # no "planck" -> shortwave
        opts.wave_lower([1000.0])
        opts.wave_upper([2000.0])

        prop = torch.zeros((1, 1, c["nlay"], 3), dtype=torch.float64)
        prop[0, 0, :, 0] = torch.tensor(np.asarray(c["dtau"], dtype=float))
        prop[0, 0, :, 1] = torch.tensor(np.asarray(c["w0"], dtype=float))
        prop[0, 0, :, 2] = torch.tensor(np.asarray(c["g"], dtype=float))

        ref = ToonMcKay89(opts).forward(
            prop, "", None,
            umu0=torch.tensor([float(c["mu"][0])], dtype=torch.float64),
            fbeam=torch.tensor([[float(c["f0"])]], dtype=torch.float64),
            albedo=torch.tensor([[float(c["w_surf"])]], dtype=torch.float64),
        )[0, 0].numpy()

        where = (f"nlay={c['nlay']} dtau0={c['dtau'][0]:.3g} "
                 f"mu0={c['mu'][0]:.3g} w0={np.mean(c['w0']):.3g} "
                 f"alb={c['w_surf']:.3g}")

        # Normalise by the column's own flux scale, not by each stream's max.
        # In an optically thin column the upward flux is ~1e-6 of the incident
        # beam, so dividing by its own maximum turns ordinary round-off
        # (~1e-14 W/m^2) into an apparent 1e-11 relative error while the
        # physical agreement is exact. The incident beam is the scale that
        # means something.
        scale = max(np.abs(ref).max(), c["f0"] * c["mu"][0], 1e-300)

        for name, mine, theirs in (("up", f_up, ref[:, 0]),
                                   ("dn", f_dn, ref[:, 1])):
            rel = np.abs(mine - theirs).max() / scale
            assert rel < 1e-13, f"{name} flux mismatch ({rel:.2e}) at {where}"


def test_pure_absorption_is_beer_lambert(probe):
    """With no scattering the answer is analytic, and nothing reflects."""
    nlay = 12
    mu0 = 0.6
    f0 = 800.0
    dtau = np.linspace(0.05, 0.4, nlay)

    (f_up, f_dn), = _run(probe, [dict(
        nlay=nlay, f0=f0, w_surf=0.0,
        dtau=dtau, w0=np.zeros(nlay), g=np.zeros(nlay),
        mu=np.full(nlay + 1, mu0))])

    # Cumulative optical depth measured from the top, then mapped back to the
    # bottom-to-top level ordering the solver returns.
    tau_td = np.concatenate([[0.0], np.cumsum(dtau[::-1])])
    expect = f0 * mu0 * np.exp(-tau_td / mu0)[::-1]

    np.testing.assert_allclose(f_dn, expect, rtol=1e-13, atol=0)
    assert np.all(f_up == 0.0), "a non-scattering, non-reflecting column cannot radiate up"


def test_incident_flux_at_the_top_is_mu_times_f0(probe):
    """The convention that everything downstream depends on: f0 is the flux
    normal to the beam, so the flux entering the top is mu0*f0. Getting this
    wrong rescales every irradiated run and still looks plausible."""
    f0, mu0 = 1361.0, 0.5
    nlay = 4
    (f_up, f_dn), = _run(probe, [dict(
        nlay=nlay, f0=f0, w_surf=0.0,
        dtau=np.full(nlay, 0.1), w0=np.full(nlay, 0.4),
        g=np.zeros(nlay), mu=np.full(nlay + 1, mu0))])

    assert f_dn[-1] == pytest.approx(f0 * mu0, rel=1e-12)


def test_slant_path_converges_on_the_uniform_branch(probe):
    """The per-level-mu branch has no pyHARP reference, so it is pinned to the
    uniform branch it must reduce to as the level angles converge."""
    nlay, mu0 = 10, 0.6
    base = dict(nlay=nlay, f0=1000.0, w_surf=0.0,
                dtau=np.full(nlay, 0.2), w0=np.full(nlay, 0.7),
                g=np.full(nlay, 0.5))

    uniform = dict(base, mu=np.full(nlay + 1, mu0))
    # A gradient small enough that the physical answer cannot really differ,
    # but not exactly equal -- which is what selects the other branch.
    slant = dict(base, mu=np.linspace(mu0, mu0 * (1 + 1e-9), nlay + 1))

    (u_up, u_dn), (s_up, s_dn) = _run(probe, [uniform, slant])

    for name, a, b in (("up", u_up, s_up), ("dn", u_dn, s_dn)):
        scale = max(np.abs(a).max(), 1e-300)
        assert np.abs(a - b).max() / scale < 1e-7, (
            f"{name}: slant-path branch does not reduce to the uniform one")


def test_a_bright_surface_sends_flux_back_up(probe):
    """The surface term is the one place the shortwave bottom boundary does
    anything, and a gas giant sets it to zero -- so it needs its own check."""
    nlay = 6
    common = dict(nlay=nlay, f0=1000.0,
                  dtau=np.full(nlay, 0.05), w0=np.full(nlay, 1e-14),
                  g=np.zeros(nlay), mu=np.full(nlay + 1, 1.0))

    (dark_up, _), (bright_up, _) = _run(
        probe, [dict(common, w_surf=0.0), dict(common, w_surf=0.9)])

    assert np.all(dark_up == 0.0)
    assert bright_up[0] > 0.0, "a reflecting surface must produce upward flux"
