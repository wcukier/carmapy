"""Generate the brightness-temperature regression reference.

Run this script once from the repo root after building carmapy, then commit
the resulting files in `tests/data/regression/spectrum/`:

    python tests/generate_spectrum_fixture.py

Re-run only when intentionally changing physics. The regression test
(`tests/regression/test_spectrum_regression.py`) compares against these
committed artifacts.

Configuration mirrors the my_first_carma tutorial exactly (dt=100,
wt_mol=mu[0], het Mg2SiO4-on-TiO2, populate_abundances_at_cloud_base).
After a 24000-step warm-up the simulation is restarted for a fine 10000-step
window at output_gap=10 (~1000 outputs), then the picaso cloud opacity is
time-averaged over the full fine window and a thermal spectrum is computed.

Reference artifacts written:
  - sim.sha256        SHA256 of sim.txt for the fast bitwise-equality path
  - Tb_reference.npz  brightness temperature on a regridded wavenumber grid
                      (R=1000) for the tolerance path
"""
import hashlib
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
REF_DIR = REPO_ROOT / "tests" / "data" / "regression" / "spectrum"

# picaso reference data (developer-set; override with env vars if needed)
os.environ.setdefault(
    "picaso_refdata",
    "/Users/wcukier/Dropbox/Research/Projects/24-Brown Dwarfs/picaso/reference",
)
os.environ.setdefault("PYSYN_CDBS", os.environ["picaso_refdata"])

sys.path.insert(0, str(REPO_ROOT / "src"))

import numpy as np
from carmapy.carmapy import Carma
from carmapy.example import example_levels
from carmapy.chemistry import populate_abundances_at_cloud_base
from carmapy.results import Results

WARMUP_NSTEP = 24000
WARMUP_GAP = 800
FINE_NSTEP = 10000
FINE_GAP = 10
DT = 100
WAVELENGTHS = np.linspace(1e-5, 2e-3, 10000)
SPECTRUM_R = 1000


def _build_canonical(name: str) -> Carma:
    P, T, kzz, mu = example_levels()
    c = Carma(name)
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.add_gas("TiO2")
    c.add_gas("Mg2SiO4")
    c.add_hom_group("TiO2", 1e-8)
    c.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2 ** (1 / 3))
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.set_physical_params(surface_grav=31600, wt_mol=mu[0])
    c.calculate_z()
    c.calculate_z(mu)
    c.extend_atmosphere(1e10)
    populate_abundances_at_cloud_base(c)
    return c


def _run_warmup_and_fine(carma: Carma) -> None:
    carma.set_stepping(dt=DT, output_gap=WARMUP_GAP, n_tstep=WARMUP_NSTEP)
    print(f"  Warm-up {WARMUP_NSTEP} steps at gap={WARMUP_GAP} ...", flush=True)
    carma.run(suppress_output=True)
    carma.restart = 1
    carma.set_stepping(dt=DT, output_gap=FINE_GAP, n_tstep=FINE_NSTEP)
    print(f"  Fine restart {FINE_NSTEP} steps at gap={FINE_GAP} ...", flush=True)
    carma.run(suppress_output=True)


def _hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _compute_brightness_T(carma: Carma, n_avg: int):
    from picaso import justdoit as jdi
    from picaso import justplotit as jpi

    carma.results.gen_picaso_atm_file()
    carma.results.gen_picaso_cloud_file(WAVELENGTHS, n_avg=n_avg)
    case = jdi.inputs(calculation="browndwarf")
    case.phase_angle(0)
    case.gravity(gravity=carma.surface_grav, gravity_unit=jdi.u.Unit("cm/(s**2)"))
    case.atmosphere(filename=f"{carma.name}/fastchem.atm", sep=r"\s+")
    case.clouds(filename=f"{carma.name}/clouds.atm", sep=r"\s+")
    case.inputs["output_dir"] = None
    opacity = jdi.opannection(wave_range=[0.1, 15])
    df = case.spectrum(opacity, full_output=True, calculation="thermal")
    Tb = jpi.brightness_temperature(
        {"wavenumber": df["wavenumber"], "thermal": df["thermal"]},
        R=SPECTRUM_R,
        plot=False,
    )
    wno, _ = jdi.mean_regrid(df["wavenumber"], df["thermal"], R=SPECTRUM_R)
    return np.asarray(wno), np.asarray(Tb)


def main() -> None:
    import tempfile

    REF_DIR.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        sim_path = Path(tmp) / "spectrum_ref"
        carma = _build_canonical(str(sim_path))
        _run_warmup_and_fine(carma)

        sim_txt = sim_path / "spectrum_ref.txt"
        digest = _hash_file(sim_txt)
        (REF_DIR / "sim.sha256").write_text(digest + "\n")
        print(f"  sim.txt sha256 = {digest}")

        carma.results = Results(carma, read_diag=False)
        n_avg = carma.results.numden.shape[-1]

        # Time-averaged cloud + gas state for the tutorial-regression tests
        numden_avg = carma.results.numden.mean(axis=-1)
        gas_avg = carma.results.gas_abund.mean(axis=-1)
        np.savez_compressed(
            REF_DIR / "numden_timeavg.npz",
            numden=numden_avg, n_avg=n_avg,
        )
        np.savez_compressed(
            REF_DIR / "gas_abund_timeavg.npz",
            gas_abund=gas_avg, n_avg=n_avg,
        )
        print(f"  numden_timeavg.npz / gas_abund_timeavg.npz: averaged over {n_avg} outputs")

        print(f"  Computing brightness T with n_avg={n_avg} ...", flush=True)
        wno, Tb = _compute_brightness_T(carma, n_avg=n_avg)

        np.savez_compressed(REF_DIR / "Tb_reference.npz", wavenumber=wno, Tb=Tb)
        print(f"  Tb_reference.npz: {len(Tb)} wavenumbers, T range "
              f"{Tb.min():.1f}-{Tb.max():.1f} K")

    print("\nDone. Commit tests/data/regression/spectrum/ to git.")


if __name__ == "__main__":
    main()
