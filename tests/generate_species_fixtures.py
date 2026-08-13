"""Generate species-coverage regression reference files.

Run once from the repo root after building carmapy, then commit the .npz files:

    python tests/generate_species_fixtures.py

After generation, measure run-to-run variability to guide tolerance selection:

    python tests/generate_species_fixtures.py --measure-error

The --measure-error flag runs every species twice and prints the observed
max/mean relative error for numden and gas_abund. Set NUMDEN_RTOL and
GAS_RTOL in tests/pathway/test_species_coverage.py to 2-3× the max.
"""
import argparse
import sys
import tempfile
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
SPECIES_DATA = REPO_ROOT / "tests" / "data" / "species"
sys.path.insert(0, str(REPO_ROOT / "src"))

WARMUP_NSTEP = 2000
WARMUP_GAP   = 200
FINE_NSTEP   = 500
FINE_GAP     = 5
DT           = 250

# H2O is excluded: it has lat_heat_e=NaN in gas_dict (condensation not
# supported in this configuration) and hangs at T=1000K (never condenses).
ALL_SPECIES = ["KCl", "ZnS", "Na2S", "MnS", "Cr",
               "Mg2SiO4", "Fe", "TiO2", "Al2O3", "SiO"]


def _build(species, path, levels):
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base

    P, T, kzz, mu = levels
    c = Carma(path)
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)
    c.add_gas(species)
    c.add_hom_group(species, 1e-8)
    populate_abundances_at_cloud_base(c)
    return c


def _run_warmup_and_fine(carma):
    carma.set_stepping(dt=DT, output_gap=WARMUP_GAP, n_tstep=WARMUP_NSTEP)
    carma.run(suppress_output=True)
    carma.restart = 1
    carma.set_stepping(dt=DT, output_gap=FINE_GAP, n_tstep=FINE_NSTEP)
    carma.run(suppress_output=True)


def _get_timeavg(carma):
    from carmapy.results import Results
    r = Results(carma, read_diag=False)
    return r.numden.mean(axis=-1), r.gas_abund.mean(axis=-1), r.numden.shape[-1]


def _rel_err(a, b, atol=1.0):
    """Max relative error between arrays a and b, gating near-zero cells."""
    mask = np.abs(b) > atol
    if not mask.any():
        return 0.0, 0.0
    rel = np.abs(a[mask] - b[mask]) / np.abs(b[mask])
    return float(rel.max()), float(rel.mean())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--measure-error", action="store_true",
        help="Run each species twice and report run-to-run variability.",
    )
    args = parser.parse_args()

    from carmapy.example import example_levels
    levels = example_levels(t=1000)

    if args.measure_error:
        print("Measuring run-to-run variability (2 runs per species)...\n")
        print(f"{'Species':<12} {'numden max-rel':>16} {'numden mean-rel':>16} "
              f"{'gas max-rel':>14} {'gas mean-rel':>14}")
        print("-" * 76)
        for species in ALL_SPECIES:
            with tempfile.TemporaryDirectory() as t1, \
                 tempfile.TemporaryDirectory() as t2:
                c1 = _build(species, t1 + "/" + species, levels)
                _run_warmup_and_fine(c1)
                nd1, ga1, _ = _get_timeavg(c1)

                c2 = _build(species, t2 + "/" + species, levels)
                _run_warmup_and_fine(c2)
                nd2, ga2, _ = _get_timeavg(c2)

            nd_max, nd_mean = _rel_err(nd1, nd2, atol=1.0)
            ga_max, ga_mean = _rel_err(ga1, ga2, atol=1e-8)
            print(f"{species:<12} {nd_max:>16.2e} {nd_mean:>16.2e} "
                  f"{ga_max:>14.2e} {ga_mean:>14.2e}")
        print("\nSet NUMDEN_RTOL = 2-3× max numden rel, GAS_RTOL = 2-3× max gas rel.")
        return

    SPECIES_DATA.mkdir(parents=True, exist_ok=True)
    for species in ALL_SPECIES:
        species_dir = SPECIES_DATA / species
        species_dir.mkdir(parents=True, exist_ok=True)
        ref_path = species_dir / f"{species}_reference.npz"

        with tempfile.TemporaryDirectory() as tmp:
            print(f"  {species}: running warm-up + fine ...", end=" ", flush=True)
            c = _build(species, tmp + "/" + species, levels)
            _run_warmup_and_fine(c)
            numden_avg, gas_avg, n_avg = _get_timeavg(c)

        np.savez_compressed(ref_path, numden=numden_avg, gas_abund=gas_avg, n_avg=n_avg)
        size_kb = ref_path.stat().st_size / 1024
        print(f"saved ({size_kb:.1f} KB, n_avg={n_avg})")

    print("\nDone. Commit tests/data/species/ to git.")
    print("Then run with --measure-error to set tight tolerances.")


if __name__ == "__main__":
    main()
