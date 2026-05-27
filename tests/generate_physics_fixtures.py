"""Generate physics-switch regression reference files.

Run once from the repo root after building carmapy, then commit the .npz files:

    python tests/generate_physics_fixtures.py

After generation, measure run-to-run variability to guide tolerance selection:

    python tests/generate_physics_fixtures.py --measure-error

The --measure-error flag runs every variant twice and prints the observed
max/mean relative error for numden and gas_abund. Set NUMDEN_RTOL and
GAS_RTOL in tests/pathway/test_physics_switches.py to 2-3× the max.
"""
import argparse
import math
import sys
import tempfile
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
PHYSICS_DATA = REPO_ROOT / "tests" / "data" / "physics"
sys.path.insert(0, str(REPO_ROOT / "src"))

WARMUP_NSTEP = 2000
WARMUP_GAP   = 200
FINE_NSTEP   = 500
FINE_GAP     = 5
DT           = 250

ALL_VARIANTS = [
    "coag_on", "coag_off",
    "het_on", "het_off",
    "nbin_20", "nbin_40",
    "high_lat_heat", "high_viscosity",
    "wind_up", "wind_down",
]

_RMU_1 = 1.7970e-6
_RMU_2 = 0.685
_RMU_3 = -0.59
_RMU_4 = 140
_THCOND_0 = 7992.77
_THCOND_1 = 38.08
_THCOND_2 = -1.2585e-4
_CP = 1.3e8


def _build(variant, path, levels):
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base

    P, T, kzz, mu = levels
    c = Carma(path)
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)

    if variant in ("nbin_20", "nbin_40"):
        c.NBIN = int(variant.split("_")[1])

    c.calculate_z(mu)

    if variant in ("coag_on", "coag_off",
                   "nbin_20", "nbin_40",
                   "high_lat_heat", "high_viscosity",
                   "wind_up", "wind_down"):
        c.add_gas("TiO2")
        if variant == "high_lat_heat":
            gas = c.gases["TiO2"]
            natural_lhe = gas.vp_tcoeff * math.log(10) * 8.314e7 / gas.wtmol_dif
            gas.lat_heat_e = natural_lhe * 100
        c.add_hom_group("TiO2", 1e-8)
        if variant == "coag_on":
            c.add_coag("Pure TiO2")

    elif variant == "het_on":
        c.add_gas("TiO2")
        c.add_gas("Mg2SiO4")
        c.add_hom_group("TiO2", 1e-8)
        c.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2 ** (1 / 3))

    elif variant == "het_off":
        c.add_gas("Mg2SiO4")
        c.add_hom_group("Mg2SiO4", 1e-8)

    if variant == "high_viscosity":
        c.set_atmospheric_parameters(
            rmu_1=_RMU_1 * 1000,
            rmu_2=_RMU_2,
            rmu_3=_RMU_3,
            rmu_4=_RMU_4,
            thcond_0=_THCOND_0,
            thcond_1=_THCOND_1,
            thcond_2=_THCOND_2,
            cp=_CP,
        )

    if variant == "wind_up":
        c.add_vertical_winds(np.full(c.NZ, +100.0))
    elif variant == "wind_down":
        c.add_vertical_winds(np.full(c.NZ, -100.0))

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


def _rel_err(a, b, atol):
    mask = np.abs(b) > atol
    if not mask.any():
        return 0.0, 0.0
    rel = np.abs(a[mask] - b[mask]) / np.abs(b[mask])
    return float(rel.max()), float(rel.mean())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--measure-error", action="store_true",
        help="Run each variant twice and report run-to-run variability.",
    )
    args = parser.parse_args()

    from carmapy.example import example_levels
    levels = example_levels(t=1000)

    if args.measure_error:
        print("Measuring run-to-run variability (2 runs per variant)...\n")
        print(f"{'Variant':<16} {'numden max-rel':>16} {'numden mean-rel':>16} "
              f"{'gas max-rel':>14} {'gas mean-rel':>14}")
        print("-" * 80)
        for variant in ALL_VARIANTS:
            with tempfile.TemporaryDirectory() as t1, \
                 tempfile.TemporaryDirectory() as t2:
                c1 = _build(variant, t1 + "/" + variant, levels)
                _run_warmup_and_fine(c1)
                nd1, ga1, _ = _get_timeavg(c1)

                c2 = _build(variant, t2 + "/" + variant, levels)
                _run_warmup_and_fine(c2)
                nd2, ga2, _ = _get_timeavg(c2)

            nd_max, nd_mean = _rel_err(nd1, nd2, atol=1.0)
            ga_max, ga_mean = _rel_err(ga1, ga2, atol=1e-8)
            print(f"{variant:<16} {nd_max:>16.2e} {nd_mean:>16.2e} "
                  f"{ga_max:>14.2e} {ga_mean:>14.2e}")
        print("\nSet NUMDEN_RTOL = 2-3× max numden rel, GAS_RTOL = 2-3× max gas rel.")
        return

    PHYSICS_DATA.mkdir(parents=True, exist_ok=True)
    for variant in ALL_VARIANTS:
        variant_dir = PHYSICS_DATA / variant
        variant_dir.mkdir(parents=True, exist_ok=True)
        ref_path = variant_dir / f"{variant}_reference.npz"

        with tempfile.TemporaryDirectory() as tmp:
            print(f"  {variant}: running warm-up + fine ...", end=" ", flush=True)
            c = _build(variant, tmp + "/" + variant, levels)
            _run_warmup_and_fine(c)
            numden_avg, gas_avg, n_avg = _get_timeavg(c)

        np.savez_compressed(ref_path, numden=numden_avg, gas_abund=gas_avg, n_avg=n_avg)
        size_kb = ref_path.stat().st_size / 1024
        print(f"saved ({size_kb:.1f} KB, n_avg={n_avg})")

    print("\nDone. Commit tests/data/physics/ to git.")
    print("Then run with --measure-error to set tight tolerances.")


if __name__ == "__main__":
    main()
