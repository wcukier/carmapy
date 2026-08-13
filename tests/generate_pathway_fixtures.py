"""Generate pathway autoregression reference files.

Run this script once from the repo root after building carmapy, then commit
the resulting .npz files:

    python tests/generate_pathway_fixtures.py

Re-run only when intentionally changing physics (e.g., updating vapor pressure
coefficients, switching coagulation kernels, modifying the nucleation scheme).
The generated files are committed to tests/data/pathways/.
"""
import sys
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

PATHWAY_DATA = REPO_ROOT / "tests" / "data" / "pathways"
PATHWAYS = ["tio2_hom", "tio2_coag", "kcl_hom", "mg2sio4_het", "multispecies"]


def _build_pathway(name, path, example_levels):
    from carmapy.carmapy import Carma
    from carmapy.chemistry import populate_abundances_at_cloud_base

    P, T, kzz, mu = example_levels
    c = Carma(path)
    c.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    c.set_atmospheric_parameters_from_defaults("Pure H2")
    c.set_stepping(dt=250, output_gap=5, n_tstep=50)
    c.add_P(P)
    c.add_T(T)
    c.add_kzz(kzz)
    c.calculate_z(mu)

    if name == "tio2_hom":
        c.add_gas("TiO2")
        c.add_hom_group("TiO2", 1e-8)

    elif name == "tio2_coag":
        c.add_gas("TiO2")
        c.add_hom_group("TiO2", 1e-8)
        c.add_coag("Pure TiO2")

    elif name == "kcl_hom":
        c.add_gas("KCl")
        c.add_hom_group("KCl", 1e-8)

    elif name == "mg2sio4_het":
        c.add_gas("TiO2")
        c.add_gas("Mg2SiO4")
        c.add_hom_group("TiO2", 1e-8)
        c.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2 ** (1 / 3))

    elif name == "multispecies":
        c.add_gas("TiO2")
        c.add_gas("KCl")
        c.add_gas("Mg2SiO4")
        c.add_hom_group("TiO2", 1e-8)
        c.add_hom_group("KCl", 1e-8)
        c.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2 ** (1 / 3))
        c.add_coag("Pure TiO2")

    populate_abundances_at_cloud_base(c)
    return c


def main():
    from carmapy.example import example_levels
    from carmapy.results import Results

    levels = example_levels(t=1000)

    for name in PATHWAYS:
        run_dir = PATHWAY_DATA / name / "run"
        run_dir.mkdir(parents=True, exist_ok=True)

        print(f"Running pathway {name!r} ...", end=" ", flush=True)
        c = _build_pathway(name, str(run_dir), levels)
        c.run(suppress_output=True)
        r = Results(c, read_diag=False)

        ref_path = PATHWAY_DATA / name / f"{name}_reference.npz"
        np.savez_compressed(
            ref_path,
            numden=r.numden[:, :, :, -1],
            gas_abund=r.gas_abund[:, :, -1],
        )
        size_kb = ref_path.stat().st_size / 1024
        print(f"saved {ref_path.name} ({size_kb:.1f} KB)")

    print("\nDone. Commit tests/data/pathways/ to git.")
    print("Remember: only re-run this script when intentionally changing physics.")


if __name__ == "__main__":
    main()
