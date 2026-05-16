"""Generate the small test fixture used by parsing tests.

Run this script once from the repo root after building carmapy:

    python tests/generate_fixture.py

Commits tests/data/mini_output/ to git after running.
The fixture is a 50-timestep TiO2 simulation (~15 KB of output).
"""
import sys
import os
import numpy as np
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

OUTPUT_PATH = REPO_ROOT / "tests" / "data" / "mini_output" / "mini"


def main():
    from carmapy.carmapy import Carma
    from carmapy.example import example_levels
    from carmapy.chemistry import populate_abundances_at_cloud_base

    print(f"Generating fixture at {OUTPUT_PATH} ...")

    P, T, kzz, mu = example_levels(t=1000)
    carma = Carma(str(OUTPUT_PATH))
    carma.set_physical_params(surface_grav=31600, wt_mol=float(np.mean(mu)))
    carma.set_atmospheric_parameters_from_defaults("Pure H2")
    carma.set_stepping(dt=250, output_gap=5, n_tstep=50)
    carma.add_gas("TiO2")
    carma.add_hom_group("TiO2", 1e-8)
    carma.add_P(P)
    carma.add_T(T)
    carma.add_kzz(kzz)
    carma.calculate_z(mu)
    populate_abundances_at_cloud_base(carma)
    carma.run(suppress_output=False)

    out = OUTPUT_PATH / "mini.txt"
    size_kb = out.stat().st_size / 1024
    print(f"Done. {out} ({size_kb:.1f} KB)")
    print("Now commit tests/data/mini_output/ to git.")


if __name__ == "__main__":
    main()
