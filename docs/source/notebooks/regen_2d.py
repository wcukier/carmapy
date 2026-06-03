import carmapy
from matplotlib import pyplot as plt
import numpy as np
import warnings
import os


carma = carmapy.Carma("2d_carmapy", is_2d=True)

carma.set_stepping(dt=5000, output_gap=99999, n_tstep=100000)

P_levs, T_levs, kzz_levs, U_levs, longitudes = carmapy.example.example_2d_levels()

carma.add_P(P_levs)
carma.add_T(T_levs)
carma.add_kzz(kzz_levs)


z_wind = int(np.argmin(np.abs(P_levs - 1e3)))  # 1 mbar = 1e3 barye
velocity_avg = float(np.mean(np.abs(U_levs[z_wind, :])))

carma.set_physical_params(
    surface_grav=10**(1.3 + 2),   # log g = 1.3 (SI) → cgs
    wt_mol=2.3,
    log_metallicity=0.0,
    velocity_avg=velocity_avg,
    r_planet=1.3,
    use_jovian_radius=True,
)
carma.set_atmospheric_parameters_from_defaults("Pure H2")


carma.calculate_z()
carma.extend_atmosphere(1e9)


carma.add_hom_group("TiO2", 1e-8)
carma.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2**(1/3))
carmapy.chemistry.populate_abundances_at_cloud_base(carma)


carma.run(suppress_output=False, nthreads=4, omp_schedule="dynamic")