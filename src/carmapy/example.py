from carmapy.carmapy import Carma
from carmapy.chemistry import populate_abundances_at_cloud_base
import os
import numpy as np

SRC = os.path.dirname(os.path.dirname(__file__))

def example_levels():
    data = np.genfromtxt(os.path.join(SRC, "example_data", "example_levels"), skip_header=1)

    P_levels   = data[:, 0]
    T_levels   = data[:, 1]
    kzz_levels = data[:, 2]
    mu_levels  = data[:, 3]


    return P_levels, T_levels, kzz_levels, mu_levels


def example_carma(name):
    P_levels, T_levels, kzz_levels, mu_levels = example_levels()

    carma = Carma(name)

    carma.set_physical_params(surface_grav=31600,
                              wt_mol=np.mean(mu_levels))


    carma.set_stepping(dt=250, output_gap=10, n_tstep=10000)


    # Optional, here to preserve ordering
    carma.add_gas("TiO2")
    carma.add_gas("Mg2SiO4")

    carma.add_hom_group("TiO2", 1e-8)
    carma.add_het_group("Mg2SiO4", "TiO2", 1e-8 * 2**(1/3))


    carma.add_P(P_levels)
    carma.add_T(T_levels)
    carma.add_kzz(kzz_levels)
    
    carma.calculate_z(mu_levels)
    carma.extend_atmosphere(1e10)

    populate_abundances_at_cloud_base(carma, metalicity=1)

    return carma