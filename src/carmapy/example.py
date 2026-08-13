"""Example simulation configs and data"""

from carmapy.carmapy import Carma
from carmapy.chemistry import populate_abundances_at_cloud_base
import os
import numpy as np

SRC = os.path.dirname(os.path.dirname(__file__))

def example_levels(t=2000) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Loads example atmospheric structure "levels."  Atmospheric levels describe
    the boundaries between atmospheric layers and thus the length of "levels" 
    arrays is one more than that of the "centers" arrays.

    The files are taken from a (default 2000 K), log g = 4.5, f_sed = 4, log [Fe/H] = 0
    sonora diamondback [#morley2024]_ run. 

    References
    ----------
    .. [#morley2024] Morley, C. V., Mukherjee, S., Marley, M. S., et al. 2024 (arXiv), 
       http://arxiv.org/abs/2402.00758

    Args
    ----
    t : int
        Temperature of the sample P-T profile to use in K.  Defaults to 2000.
        Options are 1000 or 2000.

    Returns
    -------
    P_levels : np.ndarray
        Pressure levels [barye]
    T_levels : np.ndarray) 
        Temperature levels [K]
    kzz_levels : np.ndarray
        Eddy diffusion coefficients [cm²/s]
    mu_levels : np.ndarray
        mean molecular weights [dimensionless]
    """

    if t == 1000: 
        path = "example_levels"
    elif t == 2000:
        path = "example_levels2"
    else:
        raise ValueError(
            "t must be either 1000 or 2000 (or left blank to default to 2000)")

    data = np.genfromtxt(os.path.join(SRC, "example_data", path), comments='#')

    P_levels   = data[:, 0]
    T_levels   = data[:, 1]
    kzz_levels = data[:, 2]
    mu_levels  = data[:, 3]


    return P_levels, T_levels, kzz_levels, mu_levels


def example_2d_levels() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Loads example atmospheric structure "levels" for a 2d CARMApy run.

    The bundled profile is a hot Jupiter (Teq = 1800 K, log g = 3.3 cgs,
    Rp = 1.3 Rjup, solar metallicity, no drag, no TiO/VO) derived from a GCM
    grid and latitude-averaged over |lat| <= 20°. The eddy diffusion profile
    follows Moses et al. (2021) Eq 1 with a 1e11 cm²/s ceiling.

    Returns
    -------
    P_levels : np.ndarray
        Pressure levels array (NZ) [barye]
    T_levels : np.ndarray
        Temperature levels array (NZ, NLONGITUDE) [K]
    kzz_levels : np.ndarray
        Eddy diffusion coefficients (NZ) [cm²/s]
    U_levels : np.ndarray
        Zonal wind speed (NZ, NLONGITUDE) [cm/s]
    longitudes : np.ndarray
        Longitude points (NLONGITUDE) [degrees]
    """
    levels_data = np.genfromtxt(
        os.path.join(SRC, "example_data", "2d_levels.txt"), skip_header=1
    )
    P_levels = levels_data[:, 0]    # barye
    kzz_levels = levels_data[:, 1]

    T_levels = np.genfromtxt(os.path.join(SRC, "example_data", "2d_temps.txt"))
    U_levels = np.genfromtxt(os.path.join(SRC, "example_data", "2d_winds.txt"))
    longitudes = np.genfromtxt(os.path.join(SRC, "example_data", "2d_longitudes.txt"))

    return P_levels, T_levels, kzz_levels, U_levels, longitudes




def example_carma(name):
    """
    Initializes and returns a Carma object with an example atmospheric structure
    based on a sonora diamondback run (Morley et al. 2024)/

    The returned carma object includes Pure TiO2 and Mg2SiO4 on TiO2 clouds

    Args:
        name (str): A name identifier for the CARMA model instance.

    Returns:
        Carma: An initialized and configured CARMA model object.
    """
    P_levels, T_levels, kzz_levels, mu_levels = example_levels()

    carma = Carma(name)

    carma.set_physical_params(surface_grav=31600,
                              wt_mol=np.mean(mu_levels))

    carma.set_atmospheric_parameters_from_defaults("Pure H2")
    
    carma.set_stepping(dt=250, output_gap=10, n_tstep=100)

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

    populate_abundances_at_cloud_base(carma)

    return carma