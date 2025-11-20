import numpy as np
import PyMieScatt as ps
import os

from scipy.interpolate import interp1d
from .constants import gas_dict

_SRC = os.path.dirname(os.path.dirname(__file__))
defualt_lambdas = np.linspace(1e-4, 3e-3, 1000)


def _load_cloud_indices(path):
    data = np.genfromtxt(path,
                        comments="#")
    wavelengths = data[:, 0] * 1e-4 # convert to cm
    n_interp = interp1d(wavelengths, data[:, 1])
    k_interp = interp1d(wavelengths, data[:, 2])
    return n_interp, k_interp

def gen_mie_table(species: str,
                  min_r: float = 1e-8,
                  ratio: float = 2,
                  n_r: int = 80,
                  lambdas : np.ndarray = defualt_lambdas, 
                  indices_path : str = None):
    """(Re-)generate mie tables for the specified species.  Note the units for
    the indices of refraction table has wavelength in microns, not cm.

    Parameters
    ----------
    species : str
        The cloud species to generate mie tables for
    min_r : float, optional
        The minimum particle radius to generate the table for [cm], 
        by default 1e-8
    ratio : float, optional
        The ratio between sucessive radii, by default 2
    n_r : int, optional
        The number of radii to generate the table for, by default 80
    lambdas : np.ndarray, optional
        The wavelengths to calculate mie scattering at [cm], by default 
        1000 evenly spaced wavelengths between 1e-4 and 3e-3 cm
    indices_path : str, optional
        File path that points to the complex indices of refraction of the 
        specified species.  The columns must be whitespace separated and have 
        columns of wavelength (micron), n, k, under the convention m = n + ik.
        By default will use the carmapy default files
    """
    if indices_path is None:
        indices_path = os.path.join(os.path.dirname(_SRC),
                                    "mie_tables",
                                    "indices",
                                    gas_dict[species]["opacity_files"])
    n_interp, k_interp = _load_cloud_indices(indices_path)

    rs = np.array([min_r * ratio ** ir for ir in range(len(n_r))])

    with open(os.path.join(_SRC, "mie_tables", f"{species}_user.dat"), "w+") as f:
        f.write("r[cm]\tλ[cm]\tQ_ext\tQ_sca\tg\n")
        for ibin in range(n_r):
            for il in range(len(lambdas)):
                m = n_interp(lambdas[il]) + 1j*k_interp(lambdas[il])

                mie_res = ps.AutoMieQ(m, 
                                lambdas*1e7, 
                                2*rs*1e7, 
                                asDict=True)
                f.write(f"{rs[ibin]}\t{lambdas[il]}\t{mie_res['Qext']}"
                         "\t{mie_res['Qsca']}\t{mie_res['g']}\n")

            f.flush()