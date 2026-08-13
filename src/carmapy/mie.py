import numpy as np
import os

import scipy.integrate as _integrate

# PyMieScatt does `from scipy.integrate import trapz`, which scipy 1.14 removed
# in favour of the identically-behaved `trapezoid`. Without this alias importing
# PyMieScatt raises, taking this whole module with it.
if not hasattr(_integrate, "trapz"):
    _integrate.trapz = _integrate.trapezoid

import PyMieScatt as ps

from scipy.interpolate import interp1d
from .constants import group_dict

_SRC = os.path.dirname(os.path.dirname(__file__))
defualt_lambdas = np.linspace(1e-4, 3e-3, 1000)


def default_indices_path(species):
    """Path to the bundled complex refractive index file for a species.

    The filenames vary by source ("TiO2[s].txt", "Mg2SiO4_amorph.txt",
    "Fe_complex.txt", ...), so the mapping lives in ``group_dict`` rather than
    being derivable from the species name.
    """
    try:
        filename = group_dict[species]["opacity_files"]
    except KeyError:
        raise KeyError(
            f"no bundled refractive index file for {species!r}; pass "
            f"indices_path explicitly. Known species: "
            f"{sorted(k for k, v in group_dict.items() if 'opacity_files' in v)}"
        ) from None

    return os.path.join(_SRC, "mie_tables", "indices", filename)


def _load_cloud_indices(path, clamp=False):
    """Build n(lambda) and k(lambda) interpolators from a refractive index file.

    By default the interpolators raise outside the tabulated wavelength range.
    With ``clamp=True`` they hold the endpoint values instead, which is what the
    radiative-coupling path needs -- the Sonora band set is far wider than any
    of the bundled index files.  Callers that clamp are expected to tell the
    user how far they are extrapolating; see
    :func:`carmapy.radiation.gen_band_mie_table`.

    Returns ``(n_interp, k_interp, lambda_min, lambda_max)`` with the bounds in
    cm, so callers can report the extrapolation without re-reading the file.
    """
    data = np.genfromtxt(path,
                        comments="#")
    wavelengths = data[:, 0] * 1e-4 # convert to cm

    # The files are ordered by wavelength, but don't rely on it -- interp1d
    # sorts internally and the endpoint values must match that order.
    order = np.argsort(wavelengths)
    wavelengths = wavelengths[order]
    n_vals = data[order, 1]
    k_vals = data[order, 2]

    if clamp:
        n_interp = interp1d(wavelengths, n_vals, bounds_error=False,
                            fill_value=(n_vals[0], n_vals[-1]))
        k_interp = interp1d(wavelengths, k_vals, bounds_error=False,
                            fill_value=(k_vals[0], k_vals[-1]))
    else:
        n_interp = interp1d(wavelengths, n_vals)
        k_interp = interp1d(wavelengths, k_vals)

    return n_interp, k_interp, wavelengths[0], wavelengths[-1]

def gen_mie_table(species: str,
                  min_r: float = 1e-8,
                  ratio: float = 2,
                  n_r: int = 80,
                  lambdas : np.ndarray = defualt_lambdas,
                  indices_path : str = None,
                  out_path : str = None,
                  clamp_indices : bool = False):
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
    out_path : str, optional
        Where to write the table.  By default writes
        ``src/mie_tables/{species}_user.dat`` inside the installed package,
        which is where :func:`carmapy.results._get_cloud_opacities` looks for
        it.  Pass an explicit path to keep a table with a run instead.
    clamp_indices : bool, optional
        If True, hold n and k at their endpoint values for wavelengths outside
        the refractive index file rather than raising.  Needed when ``lambdas``
        is wider than the index data; the caller should report the
        extrapolation.  Defaults to False.

    Returns
    -------
    str
        The path the table was written to.
    """
    if indices_path is None:
        indices_path = default_indices_path(species)
    n_interp, k_interp, _, _ = _load_cloud_indices(indices_path,
                                                   clamp=clamp_indices)

    rs = np.array([min_r * ratio ** ir for ir in range(n_r)])

    if out_path is None:
        out_path = os.path.join(_SRC, "mie_tables", f"{species}_user.dat")

    with open(out_path, "w+") as f:
        f.write("r[cm]\tλ[cm]\tQ_ext\tQ_sca\tg\n")
        for ibin in range(n_r):
            for il in range(len(lambdas)):
                m = n_interp(lambdas[il]) + 1j*k_interp(lambdas[il])

                mie_res = ps.AutoMieQ(m,
                                lambdas[il]*1e7,
                                2*rs[ibin]*1e7,
                                asDict=True)
                f.write(f"{rs[ibin]}\t{lambdas[il]}\t{mie_res['Qext']}"
                        f"\t{mie_res['Qsca']}\t{mie_res['g']}\n")

            f.flush()

    return out_path