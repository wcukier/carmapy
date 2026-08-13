""" Preparation of the correlated-k opacity tables used by the radiatively
coupled version of CARMA.

CARMA's Fortran engine solves the radiative transfer itself; it does not link
against pyHARP.  What it needs is the Sonora 2020 correlated-k data in a form a
Fortran program can read.  pyHARP owns the download-and-unpack pipeline and
stores the result in a TorchScript container, so the functions here use pyHARP
once, up front, to produce a plain binary file that ``carmapy.exe`` reads at
startup.

``pyharp`` (and therefore ``torch``) is an **optional** dependency, needed only
by :func:`fetch_sonora_table`.  Once the table has been exported, nothing in the
run path imports it -- :func:`read_ck_table` and :func:`sonora_bands` read the
exported file directly.

Notes
-----
The layout of the exported file mirrors the contract pyHARP itself uses, which
was read off ``pyharp/python/sonora/sonora.py`` and
``pyharp/src/opacity/multiband.cpp`` rather than from documentation:

- ``wmin``/``wmax`` are wavenumbers in cm^-1, despite pyHARP's
  ``load_sonora_window`` docstring claiming nm.
- Spectral points are ordered band-major, ``iwave = iband * ng + ig``, matching
  the ``(nband, ng) -> (nband*ng,)`` flatten in ``save_sonora_multiband``.
- ``kappa`` is stored as **ln(cm^2 / molecule)**.  ``multiband.cpp`` converts to
  an attenuation coefficient with ``1e-4 * N_A * exp(kappa) * conc``, where
  ``conc`` is a mole concentration in mol/m^3.
- The temperature axis is **absolute** temperature, not the anomaly convention
  used by pyHARP's line-by-line tables.
- pyHARP interpolates trilinearly over ``(wavenumber, ln(pressure),
  temperature)`` with ``extrapolate=false``, i.e. values are **clamped** at the
  edges of the table.  The Fortran reader must clamp too; see
  ``carma_ckopacity.F90``.
"""

import os
import re
import warnings

import numpy as np

__all__ = [
    "fetch_sonora_table",
    "export_ck_table",
    "read_ck_table",
    "sonora_bands",
    "band_wavelengths",
    "gen_band_mie_table",
    "write_optics_file",
    "SONORA_2020_SOLAR",
]

#: Default Sonora 2020 correlated-k dataset: solar metallicity, solar C/O,
#: 196 spectral windows.  This is the file
#: ``pyharp/examples/example_sonora_2020.yaml`` is built around.
SONORA_2020_SOLAR = "sonora_2020_feh+000_co_100.data.196"

#: Written as the first 8 bytes of an exported table.  Reading it back as a
#: little-endian int64 and getting something else means the file was written on
#: a machine of the opposite endianness.
_CK_MAGIC = 0x43524D41434B3031  # "CARMACK01"
_CK_VERSION = 1


#: Matches the Sonora 2020 dataset naming scheme, whose metallicity and C/O
#: fields are what select the file on Zenodo.
_SONORA_NAME = re.compile(
    r"sonora_2020_feh(?P<feh>[+-]\d{3})_co_(?P<co>\d{3})\.data\.196")


def _download_sonora_archive(name: str) -> str:
    """Download the Sonora 2020 ``.tar.gz`` for ``name`` into the cwd.

    pyHARP's ``load_sonora_data`` opens the archive from the working directory
    and does not fetch it, so this fills that gap using pyHARP's own Zenodo
    URL builder.
    """
    from pyharp.api.fetch_sonora import download_file, get_sonora_data

    match = _SONORA_NAME.fullmatch(name)
    if match is None:
        raise ValueError(
            f"cannot work out the Zenodo URL for {name!r}: expected a name "
            "like 'sonora_2020_feh+000_co_100.data.196'. Download the archive "
            "yourself and place it in the working directory as "
            f"'{name}.tar.gz' to skip this step.")

    archive = f"{name}.tar.gz"
    download_file(get_sonora_data(match["feh"], match["co"]), filename=archive)

    return archive


def fetch_sonora_table(name: str = SONORA_2020_SOLAR,
                       directory: str = ".") -> str:
    """Download and unpack a Sonora 2020 correlated-k dataset, writing the
    TorchScript container pyHARP uses.

    Requires ``pyharp`` to be installed (``pip install carmapy[radiation]``).
    This is a one-time step; :func:`export_ck_table` turns the result into the
    binary CARMA actually reads.

    Both intermediates are kept and reused: the ``.tar.gz`` is only fetched if
    it is not already in ``directory``, and nothing is done at all if the
    ``.pt`` already exists.  The archive is about 39 MB; the ``.pt`` it unpacks
    to is considerably larger.

    Parameters
    ----------
    name : str, optional
        The Sonora dataset name, without the ``.tar.gz`` extension.  Defaults
        to :data:`SONORA_2020_SOLAR`.  The metallicity and C/O fields in the
        name select which file is fetched.
    directory : str, optional
        Directory to work in.  Defaults to the current directory.

    Returns
    -------
    str
        The path to the written ``.pt`` file.
    """
    try:
        from pyharp.sonora import load_sonora_data, save_sonora_multiband
    except ImportError as exc:
        raise ImportError(
            "fetch_sonora_table requires pyharp. Install it with "
            "`pip install carmapy[radiation]`. Note that this is only needed "
            "once, to prepare the opacity tables; running CARMA does not "
            "import pyharp."
        ) from exc

    cwd = os.getcwd()
    try:
        os.chdir(directory)
        pt_path = os.path.abspath(f"{name}.pt")
        if not os.path.exists(pt_path):
            if not os.path.exists(f"{name}.tar.gz"):
                _download_sonora_archive(name)
            data = load_sonora_data(name)
            # clean=False keeps the archive, so a re-export does not re-download.
            save_sonora_multiband(name, data, clean=False)
    finally:
        os.chdir(cwd)

    return pt_path


def export_ck_table(pt_path: str, out_path: str = None) -> str:
    """Convert a pyHARP Sonora ``.pt`` container into the flat binary that
    ``carmapy.exe`` reads.

    Requires ``torch`` (to open the TorchScript container) but not the rest of
    pyHARP.  Like :func:`fetch_sonora_table`, this runs once; the resulting
    file is a static input to every coupled simulation.

    Parameters
    ----------
    pt_path : str
        Path to the ``.pt`` file written by :func:`fetch_sonora_table`.
    out_path : str, optional
        Where to write the binary.  Defaults to ``pt_path`` with the extension
        replaced by ``.ck``.

    Returns
    -------
    str
        The path to the written table.
    """
    try:
        import torch
    except ImportError as exc:
        raise ImportError(
            "export_ck_table requires torch to read the TorchScript container "
            "pyharp writes. Install it with `pip install carmapy[radiation]`."
        ) from exc

    if out_path is None:
        out_path = os.path.splitext(pt_path)[0] + ".ck"

    container = torch.jit.load(pt_path)

    def _get(name):
        return np.asarray(getattr(container, name).detach().numpy(),
                          dtype=np.float64)

    wmin = _get("wmin")            # (nband,)   [cm^-1]
    wmax = _get("wmax")            # (nband,)   [cm^-1]
    gauss_pts = _get("gauss_pts")  # (ng,)
    gauss_wts = _get("gauss_wts")  # (ng,)
    wavenumber = _get("wavenumber")  # (nband*ng,) [cm^-1]
    weights = _get("weights")      # (nband*ng,)
    pres = _get("pres")            # (npres,)   [Pa]
    temp = _get("temp")            # (ntemp,)   [K]
    kappa = _get("kappa")          # (nband*ng, npres, ntemp) [ln(cm^2/molec)]

    nband = wmin.size
    ng = gauss_pts.size
    npres = pres.size
    ntemp = temp.size

    # Guard the two shape assumptions the Fortran side bakes in: band-major
    # ordering, and kappa already flattened over (band, g).
    if wavenumber.size != nband * ng:
        raise ValueError(
            f"wavenumber has {wavenumber.size} entries but nband*ng = "
            f"{nband * ng}; the (nband, ng) flatten convention has changed")
    if kappa.shape != (nband * ng, npres, ntemp):
        raise ValueError(
            f"kappa has shape {kappa.shape}, expected "
            f"{(nband * ng, npres, ntemp)}")

    with open(out_path, "wb") as f:
        np.array([_CK_MAGIC], dtype=np.int64).tofile(f)
        np.array([_CK_VERSION, nband, ng, npres, ntemp],
                 dtype=np.int64).tofile(f)
        for arr in (wmin, wmax, gauss_pts, gauss_wts,
                    wavenumber, weights, pres, temp):
            arr.astype(np.float64).tofile(f)
        # Fortran reads this as kappa(nband*ng, npres, ntemp), so write it in
        # column-major order.
        kappa.astype(np.float64).ravel(order="F").tofile(f)

    return out_path


def read_ck_table(path: str) -> dict:
    """Read back a table written by :func:`export_ck_table`.

    Used by the validation tests, which compare the Fortran interpolation
    against the same table interpolated in Python.  Needs neither pyharp nor
    torch.

    Parameters
    ----------
    path : str
        Path to an exported ``.ck`` file.

    Returns
    -------
    dict
        Keys ``wmin``, ``wmax``, ``gauss_pts``, ``gauss_wts``, ``wavenumber``,
        ``weights``, ``pres``, ``temp``, ``kappa``.
    """
    with open(path, "rb") as f:
        magic = np.fromfile(f, dtype=np.int64, count=1)[0]
        if magic != _CK_MAGIC:
            raise ValueError(
                f"{path} is not a carmapy ck table, or was written on a "
                f"machine of the opposite endianness (magic = {magic:#x})")

        version, nband, ng, npres, ntemp = np.fromfile(f, dtype=np.int64,
                                                       count=5)
        if version != _CK_VERSION:
            raise ValueError(f"{path} is version {version}, this build reads "
                             f"version {_CK_VERSION}")

        def _read(n):
            return np.fromfile(f, dtype=np.float64, count=int(n))

        out = {
            "wmin":       _read(nband),
            "wmax":       _read(nband),
            "gauss_pts":  _read(ng),
            "gauss_wts":  _read(ng),
            "wavenumber": _read(nband * ng),
            "weights":    _read(nband * ng),
            "pres":       _read(npres),
            "temp":       _read(ntemp),
        }
        out["kappa"] = _read(nband * ng * npres * ntemp).reshape(
            (nband * ng, npres, ntemp), order="F")

    return out


def sonora_bands(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Return the band edges of an exported table, in wavenumber [cm^-1].

    These are the wavelengths the Mie tables must be regenerated on -- cloud
    optical properties are evaluated per band, not per g-point, since they vary
    smoothly across a band.

    Parameters
    ----------
    path : str
        Path to an exported ``.ck`` file.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        ``(wmin, wmax)``, each of length ``nband``, in cm^-1.
    """
    table = read_ck_table(path)
    return table["wmin"], table["wmax"]


def band_wavelengths(path: str) -> np.ndarray:
    """Return one representative wavelength per band, in cm.

    Cloud optical properties are evaluated per *band* rather than per g-point:
    the g-points within a band share a wavelength range, and Mie properties
    vary smoothly across it.  That is an 8x saving on the dominant cost of the
    cloud-optics sum and introduces no error the band discretization does not
    already have.

    The representative point is the band's mean wavenumber, converted to a
    wavelength.  The returned array is in the table's own band order (ascending
    wavenumber, hence *descending* wavelength) -- the Fortran indexes it
    directly, so the order must not be sorted.

    Parameters
    ----------
    path : str
        Path to an exported ``.ck`` file.

    Returns
    -------
    np.ndarray
        ``(nband,)`` wavelengths in cm.
    """
    wmin, wmax = sonora_bands(path)
    return 2.0 / (wmin + wmax)


def bin_radii(rmin: float, nbin: int, rmrat: float = 2.0) -> np.ndarray:
    """The CARMA bin radii for a group, in cm.

    Mirrors ``setupbins.F90``: bin masses go as ``rmrat**(j-1)``, so for
    spheres the radii go as ``rmrat**((j-1)/3)``.  Note ``rmrat`` is a *volume*
    ratio -- passing it straight to :func:`carmapy.mie.gen_mie_table`, whose
    ``ratio`` is a *radius* ratio, would give the wrong grid.

    Parameters
    ----------
    rmin : float
        Radius of the smallest bin [cm].
    nbin : int
        Number of bins.
    rmrat : float, optional
        Volume ratio between successive bins.  Defaults to 2, which is what
        ``carma_carmapy.F90`` passes to ``CARMAGROUP_Create``.

    Returns
    -------
    np.ndarray
        ``(nbin,)`` radii in cm.
    """
    return rmin * rmrat ** (np.arange(nbin) / 3.0)


def _species_of(group_name: str) -> str:
    """The condensing species of a group, from its CARMApy name.

    Group names are ``"Pure {species}"`` or ``"{mantle} on {seed}"``; in both
    cases the optical properties are those of the outer material.  Same
    convention as :func:`carmapy.results._get_cloud_opacities`.
    """
    parts = group_name.split()
    return parts[-1] if parts[0] == "Pure" else parts[0]


def gen_band_mie_table(group_name: str,
                       rmin: float,
                       nbin: int,
                       ck_path: str,
                       out_path: str,
                       rmrat: float = 2.0,
                       indices_path: str = None,
                       overwrite: bool = False) -> np.ndarray:
    """Generate a group's Mie table on the CARMA bin grid and the ck band set.

    Unlike the bundled tables -- which are on their own radius grid and must be
    resampled with a spline -- this writes the table on *exactly* the radii and
    wavelengths the run uses, so the Fortran needs no interpolation at all.

    The Sonora band set is much wider than any bundled refractive index file
    (196 bands span roughly 0.26-325 um), so n and k are clamped at the ends of
    the index data.  How far that clamping reaches is reported as a warning
    rather than left implicit: for the far-infrared bands it is harmless, since
    the grains are deep in the Rayleigh limit and the extinction is collapsing
    anyway, but it is a real extrapolation and worth seeing.

    Parameters
    ----------
    group_name : str
        CARMApy group name, e.g. ``"Mg2SiO4 on TiO2"``.
    rmin, nbin, rmrat
        The group's bin grid; see :func:`bin_radii`.
    ck_path : str
        Path to an exported ``.ck`` file, which fixes the band set.
    out_path : str
        Where to write the table.
    indices_path : str, optional
        Override the bundled refractive index file.
    overwrite : bool, optional
        Recompute even if ``out_path`` already exists.  Defaults to False, so
        repeated runs reuse the table.

    Returns
    -------
    np.ndarray
        ``(nbin, nband, 3)`` of ``Q_ext``, ``Q_sca``, ``g``, in band order.
    """
    from . import mie

    species = _species_of(group_name)
    lambdas = band_wavelengths(ck_path)
    radii = bin_radii(rmin, nbin, rmrat)

    if indices_path is None:
        indices_path = mie.default_indices_path(species)

    _, _, lam_lo, lam_hi = mie._load_cloud_indices(indices_path, clamp=True)
    lo, hi = lambdas.min(), lambdas.max()

    tol = 1e-6
    reach = []
    if lo < lam_lo * (1 - tol):
        reach.append(f"{lam_lo / lo:.2g}x into the blue")
    if hi > lam_hi * (1 + tol):
        reach.append(f"{hi / lam_hi:.2g}x into the red")
    if reach:
        nout = int(np.count_nonzero((lambdas < lam_lo * (1 - tol))
                                    | (lambdas > lam_hi * (1 + tol))))
        warnings.warn(
            f"{species}: refractive index data covers "
            f"{lam_lo * 1e4:.3g}-{lam_hi * 1e4:.4g} um, but the band set spans "
            f"{lo * 1e4:.3g}-{hi * 1e4:.4g} um. Holding n and k at their "
            f"endpoint values for {nout} of {lambdas.size} bands, reaching "
            f"{' and '.join(reach)}.",
            stacklevel=2)

    if overwrite or not os.path.exists(out_path):
        # gen_mie_table's `ratio` is a radius ratio; rmrat is a volume ratio.
        mie.gen_mie_table(species,
                          min_r=rmin,
                          ratio=rmrat ** (1.0 / 3.0),
                          n_r=nbin,
                          lambdas=lambdas,
                          indices_path=indices_path,
                          out_path=out_path,
                          clamp_indices=True)

    # The writer loops radius-major with `lambdas` innermost, in the order
    # given, so the file can be reshaped without sorting -- which matters,
    # since band order is not wavelength order.
    raw = np.genfromtxt(out_path, delimiter="\t", names=True)
    nband = lambdas.size
    if raw.size != nbin * nband:
        raise ValueError(
            f"{out_path} has {raw.size} rows, expected {nbin * nband} "
            f"({nbin} bins x {nband} bands). Delete it to regenerate.")

    table = np.empty((nbin, nband, 3))
    table[:, :, 0] = raw["Q_ext"].reshape(nbin, nband)
    table[:, :, 1] = raw["Q_sca"].reshape(nbin, nband)
    table[:, :, 2] = raw["g"].reshape(nbin, nband)

    # Cheap guard against a stale or mismatched cache.
    if not np.allclose(raw["rcm"].reshape(nbin, nband)[:, 0], radii,
                       rtol=1e-9):
        raise ValueError(
            f"{out_path} is on a different radius grid than this group "
            f"(rmin={rmin:g}, nbin={nbin}, rmrat={rmrat:g}). Delete it to "
            f"regenerate.")

    return table


def write_optics_file(path: str, tables: dict) -> str:
    """Write the per-group cloud optics CARMA reads at startup.

    One row per (group, band, bin), carrying extinction efficiency, single
    scattering albedo and asymmetry factor -- the ``(NWAVE, NBIN)`` arrays
    ``CARMAGROUP_Create`` already accepts.  Because the tables are on the run's
    own bin grid and band set, these are static for the whole run; the only
    per-step work in the Fortran is the number-density-weighted sum over bins.

    Parameters
    ----------
    path : str
        Where to write the file.
    tables : dict
        Maps group name to the ``(nbin, nband, 3)`` array returned by
        :func:`gen_band_mie_table`.  Iteration order sets the group index, so
        it must match the order groups are written to ``groups.txt``.

    Returns
    -------
    str
        The path written.
    """
    names = list(tables)
    shapes = {tables[n].shape[:2] for n in names}
    if len(shapes) != 1:
        raise ValueError(f"groups disagree on (nbin, nband): {shapes}")
    nbin, nband = shapes.pop()

    with open(path, "w+") as f:
        f.write(f"{len(names)}\t{nband}\t{nbin}\n")
        f.write("igroup\tiwave\tibin\tqext\tssa\tasym\n")
        for igroup, name in enumerate(names):
            table = tables[name]
            qext = table[:, :, 0]
            qsca = table[:, :, 1]
            asym = table[:, :, 2]
            # Q_ext is zero only where the Mie call underflowed; leaving ssa at
            # zero there is right, since nothing is being scattered either.
            ssa = np.divide(qsca, qext, out=np.zeros_like(qsca),
                            where=qext > 0.0)
            for iwave in range(nband):
                for ibin in range(nbin):
                    f.write(f"{igroup + 1}\t{iwave + 1}\t{ibin + 1}\t"
                            f"{qext[ibin, iwave]:.17e}\t"
                            f"{ssa[ibin, iwave]:.17e}\t"
                            f"{asym[ibin, iwave]:.17e}\n")

    return path
