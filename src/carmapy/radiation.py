""" Preparation of the correlated-k opacity tables used by the radiatively
coupled version of CARMA.

pyHARP downloads and unpacks the Sonora 2020 correlated-k data into a
TorchScript container; the functions here convert that once into a plain
binary that ``carmapy.exe`` reads at startup.  ``pyharp`` (and ``torch``) are
optional dependencies, needed only to fetch and export the table.

Notes
-----
The exported layout follows pyHARP's (``pyharp/python/sonora/sonora.py`` and
``pyharp/src/opacity/multiband.cpp``):

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
    "Radiation",
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

    Needs neither pyharp nor torch.

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

    Cloud optical properties are evaluated per band rather than per g-point,
    at the band's mean wavenumber.  The array is in the table's band order
    (descending wavelength), which the Fortran indexes directly, so it must
    not be sorted.

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

    The table is on exactly the radii and wavelengths the run uses, so the
    Fortran does no interpolation.  The Sonora bands (~0.26-325 um) are wider
    than any bundled refractive index file, so n and k are held at their
    endpoint values outside it, with a warning saying how far.

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

    # Rows are radius-major with lambdas innermost, in band order.
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
    ``CARMAGROUP_Create`` accepts.

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
            # ssa = 0 where Q_ext underflowed to zero.
            ssa = np.divide(qsca, qext, out=np.zeros_like(qsca),
                            where=qext > 0.0)
            for iwave in range(nband):
                for ibin in range(nbin):
                    f.write(f"{igroup + 1}\t{iwave + 1}\t{ibin + 1}\t"
                            f"{qext[ibin, iwave]:.17e}\t"
                            f"{ssa[ibin, iwave]:.17e}\t"
                            f"{asym[ibin, iwave]:.17e}\n")

    return path


# LLM-GENERATED — UNVERIFIED (tier: simulation)
# Remove this banner only after completing the checks in Assistant/POLICY.md §2.
class Radiation:
    """Settings for a radiatively coupled temperature profile.

    Built by :meth:`carmapy.Carma.set_radiation`, which documents the
    physics, and held as ``Carma.radiation``; ``None`` there is a
    fixed-temperature run.  The constructor checks the ranges but not that
    the ck table exists, so a run whose table has since moved still loads.

    Parameters
    ----------
    ck_table_path : str
        Correlated-k table written by :func:`export_ck_table`.
    t_int : float
        Internal temperature [K].
    t_irr, t_star, mu0, w_surf, cloud_rad : optional
        As in :meth:`carmapy.Carma.set_radiation`.
    self_consistent_kzz, kzz_mixl_scale, accel : optional
        As in :meth:`carmapy.Carma.set_radiation`.
    gap_max : int, optional
        ``rad_gap_max`` of :meth:`carmapy.Carma.set_radiation`.
    dT_tol : float, optional
        ``rad_dT_tol`` of :meth:`carmapy.Carma.set_radiation`.
    dtau_tol : float, optional
        ``rad_dtau_tol`` of :meth:`carmapy.Carma.set_radiation`.
    """

    ck_table_path: str
    """ Correlated-k table read by the radiative transfer """

    t_int: float
    """ Internal temperature [K] """

    t_irr: float
    """ Irradiation temperature [K]; 0 is an isolated object """

    t_star: float
    """ Blackbody temperature setting the incident spectrum [K]; 0 uses t_irr """

    mu0: float
    """ Cosine of the incidence angle """

    w_surf: float
    """ Surface albedo """

    cloud_rad: bool
    """ Whether the clouds contribute opacity """

    self_consistent_kzz: bool
    """ Whether Kzz is recomputed from mixing length theory on each solve """

    kzz_mixl_scale: float
    """ Multiplier on the mixing length """

    accel: float
    """ Speed-up of the radiative clock relative to the microphysical one """

    gap_max: int
    """ Maximum number of timesteps between full radiative solves """

    dT_tol: float
    """ Accumulated max|dT| since the last solve that forces a new one [K] """

    dtau_tol: float
    """ Fractional change in column cloud cross section that forces a solve """

    #: The &radiation entries written for an uncoupled run.
    UNCOUPLED_NML = {
        "do_radiation":         0,
        "ck_table_file":        "",
        "t_int":                0.,
        "t_irr":                0.,
        "t_star":               0.,
        "rad_mu0":              0.5,
        "rad_w_surf":           0.,
        "rad_accel":            1.,
        "rad_dT_tol":           1.,
        "rad_dtau_tol":         0.02,
        "rad_gap_max":          100,
        "kzz_mode":             1,
        "kzz_mixl_scale":       1.,
        "cloud_rad":            1,
    }

    def __init__(self,
                 ck_table_path: str,
                 t_int: float,
                 t_irr: float = 0.0,
                 t_star: float = 0.0,
                 mu0: float = 0.5,
                 w_surf: float = 0.0,
                 cloud_rad: bool = True,
                 self_consistent_kzz: bool = True,
                 kzz_mixl_scale: float = 1.0,
                 accel: float = 1.0,
                 gap_max: int = 100,
                 dT_tol: float = 1.0,
                 dtau_tol: float = 0.02):
        if t_int <= 0:
            raise ValueError("t_int must be positive")

        if t_irr < 0:
            raise ValueError("t_irr must be non-negative (0 = isolated object)")

        if t_star < 0:
            raise ValueError("t_star must be non-negative (0 = shape like t_irr)")

        if not 0 < mu0 <= 1:
            raise ValueError("mu0 is a cosine and must be in (0, 1]")

        if not 0 <= w_surf < 1:
            raise ValueError("w_surf is an albedo and must be in [0, 1)")

        if gap_max < 1:
            raise ValueError("gap_max must be at least 1")

        if kzz_mixl_scale <= 0:
            raise ValueError("kzz_mixl_scale must be positive, got "
                             f"{kzz_mixl_scale}")

        self.ck_table_path = ck_table_path
        self.t_int = t_int
        self.t_irr = t_irr
        self.t_star = t_star
        self.mu0 = mu0
        self.w_surf = w_surf
        self.cloud_rad = cloud_rad
        self.self_consistent_kzz = self_consistent_kzz
        self.kzz_mixl_scale = kzz_mixl_scale
        self.accel = accel
        self.gap_max = gap_max
        self.dT_tol = dT_tol
        self.dtau_tol = dtau_tol

    def to_nml(self) -> dict:
        """The &radiation entries for this coupling, without the adiabat's."""
        return {
            "do_radiation":         1,
            "ck_table_file":        os.path.abspath(self.ck_table_path),
            "t_int":                self.t_int,
            "t_irr":                self.t_irr,
            "t_star":               self.t_star,
            "rad_mu0":              self.mu0,
            "rad_w_surf":           self.w_surf,
            "rad_accel":            self.accel,
            "rad_dT_tol":           self.dT_tol,
            "rad_dtau_tol":         self.dtau_tol,
            "rad_gap_max":          self.gap_max,
            "kzz_mode":             int(self.self_consistent_kzz),
            "kzz_mixl_scale":       self.kzz_mixl_scale,
            "cloud_rad":            int(self.cloud_rad),
        }

    @classmethod
    def from_nml(cls, rad) -> "Radiation":
        """Rebuild from a parsed &radiation group, or return None when it has
        radiation off.  Keys are matched case-insensitively, as f90nml reads
        them."""
        rad = {k.lower(): v for k, v in rad.items()}
        if not rad.get("do_radiation", 0):
            return None
        return cls(rad["ck_table_file"],
                   rad["t_int"],
                   t_irr=rad["t_irr"],
                   t_star=rad["t_star"],
                   mu0=rad["rad_mu0"],
                   w_surf=rad["rad_w_surf"],
                   cloud_rad=bool(rad["cloud_rad"]),
                   self_consistent_kzz=bool(rad["kzz_mode"]),
                   kzz_mixl_scale=rad["kzz_mixl_scale"],
                   accel=rad["rad_accel"],
                   gap_max=int(rad["rad_gap_max"]),
                   dT_tol=rad["rad_dt_tol"],
                   dtau_tol=rad["rad_dtau_tol"])
