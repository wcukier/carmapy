import numpy as np
from carmapy.constants import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button, Slider
import matplotlib as mpl
from itertools import cycle
from scipy.signal import periodogram
from scipy.interpolate import interp1d
from carmapy.chemistry import get_fastchem_abundances
import os
# import PyMieScatt as ps
from scipy.interpolate import RectBivariateSpline

from collections.abc import Callable

# petroff 10 color cycle
petroff10 = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6",
              "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]

_SRC = os.path.dirname(os.path.dirname(__file__))

MICRON_TO_CM = 1e-4

class Results:
    """ An object which stores the results of a CARMA simulation.

    Parameters
    ----------
    carma : Carma
        The carma simulation to load results from
    read_diag : boolean, optional
        If true reads in the microphysical rates and core mass fraction.
        Defaults to False.

    Notes
    ------
    The ``results.clouds`` dictionary is indexed by names of the cloud
    species (see ``results.group_names``) and stores the following values:

        - ``results.clouds['r']`` is a 1D numpy array of length NBIN which 
          stores the radii of the cloud particle bins in cm
        - ``results.clouds['rmass']`` is a 1D numpy array of length NBIN which 
          stores the mass of the cloud particle bins in g
        - ``results.clouds['numden']`` is a 3D numpy array of shape 
          (NZ, NBIN, NT) which stores the number density of each cloud species 
          in each size bin at each time step in units of cm⁻³
        - ``results.clouds['coremass_frac']`` is a 3D numpy array of 
           shape (NZ, NBIN, NT) which stores the fraction of mass in each
           bin stored in the particle core.
           Only stored if ``read_diag=True``.
        - ``results.clouds['nuc_gain_rate']`` is a 3D numpy array of 
           shape (NZ, NBIN, NT) which stores the rate at which particles
           are added to the bin via nucleation [particles/s/cm³].
           Only stored if ``read_diag=True``.
         - ``results.clouds['nuc_loss_rate']`` is a 3D numpy array of 
           shape (NZ, NBIN, NT) which stores the rate at which particles are 
           removed from a bin due to nucleation [particles/s/cm³].
           Only stored if ``read_diag=True``.
        - ``results.clouds['grow_gain_rate']`` is a 3D numpy array of 
           shape (NZ, NBIN, NT) which stores the rate at which 
           particles are added to a bin due to condensational growth 
           [particles/s/cm³].
           Only stored if ``read_diag=True``.
        - ``results.clouds['grow_loss_rate']`` is a 3D numpy array of 
           shape (NZ, NBIN, NT) which stores the rate at which 
           particles are removed from a bin due to condensational growth
           [particles/s/cm³].
           Only stored if ``read_diag=True``.
        - ``results.clouds['evap_gain_rate']`` is a 3D numpy array of 
           shape (NZ, NBIN, NT) which stores the rate at which particles are 
           added to a bin due to evaporation [particles/s/cm³]   .
           Only stored if ``read_diag=True``.
        - ``results.clouds['evap_loss_rate']`` is a 3D numpy array of 
           shape (NZ, NBIN, NT) which stores the rate at which particles are 
           removed from a bin due to evaporation [particles/s/cm³].
           Only stored if ``read_diag=True``.         

    In a 2-D CARMApy run two additional items are saved:
        - ``results.gases_2d`` stores the same info as ``results.gases`` but
          stores (NZ, NLONGITUDE) arrays instead of (NZ, NT) arrays
        - ``results.longitude_map`` takes in a 3-D array of shape (NZ, NBIN, NT)
          and transforms it to an array of shape `(NZ, NBIN, NLONGITUDE)` where 
          each longitude bin in the average of all timesteps corresponding to 
          that longitude.  This function is designed to work on the 
          ``results.clouds["numden"]`` array as well as any of the microphysical
          rates arrays in ``results.clouds``.
    """



    carma : "Carma"
    """ The carma simulation that generated these results """

    # rmass : np.ndarray
    # """ Mass of particles in each radius bin [g] (NBIN x NGROUP) """

    # r : np.ndarray
    # """ Radius of particules in each radius bin [cm] (NBIN x NGROUP) """

    # gas_abund : np.ndarray
    # """ The number mixing rations of the gases (NZ x NGAS x NT) """

    sat_vp : np.ndarray
    """ The saturation mixing ratios of the gases (NZ, NGAS, NT, [NLONG]) """

    ts : np.ndarray
    """ Time elapsed since start of simulation at each step [s] (NT)"""

    P : np.ndarray
    """ Pressure centers alias [barye] (NZ) """

    Z : np.ndarray
    """ Altitude centers alias [cm] (NZ) """

    T : np.ndarray
    """ Temperature centers alias [K] (NZ, [NLONG]) """

    T_history : np.ndarray
    """ Temperature centers at each output of a radiatively coupled run [K]
    (NZ, NT). None if the run was not radiatively coupled. """

    dTdt : np.ndarray
    """ Radiative heating rate at each output [K/s] (NZ, NT).  None if the run
    was not radiatively coupled. """

    flux_net : np.ndarray
    """ Net upward longwave flux at the cell levels [W/m²] (NZ+1, NT).  None if
    the run was not radiatively coupled. """

    F_TOA : np.ndarray
    """ Emergent flux at the top of the atmosphere [W/m²] (NT).  Converges to
    ``sigma * teff**4``.  None if the run was not radiatively coupled. """

    convective : np.ndarray
    """ Which layers the convective adjustment left neutrally stable, at each
    output (NZ, NT), bool.  There may be more than one contiguous zone; see
    ``convective_zones``.  None if the run was not radiatively coupled. """

    nzone : np.ndarray
    """ Number of contiguous convective zones at each output (NT).  None if the
    run was not radiatively coupled. """

    conv_resid : np.ndarray
    """ Largest fractional superadiabaticity the convective adjustment left
    behind at each output (NT).  The adjustment stops after a fixed number of
    sweeps and relaxes the remainder over subsequent steps, so a small non-zero
    value is normal; a growing one means the profile is moving faster than the
    adjustment can follow.  None if the run was not radiatively coupled. """

    rad_interval : np.ndarray
    """ Realized number of timesteps per full radiative solve (NT).  1 means the
    adaptive cadence is saving nothing.  None if the run was not radiatively
    coupled. """

    group_names : list[str]
    """ The name of each of the simulated groups """

    gas_names :  list[str]
    """ The name of each of the simulated gases """

    dt_timestep : float
    """ The length of time between each output [s] """

    path : str
    """ The path to the CARMA output files """
    
    gases : dict[str, np.ndarray]
    """ The number mixing ratio for each of the gases """

    gases_2d : dict[str, np.ndarray]
    """ The number mixing ratio for each of the gases by longitude """

    clouds : dict[str, dict[str, np.ndarray]]
    """ A dictionary storing results for each cloud species (see notes) """

    longitude_map : Callable
    """ A function to transform arrays to be by longitude (see notes) """

    def __init__(self, carma: "Carma", read_diag=False) -> None:
        path = carma.name
        path_end = os.path.basename(path)
        file_path = os.path.join(path, f"{path_end}.txt")

        f = open(file_path)

        (NZ,
         NGROUP,
         NELEM,
         NBIN,
         NGAS,
         nstep,
         iskip) = np.array(f.readline().split(), dtype=int)

        
        if ((NZ != carma.NZ) +
            (NGROUP != len(carma.groups))+
            (NELEM != len(carma.elems))+
            (NBIN != carma.NBIN) +
            (NGAS != len(carma.gases))+
            (nstep - 1 != carma.n_tstep)+
            (iskip != carma.output_gap)
        ):
            raise ValueError(f"Output file inconsistent with carma run")
        
        r = np.zeros((NBIN, NGROUP))
        rmass = np.zeros((NBIN, NGROUP))
        
        for i in range(NGROUP):
            for j in range(NBIN):
                (_,
                 _,
                 r[j, i],
                 rmass[j, i],
                 _,
                 _,
                 _) = np.array(f.readline().split(), dtype=float)
            
        
        kzz = np.zeros(NZ)
        P = np.zeros(NZ)
        T = np.zeros(NZ)
        Z = np.zeros(NZ)

        for i in range(NZ):

            (_,
             Z[i],
             _,
             P[i],
             T[i],
             kzz[i]) = np.array(f.readline().split(), dtype=float)

        f.readline()
        f.readline()

        for j in range(NBIN):
            for i in range(NZ):
                f.readline()

        NT = int(nstep/iskip)

        numden = np.zeros((NZ, NELEM, NBIN, NT))
        gas_abund = np.zeros((NZ, NGAS, NT))
        sat_vp = np.zeros((NZ, NGAS, NT))
        ts = np.zeros(NT)
        if carma.is_2d: step = np.zeros(NT, dtype=int)

        block_rows = NBIN * NZ
        gas_idx = NELEM + 2 + 2 * np.arange(NGAS)
        sat_idx = NELEM + 3 + 2 * np.arange(NGAS)
        ncols_block = None  # detected from first data line

        for it in range(NT):
            header = f.readline()
            if not header:
                break

            if carma.is_2d:
                parts = header.split()
                if len(parts) < 5:
                    break
                ts[it] = float(parts[0])
                step[it] = int(float(parts[3]))
            else:
                ts[it] = float(header)

            lines = [f.readline() for _ in range(block_rows)]
            if not lines[-1]:
                break
            if ncols_block is None:
                ncols_block = len(lines[0].split())
            flat = np.fromstring(' '.join(lines), sep=' ', dtype=float)
            if flat.size != block_rows * ncols_block:
                break
            block = flat.reshape(NBIN, NZ, ncols_block)

            # numden[iz, ielem, ibin, it] = block[ibin, iz, 2+ielem]
            numden[:, :, :, it] = block[:, :, 2:2+NELEM].transpose(1, 2, 0)
            # gas_abund/sat_vp don't depend on ibin in the original loop;
            # all NBIN rows for the same (iz,it) carry identical values.
            gas_abund[:, :, it] = block[0, :, gas_idx].T
            sat_vp[:, :, it]    = block[0, :, sat_idx].T

        numden_groups = np.zeros((NZ, NGROUP, NBIN, NT))
        
        
        for i, key in enumerate(carma.groups.keys()):
            group = carma.groups[key]
            if group.mantle:
                g = group.igroup-1
                e = group.mantle.ielem-1
                numden_groups[:, g, :, :] = numden[:, e, :, :]
            else:
                g = group.igroup-1
                e = group.core.ielem-1
                numden_groups[:, g, :, :] = numden[:, e, :, :]

            if np.any( numden_groups[:, group.igroup-1, :, :] < 0):
                raise RuntimeError("Error in reading in number densities")
        
        f.close()
        
        n_tstep = it + 1

        self.carma = carma
        self.rmass = rmass
        self.r = r
        self.numden = numden_groups[:,:,:,:it+1]
        self.gas_abund = gas_abund[:,:,:it+1]
        self.sat_vp = sat_vp[:,:,:it+1]
        self.ts = ts[:it+1]
        self.P = P
        self.Z = Z
        self.T = T
        self.group_names = list(carma.groups.keys())
        self.gas_names = list(carma.gases.keys())
        self.dt_timestep = carma.dt * carma.output_gap
        self.path = path


        self.gases = {}
        for i in range(1, len(self.gas_names)):
            self.gases[self.gas_names[i]] = self.gas_abund[:, i, :]

        self.clouds = {}
        for i in range(len(self.group_names)):
            self.clouds[self.group_names[i]] = {
                "numden": self.numden[:, i, :, :],
                "r": self.r[:, i] * MICRON_TO_CM,
                "r_mass": self.rmass[:, i]
            }
        if self.carma.is_2d: #TODO: pre nanmean() these.  Maybe have a window for a set number of timesteps
            _, counts = np.unique(step, return_counts=True)
            self.gases_2d = {}

            for i in range(1, len(self.gas_names)):
                self.gases_2d[self.gas_names[i]] = (np.zeros((NZ, 
                                                            carma.NLONGITUDE,
                                                            np.max(counts))) 
                                                  * np.nan)

                index = np.zeros(carma.NLONGITUDE, dtype=int)
                for it in range(n_tstep):
                    self.gases_2d[self.gas_names[i]][
                        :, step[it], index[step[it]]] =  gas_abund[:, i, it] 
                    index[step[it]] += 1 
                self.gases_2d[self.gas_names[i]] = np.nanmean(
                                        self.gases_2d[self.gas_names[i]], 
                                        axis=2)

            for i in range(len(self.group_names)):
                self.clouds[self.group_names[i]]["numden_2d"] = np.zeros(
                                        (NZ, 
                                        NBIN, 
                                        carma.NLONGITUDE, 
                                        np.max(counts))) * np.nan

                index = np.zeros(carma.NLONGITUDE, dtype=int)
                for it in range(n_tstep):
                    self.clouds[self.group_names[i]]["numden_2d"][
                      :, :, step[it], index[step[it]]] = (
                          self.numden[:, i, :, it])
                    
                    index[step[it]] += 1 
                
                self.clouds[self.group_names[i]]["numden_2d"] = np.nanmean(
                    self.clouds[self.group_names[i]]["numden_2d"], axis=3)
                

            def longitude_map(arr):
                index = np.zeros(carma.NLONGITUDE, dtype=int)
                temp_array = np.zeros((NZ, 
                                       NBIN, 
                                       carma.NLONGITUDE, 
                                       np.max(counts))) * np.nan
                for it in range(n_tstep):
                    temp_array[:, :, step[it], index[step[it]]] = arr[:, :, it]
                    index[step[it]] += 1 
                
                return np.nanmean(temp_array, axis=3)
            self.longitude_map = longitude_map



        self._read_temperature(path, path_end, NZ, n_tstep)

        if not read_diag: return

        hom_nuc_gain_rates = np.zeros((NZ, NBIN, NELEM, n_tstep))
        het_nuc_gain_rates = np.zeros((NZ, NBIN, NELEM, n_tstep))
        grow_gain_rates    = np.zeros((NZ, NBIN, NELEM, n_tstep))
        evap_gain_rates    = np.zeros((NZ, NBIN, NELEM, n_tstep))

        nuc_loss_rates     = np.zeros((NZ, NBIN, NGROUP, n_tstep))
        grow_loss_rates    = np.zeros((NZ, NBIN, NGROUP, n_tstep))
        evap_loss_rates    = np.zeros((NZ, NBIN, NGROUP, n_tstep))
        coremass_frac      = np.zeros((NZ, NBIN, NGROUP, n_tstep))


        file_path = os.path.join(path, f"rates_{path_end}.txt")
        f = open(file_path)

        rows_per_step = NBIN * NZ * (NELEM + NGROUP)
        ncols_rate = None

        for it in range(n_tstep):
            t_step = f.readline()
            if not t_step:
                break

            lines = [f.readline() for _ in range(rows_per_step)]
            if not lines[-1]:
                break
            if ncols_rate is None:
                ncols_rate = len(lines[0].split())
            arr = np.fromstring(' '.join(lines), sep=' ', dtype=float)
            if arr.size != rows_per_step * ncols_rate:
                break
            arr = arr.reshape(NBIN, NZ, NELEM + NGROUP, ncols_rate)

            elem_rows  = arr[:, :, :NELEM, :]              # (NBIN, NZ, NELEM, 7)
            group_rows = arr[:, :, NELEM:, :]              # (NBIN, NZ, NGROUP, 7)

            # Target axis order is (NZ, NBIN, N{ELEM|GROUP}, NT) — transpose (1,0,2)
            hom_nuc_gain_rates[:, :, :, it] = elem_rows[..., 3].transpose(1, 0, 2)
            het_nuc_gain_rates[:, :, :, it] = elem_rows[..., 4].transpose(1, 0, 2)
            grow_gain_rates[:, :, :, it]    = elem_rows[..., 5].transpose(1, 0, 2)
            evap_gain_rates[:, :, :, it]    = elem_rows[..., 6].transpose(1, 0, 2)

            nuc_loss_rates[:, :, :, it]  = group_rows[..., 3].transpose(1, 0, 2)
            grow_loss_rates[:, :, :, it] = group_rows[..., 4].transpose(1, 0, 2)
            evap_loss_rates[:, :, :, it] = group_rows[..., 5].transpose(1, 0, 2)
            coremass_frac[:, :, :, it]   = group_rows[..., 6].transpose(1, 0, 2)

        # The Fortran diagnostic collectors integrate each production/loss term
        # over one model timestep: production terms are multiplied by dtime and
        # loss terms by (pc * dtime) before being summed over substeps (see
        # newstate.F90). The raw values are therefore per-cm3 number changes per
        # dt step, not rates. Divide by carma.dt to recover true rates in
        # particles cm^-3 s^-1. coremass_frac is a dimensionless fraction and is
        # left untouched.
        hom_nuc_gain_rates /= carma.dt
        het_nuc_gain_rates /= carma.dt
        grow_gain_rates    /= carma.dt
        evap_gain_rates    /= carma.dt
        nuc_loss_rates     /= carma.dt
        grow_loss_rates    /= carma.dt
        evap_loss_rates    /= carma.dt


        for i, key in enumerate(carma.groups.keys()):
            group = carma.groups[key]

            self.clouds[key]["nuc_loss_rate"] = nuc_loss_rates[:, :, i, :]
            self.clouds[key]["grow_loss_rate"] = grow_loss_rates[:, :, i, :]
            self.clouds[key]["evap_loss_rate"] = evap_loss_rates[:, :, i, :]
            self.clouds[key]["coremass_frac"] = coremass_frac[:, :, i, :]

            if group.mantle:
                e = group.mantle.ielem-1

            else:
                e = group.core.ielem-1

            self.clouds[key]["nuc_gain_rate"] = (
                            hom_nuc_gain_rates[:, :, e, :] 
                            + het_nuc_gain_rates[:, :, e, :])
                
            self.clouds[key]["grow_gain_rate"] = grow_gain_rates[:, :, e, :]
            self.clouds[key]["evap_gain_rate"] = evap_gain_rates[:, :, e, :]
        
        f.close()

    def _read_temperature(self, path, path_end, NZ, n_tstep) -> None:
        """Read the evolving column written by a radiatively coupled run.

        Leaves every attribute None when the file is absent, which is the case
        for any run without ``Carma.set_radiation()``.
        """
        self.T_history = None
        self.dTdt      = None
        self.flux_net  = None
        self.F_TOA     = None
        self.convective = None
        self.nzone     = None
        self.conv_resid = None
        self.P_levels_rad = None
        self.rad_interval = None
        self.teff      = None

        file_path = os.path.join(path, f"temperature_{path_end}.txt")
        if not os.path.exists(file_path):
            return

        with open(file_path) as f:
            header = f.readline().split()
            if not header:
                return
            nz_file, _, iskip = (int(x) for x in header[:3])
            self.teff  = float(header[3])

            if nz_file != NZ:
                raise ValueError("temperature output inconsistent with carma run")

            T_history = np.zeros((NZ, n_tstep))
            dTdt      = np.zeros((NZ, n_tstep))
            convective = np.zeros((NZ, n_tstep), dtype=bool)
            flux_net  = np.zeros((NZ + 1, n_tstep))
            F_TOA     = np.zeros(n_tstep)
            nzone     = np.zeros(n_tstep, dtype=int)
            conv_resid = np.zeros(n_tstep)
            nsolve    = np.zeros(n_tstep)

            # A run killed mid-output leaves a partial record; nt counts only
            # the records that read back whole.
            nt = 0
            for it in range(n_tstep):
                scalars = f.readline().split()
                if len(scalars) < 5:
                    break
                _, F_TOA[it], nsolve[it], zones, conv_resid[it] = (
                    float(x) for x in scalars)
                nzone[it] = int(zones)

                lines = [f.readline() for _ in range(2 * NZ + 1)]
                if not lines[-1]:
                    break
                block = np.fromstring(' '.join(lines[:NZ]), sep=' ')
                if block.size != 5 * NZ:
                    break

                levels = np.fromstring(' '.join(lines[NZ:]), sep=' ')
                if levels.size != 3 * (NZ + 1):
                    break

                block = block.reshape(NZ, 5)
                convective[:, it] = block[:, 4] > 0.5
                T_history[:, it] = block[:, 2]
                dTdt[:, it]      = block[:, 3]

                levels = levels.reshape(NZ + 1, 3)
                flux_net[:, it]  = levels[:, 2]
                # Fixed for the run; the levels the convective zones are
                # bounded by, which nothing else on Results carries.
                self.P_levels_rad = levels[:, 1]
                nt += 1

        self.T_history = T_history[:, :nt]
        self.dTdt      = dTdt[:, :nt]
        self.flux_net  = flux_net[:, :nt]
        self.F_TOA     = F_TOA[:nt]
        self.convective = convective[:, :nt]
        self.nzone     = nzone[:nt]
        self.conv_resid = conv_resid[:nt]

        # Steps per full radiative solve. A value of 1 means the adaptive
        # cadence is doing nothing; the ceiling is the run's rad_gap_max.
        with np.errstate(divide='ignore', invalid='ignore'):
            self.rad_interval = np.where(nsolve[:nt] > 0,
                                         iskip / np.maximum(nsolve[:nt], 1),
                                         np.inf)

    def convective_zones(self, it: int = -1) -> list[tuple[float, float]]:
        """The convective zones at one output, bottom-up.

        Returns ``(p_bottom, p_top)`` in barye for each contiguous run of
        convective layers, deepest zone first.  Empty when the column is
        entirely radiative.

        There is deliberately no single "the RCB": a cloud deck can drive a
        detached convective layer aloft over a stable region, and reducing that
        to one boundary pressure would discard it.

        Parameters
        ----------
        it : int, optional
            Which output to report on, by default the last.

        Raises
        ------
        ValueError
            If the run was not radiatively coupled.
        """
        if self.convective is None:
            raise ValueError("the run was not radiatively coupled")

        mask = self.convective[:, it]

        # Edges of each run of True, found on the mask padded with False so a
        # zone touching the bottom or top of the column is not missed.
        edges = np.flatnonzero(np.diff(np.concatenate(([False], mask, [False]))))
        starts, ends = edges[::2], edges[1::2]

        pl = self.P_levels_rad

        return [(pl[s], pl[e]) for s, e in zip(starts, ends)]

    def final_T_levels(self) -> np.ndarray:
        """The level temperatures the run ended on [K] (NZ+1).

        For an uncoupled run this is ``carma.T_levels`` unchanged.  For a
        radiatively coupled one the input profile is only the starting guess,
        so the levels are rebuilt from the last output the same way
        ``carma_rce`` does -- linear in ``ln p`` with the endpoints
        extrapolated.
        """
        carma = self.carma

        if self.T_history is None:
            return carma.T_levels

        if carma.is_2d:
            raise NotImplementedError(
                "radiative coupling is 1-D only; T_levels for a 2-D run has a "
                "longitude axis the coupled column does not")

        t  = self.T_history[:, -1]
        lp  = np.log(carma.P_centers)
        lpl = np.log(carma.P_levels)

        tl = np.zeros(carma.NZ + 1)

        w = (lpl[1:-1] - lp[1:]) / (lp[:-1] - lp[1:])
        tl[1:-1] = t[1:] + w * (t[:-1] - t[1:])

        w = (lpl[0] - lp[0]) / (lp[0] - lp[1])
        tl[0] = t[0] + w * (t[0] - t[1])

        w = (lpl[-1] - lp[-1]) / (lp[-1] - lp[-2])
        tl[-1] = t[-1] + w * (t[-1] - t[-2])

        return tl

    def final_z_levels(self) -> np.ndarray:
        """The altitude levels the run ended on [cm] (NZ+1).

        The Fortran regrids the column hydrostatically whenever ``t(:)`` moves,
        so a coupled run's layer thicknesses are not the ones ``calculate_z``
        produced at setup.  On an ``I_LOGP`` grid the altitudes are a function
        of pressure alone and do not move, matching ``rce_regrid_z``.
        """
        carma = self.carma

        if self.T_history is None or carma.igridv == I_LOGP:
            return carma.z_levels

        tl = self.final_T_levels()
        H = k_B * tl / (carma.wt_mol * PROTON_MASS * carma.surface_grav)

        z_levels = np.zeros(carma.NZ + 1)
        z_levels[1:] = np.cumsum(
            H[1:] * np.log(carma.P_levels[:-1] / carma.P_levels[1:]))

        return z_levels

    def plot_temperature_evolution(self, nprofile=8, **kwargs):
        """Plots how the P-T profile evolved over a radiatively coupled run.

        Draws a sample of profiles from the first output to the last with the
        final convective zones shaded, the emergent flux relative to
        ``sigma * Teff**4`` -- the quantity equilibrium drives to 1 -- and where
        the convective zones sat over the run, which is how a cloud deck
        reshaping the profile shows itself.

        Parameters
        ----------
        nprofile : int, optional
            Number of profiles to draw between the first and last output,
            by default 8

        Raises
        ------
        RuntimeError
            If the run was not radiatively coupled.
        """
        if self.T_history is None:
            raise RuntimeError(
                "No temperature output; this run was not configured with "
                "Carma.set_radiation()")

        plt.close()
        fig, (ax, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

        nt = self.T_history.shape[1]
        picks = np.unique(np.linspace(0, nt - 1, min(nprofile, nt)).astype(int))
        colors = plt.cm.viridis(np.linspace(0, 1, len(picks)))

        for c, it in zip(colors, picks):
            ax.plot(self.T_history[:, it], self.P * 1e-6, color=c,
                    label=f"{self.ts[it]:.3g} s", **kwargs)

        for i, (p_bot, p_top) in enumerate(self.convective_zones()):
            ax.axhspan(p_bot * 1e-6, p_top * 1e-6, color="k", alpha=0.08,
                       lw=0, label="convective" if i == 0 else None)

        ax.set_yscale("log")
        ax.invert_yaxis()
        ax.set_xlabel("Temperature [K]")
        ax.set_ylabel("Pressure [bar]")
        ax.legend(fontsize="small")

        sigma_teff4 = 5.670374419e-8 * self.teff ** 4          # W/m^2
        ax2.plot(self.ts, self.F_TOA / sigma_teff4, color=petroff10[0])
        ax2.axhline(1, color="k", ls=":", lw=1)
        ax2.set_xlabel("Time [s]")
        ax2.set_ylabel(r"$F_{\rm TOA} / \sigma T_{\rm eff}^4$")

        ax3.pcolormesh(self.ts[:nt], self.P * 1e-6,
                       self.convective.astype(float),
                       cmap="Greys", vmin=0, vmax=1.4, shading="nearest")
        ax3.set_yscale("log")
        ax3.invert_yaxis()
        ax3.set_xlabel("Time [s]")
        ax3.set_ylabel("Pressure [bar]")
        ax3.set_title("convective zones", fontsize="small")

        fig.tight_layout()
        return fig, (ax, ax2, ax3)

    def plot_toa_gas(self, skip_gases = [0], burn_in = 20, **kwargs):
        """Plots the gas abundances at the top of the atmosphere.  Useful for 
        determining whether or not the simulation has converged. ``**kwargs``
        are passed to pyplot

        Parameters
        ----------
        skip_gases : list, optional
            A list of gases, by index, to skip plotting, by default [0] (H2O)
        burn_in : int, optional
            The number of timesteps to exclude from the start of the simulation,
            by default 20
        """
        plt.close()
        fig, ax = plt.subplots()
        ax.set_prop_cycle(mpl.cycler(color=petroff10))
        j = 0
        for i, gas in enumerate(list(self.gas_names)): #TODO get this from header file
            if i not in skip_gases:
                xs = np.arange(burn_in, 
                               len(self.gas_abund[-1, i, :])) * self.dt_timestep
                
                min_abund = np.min(self.gas_abund[-1, i, burn_in:])
                max_abund = np.max(self.gas_abund[-1, i, burn_in:])
                d_abund = max_abund - min_abund + 1e-100
                
                ax.plot(xs, 
                        ((self.gas_abund[-1, i, burn_in:]-min_abund)/d_abund  
                         + j*1.2), 
                         label=f"{gas} (±{d_abund/max_abund:.2e})", 
                         **kwargs)
                

                j -= 1
        plt.xlabel("Time [s]")
        plt.ylabel("Relative gas bundance (offset)")
        plt.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout()


    def plot_numdens(self, nlevels=11, min_order = -10, **kwargs) -> None:
        """Interactively plots the number densities of the cloud species. 
        Densities are normalized to the peak number density in the plot.  
        Interactivly change the cloud species which is plotted and the timestep 
        using the sliders

        Note
        ----
        If in a notebook, requires ``%matplotlib ipympl`` to have been invoked

        Parameters
        ----------
        nlevels : int, optional
            The number of contours to plot, by default 11
        min_order : int, optional
            The number of orders of magnitude below the peak density to show,
            by default -10

        """
        plt.close()

        # The parametrized function to be plotted
        def f(it, ig):
            return np.log10(self.numden[:, ig, :, it] * self.rmass[:, ig]+1e-99)

        t = np.linspace(0, 1, 1000)

        # Define initial parameters"
        t_init = len(self.ts)-1

        # Create the figure and the line that we will manipulate
        fig, ax = plt.subplots()
        ax.set_prop_cycle(mpl.cycler(color=petroff10))

        plt.xlabel("Log₁₀ Radius [μm]")
        plt.ylabel("Log₁₀ Pressure (barye)")

        r = self.r
        P = self.P

        cmap = ax.contourf(np.log10(r[:, 0]),
                           np.log10(P),
                            f(0, 0)- np.max(f(0, 0)),
                           levels=np.linspace(min_order, 0, nlevels),
                           extend="min",
                           **kwargs)
        plt.gca().invert_yaxis()
        # ax.set_xlabel('Time [s]')
        cbar = plt.colorbar(cmap, label="Normalized Log Mass Density")
        # ax.title(group_names[0])

        # plt.gca().axesPatch.set_alpha(0.0)

        text = plt.text((np.log10(r[0, 0])+np.log10(r[-1, 0]))/2,
                        np.log10(P[-1]) - (np.log10(P[0]) -np.log10(P[-1]))*.05,
                        self.group_names[0],
                        ha="center")

        # adjust the main plot to make room for the sliders
        fig.subplots_adjust(left=0.25, bottom=0.25)

        # Make a horizontal slider to control the frequency.
        ax_t = fig.add_axes([0.25, 0.1, 0.65, 0.03])
        t_slider = Slider(
            ax=ax_t,
            label='Time [s]',
            valmin=0,
            valmax=(self.numden.shape[3]-1)*self.dt_timestep,
            valinit=0,
            valstep=self.dt_timestep
        )



        # Make a vertically oriented slider to control the amplitude
        ax_group = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
        group_slider = Slider(
            ax=ax_group,
            label='Group',
            valmin=0,
            valmax=(self.numden.shape[1]-1),
            valinit=0,
            valstep=1,
            orientation="vertical",
            # valfmt=f"{group_names[%d]}"
        )

        # # The function to be called anytime a slider's value changes
        def update(val):
            t_step = int(t_slider.val/self.dt_timestep)
            vals = f(t_step, group_slider.val)
            maxv = np.ceil(np.max(vals))
            # maxv = -10

            cmap = ax.contourf(np.log10(r[:, group_slider.val]),
                               np.log10(P),
                               vals-maxv,
                               levels=np.linspace(min_order, 0, nlevels),
                               extend="min",
                               **kwargs)
            # ax.title(group_names[group_slider.val])
            text.set_text(self.group_names[group_slider.val])
            # cbar.clim(maxv-6, maxv)
            cbar.update_normal(cmap)

            fig.canvas.draw_idle()


        # # register the update function with each slider
        t_slider.on_changed(update)
        group_slider.on_changed(update)


    def _plot_abundance_profile(self, skip_gases = [0], **kwargs):
        plt.close()
        lines = []
        fig, ax = plt.subplots()
        ax.set_prop_cycle(mpl.cycler(color=petroff10))

        for i, gas in enumerate(self.gas_names):
            if i not in skip_gases:
                zs = self.Z
                lines.append(ax.plot((self.gas_abund[:, i, -1]/1e6 * self.P), self.P, label=gas)[0])
        plt.yscale("log")
        plt.xscale("log")
        plt.xlabel("Partial Pressure [barye]")
        plt.ylabel("Pressure [barye]")
        plt.gca().invert_yaxis()
        plt.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout()

        fig.subplots_adjust(left=0.25, bottom=0.25)

        axslider = fig.add_axes([0.25, 0.1, 0.65, 0.03])
        slider = Slider(
            ax=axslider,
            label='Time [s]',
            valmin=0,
            valmax=(self.numden.shape[3]-1)*self.dt_timestep,
            valinit=0,
            valstep=self.dt_timestep
        )


        def update2(val):
            t_step = int(slider.val/self.dt_timestep)
            for i, gas in enumerate(self.gas_names):
                if i != 0:
                    lines[i-1].set_xdata((self.gas_abund[:, i, t_step]/1e6 * self.P))
                    fig.canvas.draw_idle()


        slider.on_changed(update2)
        plt.show()
        
    def _plot_saturation(self, skip_gases=[0]):
        plt.close()
        lines = []
        # plt.style.use("petroff10")
        fig, ax = plt.subplots()
        ax.set_prop_cycle(mpl.cycler(color=petroff10))

        for i, gas in enumerate(self.gas_names):
            if i not in skip_gases:
                zs = self.Z
                lines.append(ax.plot(((self.gas_abund[:, i, -1]/1e6 * self.P)
                                     /self.sat_vp[:, i, -1]),
                                     self.P,
                                     label=gas)[0])
        plt.yscale("log")
        plt.xscale("log")
        plt.gca().invert_yaxis()
        plt.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout()
        fig.subplots_adjust(left=0.25, bottom=0.25)

        axfreq = fig.add_axes([0.25, 0.1, 0.65, 0.03])
        slider = Slider(
            ax=axfreq,
            label='Time Step',
            valmin=0,
            valmax=(self.numden.shape[3]-1)*self.dt_timestep,
            valinit=0,
            valstep=self.dt_timestep
        )


        def update2(val):
            t_step = int(slider.val / self.dt_timestep)
            for i, gas in enumerate(self.gas_names):
                if i != 0:
                    lines[i-1].set_xdata((self.gas_abund[:, i, slider.val]/1e6
                                        * self.P /self.sat_vp[:, i, slider.val]))

                    fig.canvas.draw_idle()


        slider.on_changed(update2)

        plt.show()

    def plot_condensation_curves(self, ax=None, t_step=-1, skip_gases=None, **kwargs):
        """Plots condensation curves for each gas species in the CARMA model
        overlaid on the model's P-T profile.

        For each gas the condensation curve is the locus of (T, P) points at
        which the gas partial pressure equals its saturation vapour pressure
        (i.e. saturation ratio = 1).  The curve is computed analytically from
        the vapour-pressure coefficients stored in ``constants.gas_dict`` using
        the deep-atmosphere gas mixing ratio from the simulation.  Species whose
        vapour-pressure formula does not have a temperature coefficient (e.g. H2O,
        which uses a custom Fortran routine) are plotted using the saturation
        vapour pressures stored directly in ``results.sat_vp``.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on.  A new figure is created when *None*.
        t_step : int, optional
            Time-step index used to read gas abundances, by default ``-1``
            (last recorded step).
        skip_gases : list[int], optional
            List of gas indices to exclude.  Defaults to ``[]`` (plot all
            species).
        **kwargs
            Extra keyword arguments forwarded to ``ax.plot`` for the P-T
            profile line only.

        Returns
        -------
        fig : matplotlib.figure.Figure
        ax : matplotlib.axes.Axes
        """
        if skip_gases is None:
            skip_gases = []

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        ax.set_prop_cycle(mpl.cycler(color=petroff10))

        P = self.P      # [barye], shape (NZ,)
        # In a radiatively coupled run self.T is only the starting profile, so
        # the condensation curves have to be read against the evolved one.
        T = self.T if self.T_history is None else self.T_history[:, t_step]
        if T.ndim > 1:
            T = T[:, 0]  # use first longitude for 2D runs

        log_met = getattr(self.carma, "log_metallicity", 0.0)

        # P-T profile plotted in black so the colour cycle is free for species
        ax.plot(T, P / 1e6, color="black", lw=2, label="P-T profile", **kwargs)

        # Pressure grid for computing smooth condensation curves
        P_curve = np.logspace(np.log10(P.min()), np.log10(P.max()), 300)

        for i, gas_name in enumerate(self.gas_names):
            if i in skip_gases:
                continue

            # Deep abundance: max over altitude gives uninhibited mixing ratio
            x_gas = np.max(self.gas_abund[:, i, t_step]) / 1e6

            if x_gas <= 0:
                continue

            entry = gas_dict.get(gas_name, {})
            tcoeff = entry.get("vp_tcoeff", 0)

            if tcoeff == 0:
                # No analytic inversion available (e.g. H2O Murphy 2005).
                # Use stored sat_vp at each model level: the condensation
                # pressure for a given T is sat_vp / x_gas.
                order = np.argsort(T)
                P_cond = self.sat_vp[:, i, t_step] / x_gas  # barye
                ax.plot(T[order], P_cond[order] / 1e6, label=gas_name)
            else:
                offset    = entry.get("vp_offset",    0.0)
                metcoeff  = entry.get("vp_metcoeff",  0.0)
                logpcoeff = entry.get("vp_logpcoeff", 0.0)

                # Condensation condition: sat_vp(T, P) = x_gas * P
                #   1e6 * 10^(offset - tcoeff/T - metcoeff*log_met
                #              - logpcoeff*log10(P*1e-6)) = x_gas * P
                # Solving for T:
                #   T = tcoeff / (6 + offset - metcoeff*log_met
                #                 - logpcoeff*log10(P*1e-6)
                #                 - log10(x_gas * P))
                denom = (
                    6.0
                    + offset
                    - metcoeff  * log_met
                    - logpcoeff * np.log10(P_curve * 1e-6)
                    - np.log10(x_gas * P_curve)
                )
                with np.errstate(invalid="ignore", divide="ignore"):
                    T_cond = np.where(denom > 0, tcoeff / denom, np.nan)

                ax.plot(T_cond, P_curve / 1e6, label=gas_name)

        ax.set_yscale("log")
        ax.invert_yaxis()
        ax.set_xlabel("Temperature [K]")
        ax.set_ylabel("Pressure [bar]")
        ax.legend(bbox_to_anchor=(1.0, 1.0), loc="upper left")
        fig.tight_layout()

        return fig, ax

    def gen_picaso_atm_file(self, file_path: str = None,
                            longitude: int = None,
                            suppress_output: bool = False) -> None:
        """Generates an atmosphere file for use in picaso.
        Picaso is a python package which can calculate spectra and is available
        at https://github.com/natashabatalha/picaso.

        Parameters
        ----------
        file_path : str, optional
            The path to save the file to.  If not provided, creates a file
            located in the directory storing the carma simulation.
        longitude : int
            Longitude index to use.  Required for 2-D runs; must not be set
            for 1-D runs.
        suppress_output : bool, optional
            If True, suppress the "Wrote file" confirmation message.

        """
        self.carma._citation["picaso"] = True

        if self.carma.is_2d and longitude is None:
            raise ValueError("longitude is required for 2-D runs")
        if not self.carma.is_2d and longitude is not None:
            raise ValueError("longitude cannot be set for 1-D runs")

        if not file_path: file_path = os.path.join(self.path, 'fastchem.atm')


        # species = ['H2O1', 'C1H4', 'C1O1', 'C1O2', 'Na', 'K', 'H2S1', 'C1H1N1_1', 'O2S1', 'H', 'H2', 'He', 'H1-', 'H1+', 'e-']
        species = ['H1O1','H2','H2O1','H','O','C1H1','C','C1H2','C1H3','C1H4',
                   'C1O1','C1O2','O2','N','H1N1','C1N1','C1H1N1_1','N1O1',
                   'H2N1','N2','H3N1','H1S1','H2O4S1','H2S1','S','S2','O1S1',
                   'C1S1','C1O1S1','C1S2','N1S1','O2S1','S4','S8','S3','O1S2',
                   'O1Ti1','Ti','O2Ti1','H1Ti1','C2Ti1','N1Ti1','O1V1','V','He',
                   'Na','K', 'H1-', 'H1+', 'e-', 'Fe1H1']
        # species_labels = ['H2O', 'CH4', 'CO', 'CO2', 'Na', 'K', 'H2S', 'HCN', 'SO2', 'H', 'H2', 'He', 'H-', 'H+', 'e-']
        species_labels = ['OH','H2','H2O','H','O','CH','C','CH2','CH3','CH4',
                  'CO','CO2','O2',
                  'N','NH','CN','HCN','NO','NH2','N2','NH3',
                  'SH','H2SO4','H2S','S','S2','SO','CS','COS','CS2','NS','SO2','S4','S8',
                  'S3','S2O',
                  'TiO','Ti','TiO2','TiH','TiC2','TiN','VO','V','He','Na','K','H1-', 'H1+', 'e-', 'FeH']
        # species_labels = ['H2', 'He', 'CH4', 'CO', 'CO2', 'NH3', 'N2', 'H2O', 'TiO', 'VO', 'FeH', 'H2S']
        # species = ["H2", "He", "C1H4", "C1O1", "C1O2", "H3N1", "N2", "H2O1", "O1Ti1", "O1V1", "Fe1H1", "H2S1"]

        # species_dict = dict(zip(species, species_labels)) #used to convert FROM HILL notation to readable

        metallicity = 10**self.carma.log_metallicity

        T = (self.carma.T_levels[:, longitude] if self.carma.is_2d
             else self.final_T_levels())

        data = get_fastchem_abundances(T, self.carma.P_levels, species, metallicity)
        data = np.vstack((self.carma.P_levels/BAR_TO_BARYE, T, data))

        header = "pressure\ttemperature"
        for s in species_labels:
            header += '\t'
            header += s

        np.savetxt(file_path, np.transpose(data), header=header, fmt="%.18e", comments="")
        if not suppress_output: print(f"Wrote file: {file_path}")

    def gen_picaso_cloud_file(self,
                              wavelengths: np.ndarray,
                              file_path: str = None,
                              mie_table_path: str = None,
                              skip_groups=[],
                              n_avg: int = 20,
                              longitude: int = None,
                              suppress_output: bool = False):
        """Generate cloud opacity tables for use in picaso. Picaso is a python
        package which can calculate spectra and is available at
        https://github.com/natashabatalha/picaso.

        Parameters
        ----------
        wavelengths : np.ndarray
            The wavelengths at which to calculate to opacities [cm]
        file_path : str, optional
            The path to save the file to.  If not provided, creates a file
            located in the directory storing the carma simulation.
        mie_table_path : str
            The directory in which the mie tables of the species are located.
            It is assumed that each cloud species used has a .dat file in that
            directory named with the name of the cloud species (ie 'Mg2SiO4 on
            TiO2.dat'), the first row of the table is a header row, and the
            columns are radius[cm], wavelength[cm], extinction efficiency
            (Q_ext), scattering efficiency (Q_sca), and asymmetry factor (g).
            The  columns must be in that order and separated by whitespace.  If
            no path is provided, the carma default tables will be used but these
            might be less accurate, even for the same species, as they might not
            be calculated at the same particle size and wavelengths.
        skip_groups : list, optional
            The indices of any cloud groups to exclude from the opacity
            calculation, by default None
        n_avg : int, optional
            The number of timesteps to average the cloud size distribution over
            while generating the cloud file.  Defaults to 20.
        longitude : int
            Longitude index to use.  Required for 2-D runs; must not be set
            for 1-D runs.
        suppress_output : bool, optional
            If True, suppress the "Wrote file" confirmation message.
        """

        self.carma._citation["picaso"] = True

        if self.carma.is_2d and longitude is None:
            raise ValueError("longitude is required for 2-D runs")
        if not self.carma.is_2d and longitude is not None:
            raise ValueError("longitude cannot be set for 1-D runs")

        if not file_path: file_path = os.path.join(self.path, 'clouds.atm')


        carma = self.carma
        beta_exts = []
        beta_scas = []
        g_avgs     = []

        P = carma.results.P
        idx = np.argmin(np.abs(P - 1e4))

        for i in range(len(carma.groups)):
            if i in skip_groups: continue
            beta_ext, beta_sca, g_avg = _get_cloud_opacities(carma,
                                                             i,
                                                             wavelengths,
                                                             mie_table_path,
                                                             n_avg=n_avg,
                                                             longitude=longitude)

            beta_exts.append(beta_ext)
            beta_scas.append(beta_sca)
            g_avgs.append(g_avg)

        beta_ext = np.sum(np.array(beta_exts), axis=0)
        beta_sca = np.sum(np.array(beta_scas), axis=0)
        g_avg = np.sum(np.array(g_avgs) * np.array(beta_scas), axis=0) / (beta_sca + 1e-100)
        g_avg = np.where(beta_sca==0, 0, g_avg)

        ssas = beta_sca/(beta_ext + 1e-100)
        ssas = np.where(beta_ext==0, 0, ssas)

        z_levels = self.final_z_levels()
        dz = np.abs(z_levels[1:] - z_levels[:-1])

        d_tau = beta_ext * dz[:, np.newaxis]

        with open(file_path, "w+") as f:
            f.write("nlayer\tnwave\tpressure\twavenumber\tw0\tg0\topd\n")
            for iz in range(carma.NZ):
                for iwave in range(len(wavelengths)):
                    f.write(f"{iz+1}\t"
                    +f"{iwave+1}\t"
                    +f"{carma.P_centers[iz]/BAR_TO_BARYE}\t"
                    +f"{1/wavelengths[iwave]}\t"
                    +f"{ssas[iz, iwave]}\t"
                    +f"{g_avg[iz, iwave]}\t"
                    +f"{d_tau[iz, iwave]}\n")

        if not suppress_output: print(f"Wrote file: {file_path}")





def _get_cloud_opacities(carma,
                         i,
                         wavelengths,
                         mie_table_path = None,
                         min_columnden = 1e-25,
                         n_avg = 20,
                         longitude = None):

    name = carma.results.group_names[i]

    
    if name.split()[0] == "Pure":
        species = name.split()[-1]
    else:
        species = name.split()[0]

    if mie_table_path:
        data = np.genfromtxt(mie_table_path, delimiter='\t', names=True)
    else:
        try:
            data = np.genfromtxt(
                os.path.join(_SRC, "mie_tables", f'{species}_user.dat'),
                 delimiter='\t', 
                 names=True)
        except: #TODO: handle this exception properly
            data = np.genfromtxt(
                os.path.join(_SRC, "mie_tables", f'{name}.dat'),
                delimiter='\t', 
                names=True)

         

    r = data["rcm"]
    λ = data["λcm"]

    r_unique = np.unique(r)
    λ_unique = np.unique(λ)

    r_idx = np.searchsorted(r_unique, r)
    λ_idx = np.searchsorted(λ_unique, λ)

    Qext = np.empty((len(r_unique), len(λ_unique)))
    Qext[r_idx, λ_idx] = data["Q_ext"]
    Qext_interp = RectBivariateSpline(r_unique, λ_unique, Qext)

    Qsca = np.empty((len(r_unique), len(λ_unique)))
    Qsca[r_idx, λ_idx] = data["Q_sca"]
    Qsca_interp = RectBivariateSpline(r_unique, λ_unique, Qsca)

    g = np.empty((len(r_unique), len(λ_unique)))
    g[r_idx, λ_idx] = data["g"]
    g_interp = RectBivariateSpline(r_unique, λ_unique, g)


    if longitude is not None:
        # numden_2d is already the mean over all timesteps at each longitude
        numdens = carma.results.clouds[name]["numden_2d"][:, :, longitude]
    else:
        numdens = np.mean(carma.results.numden[:, i, :, -n_avg:], axis=2)
    columndens = np.sum(numdens, axis=0)

    weighted_qext = np.zeros((carma.NZ, carma.NBIN, len(wavelengths)))
    weighted_qsca = np.zeros((carma.NZ, carma.NBIN, len(wavelengths)))
    weighted_g    = np.zeros((carma.NZ, carma.NBIN, len(wavelengths)))


    for ibin in range(carma.NBIN):
        if columndens[ibin] < min_columnden: continue
        for ilambda in range(len(wavelengths)):
            # print(2 * np.pi * carma.results.r[ibin, i]/ wavelengths[ilambda])
            #

            r = carma.results.r[ibin, i] * 1e-4 # micron to cm
            wavelength = wavelengths[ilambda]

            weight_term =  np.pi * r**2 * numdens[:, ibin]



            weighted_qext[:, ibin, ilambda] = weight_term * Qext_interp.ev(r, wavelength)
            weighted_qsca[:, ibin, ilambda] = weight_term * Qsca_interp.ev(r, wavelength)
            weighted_g[:, ibin, ilambda] = weight_term * Qsca_interp.ev(r, wavelength) * g_interp.ev(r, wavelength)

    
    beta_ext = np.sum(weighted_qext, axis=1)
    beta_sca = np.sum(weighted_qsca, axis=1)
    g_avg = np.sum(weighted_g, axis=1) / (beta_sca + 1e-100)
    g_avg = np.where(beta_sca==0, 0, g_avg)
    
    # for ilambda in range(len(wavelengths)):
    #     m = n_interp(wavelengths[ilambda]) + 1j*k_interp(wavelengths[ilambda])
    #     mie_res = ps.Mie_SD(m, wavelengths[ilambda], carma.results.r[:, i]*2e7, numdens)


    return beta_ext, beta_sca, g_avg


        


