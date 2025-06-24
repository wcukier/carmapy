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
import PyMieScatt as ps
from tqdm import tqdm
from scipy.interpolate import RectBivariateSpline

# petroff 10 color cycle
petroff10 = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]
SRC = os.path.dirname(__file__)


class Results:
    
    def __init__(self, carma):
        path = carma.name
        path_end = os.path.basename(path)
        file_path = os.path.join(path, f"bd_{path_end}.txt")
        
        f = open(file_path)
        NZ, NGROUP, NELEM, NBIN, NGAS, nstep, iskip = np.array(f.readline().split(),
                                                            dtype=int)

        
        if ((NZ != carma.NZ) + 
            (NGROUP != len(carma.groups))+
            (NELEM != len(carma.elems))+
            (NBIN != carma.NBIN) +
            (NGAS != len(carma.gasses))+
            (nstep - 1 != carma.n_tstep)+
            (iskip != carma.output_gap)
        ):
            raise ValueError(f"Output file inconsistent with carma run")
        
        r = np.zeros((NBIN, NGROUP))
        rmass = np.zeros((NBIN, NGROUP))
        
        for i in range(NGROUP):
            for j in range(NBIN):
                _, _, r[j, i], rmass[j, i], _, _, _ = np.array(f.readline().split(), dtype=float)
            
        
        kzz = np.zeros(NZ)
        P = np.zeros(NZ)
        T = np.zeros(NZ)
        Z = np.zeros(NZ)

        for i in range(NZ):
            _, Z[i], _, P[i], T[i], kzz[i] = np.array(f.readline().split(), 
                                                    dtype=float)

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
        
        for it in range(NT):
            t_step = f.readline()
            if t_step:
                ts[it] = t_step
                for ibin in range(NBIN):
                    for iz in range(NZ):
                        line = np.array(f.readline().split(), dtype=float)
                        for ielem in range(NELEM):
                            numden[iz, ielem, ibin, it] = line[ielem+2]
                        for igas in range(NGAS):
                            gas_abund[iz, igas, it] = line[NELEM + 2+ 2*igas]
                            sat_vp[iz, igas, it] = line[NELEM + 3 + 2*igas]
            else:
                break   
            
        numden_groups = np.zeros((NZ, NGROUP, NBIN, NT))
        
        
        for i, key in enumerate(carma.groups.keys()):
            group = carma.groups[key]
            if group.mantle:
                numden_groups[:, group.igroup-1, :, :] = numden[:, group.mantle.ielem-1, :, :]
            else:
                numden_groups[:, group.igroup-1, :, :] = numden[:, group.core.ielem-1, :, :]
            if np.any( numden_groups[:, group.igroup-1, :, :] < 0):
                print(key, group.mantle, np.max( numden[:, group.core.ielem-1, :, :]))
                raise
        
        
        self.carma = carma
        self.rmass = rmass
        self.r = r
        self.numden = numden_groups[:,:,:,:it]
        self.gas_abund = gas_abund[:,:,:it]
        self.sat_vp = sat_vp[:,:,:it]
        self.ts = ts[:it]
        self.P = P
        self.Z = Z
        self.T = T
        self.group_names = list(carma.groups.keys())
        self.gas_names = list(carma.gasses.keys())
        self.dt_timestep = carma.dt * carma.output_gap 
        self.path = path
                        
                        
    def plot_toa_gas(self, skip_gasses = [0], burn_in = 20, **kwargs):
        plt.close()
        fig, ax = plt.subplots()
        ax.set_prop_cycle(mpl.cycler(color=petroff10))
        j = 0
        for i, gas in enumerate(list(self.gas_names)): #TODO get this from header file
            if i not in skip_gasses:
                xs = np.arange(burn_in, len(self.gas_abund[-1, i, :]))*self.dt_timestep
                ax.plot(xs, (self.gas_abund[-1, i, burn_in:]/np.max(self.gas_abund[-1, i, burn_in:])  + j), label=gas, **kwargs)
                j -= 1
        plt.xlabel("Time [s]")
        plt.ylabel("Relative gas bundance (offset)")
        plt.legend(bbox_to_anchor=(1, 1))
        plt.tight_layout()


    def plot_numdens(self, nlevels=11, min_order = -10, **kwargs):
        plt.close()

        # The parametrized function to be plotted
        def f(it, ig):
            return np.log10(self.numden[:, ig, :, it] * self.rmass[:, ig]+1e-100)

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
                        np.log10(P[-1]) -  (np.log10(P[0]) -np.log10(P[-1]))*.05,  
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
            
            
    def plot_abundance_profile(self, skip_gasses = [0], **kwargs):
        plt.close()
        lines = []
        fig, ax = plt.subplots()
        ax.set_prop_cycle(mpl.cycler(color=petroff10))

        for i, gas in enumerate(self.gas_names):
            if i not in skip_gasses:
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
        
    def plot_saturation(self, skip_gasses=[0]):
        plt.close()
        lines = []
        # plt.style.use("petroff10")
        fig, ax = plt.subplots()
        ax.set_prop_cycle(mpl.cycler(color=petroff10))

        for i, gas in enumerate(self.gas_names):
            if i not in skip_gasses:
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

    def is_prob_converged(self, burn_in=450, thresh=1e-15):
        gas_abundances = self.gas_abund[-1, :, burn_in:]


        for i in range(len(gas_abundances)):
            periods, strength = periodogram(gas_abundances[i, :] - np.mean(gas_abundances[i, :]), 1)

            if np.max(strength) > thresh:
                long_period = max(10, 1/periods[np.argmax(strength)])
            else:
                long_period = 20

            if (long_period * 2 > len(gas_abundances[0, :])): return False

            long_period = int(3*long_period+1)

            if (long_period * 2 > len(gas_abundances[0, :])): long_period = int(len(gas_abundances[0, :])/2)
        

            mins = np.array([np.min(gas_abundances[i, j:j+long_period]) for j in range(len(gas_abundances[i, :]) - 2*long_period, len(gas_abundances[i, :]) - long_period)])
            maxes = np.array([np.max(gas_abundances[i, j:j+long_period]) for j in range(len(gas_abundances[i, :]) - 2*long_period, len(gas_abundances[i, :]) - long_period)])
            means = np.array([np.mean(gas_abundances[i, j:j+long_period]) for j in range(len(gas_abundances[i, :]) - 2*long_period, len(gas_abundances[i, :]) - long_period)])

            if (np.max(np.abs(mins - np.mean(mins)))/np.mean(mins) > 0.5): return False
            if (np.max(np.abs(means - np.mean(means)))/np.mean(means) > 0.5): return False
            if (np.max(np.abs(maxes - np.mean(maxes)))/np.mean(maxes) > 0.5): return False

            # print(f"{self.gas_names[i]}\t {np.max(np.abs(mins - np.mean(mins)))/np.mean(mins):.3e}\t {np.max(np.abs(maxes - np.mean(maxes)))/np.mean(maxes):.3e}\t {np.max(np.abs(means - np.mean(means)))/np.mean(means):.3e}")
        return True
    
    def gen_picaso_atm_file(self):
        # species = ['H2O1', 'C1H4', 'C1O1', 'C1O2', 'Na', 'K', 'H2S1', 'C1H1N1_1', 'O2S1', 'H', 'H2', 'He', 'H1-', 'H1+', 'e-']
        species = ['H1O1','H2','H2O1','H','O','C1H1','C','C1H2','C1H3','C1H4',
                  'C1O1','C1O2','O2',
                  'N','H1N1','C1N1','C1H1N1_1','N1O1','H2N1','N2','H3N1',
                  'H1S1','H2O4S1','H2S1','S','S2','O1S1','C1S1','C1O1S1','C1S2','N1S1','O2S1','S4','S8',
                  'S3','O1S2',
                  'O1Ti1','Ti','O2Ti1','H1Ti1','C2Ti1','N1Ti1','O1V1','V','He','Na','K', 'H1-', 'H1+', 'e-']
        # species_labels = ['H2O', 'CH4', 'CO', 'CO2', 'Na', 'K', 'H2S', 'HCN', 'SO2', 'H', 'H2', 'He', 'H-', 'H+', 'e-']
        species_labels = ['OH','H2','H2O','H','O','CH','C','CH2','CH3','CH4',
                  'CO','CO2','O2',
                  'N','NH','CN','HCN','NO','NH2','N2','NH3',
                  'SH','H2SO4','H2S','S','S2','SO','CS','COS','CS2','NS','SO2','S4','S8',
                  'S3','S2O',
                  'TiO','Ti','TiO2','TiH','TiC2','TiN','VO','V','He','Na','K','H1-', 'H1+', 'e-']
        species_dict = dict(zip(species, species_labels)) #used to convert FROM HILL notation to readable

        data = get_fastchem_abundances(self.carma.T_levels, self.carma.P_levels, species, 1)
        data = np.vstack((self.carma.P_levels/BAR_TO_BARYE, self.carma.T_levels, data))

        header = "pressure\ttemperature"
        for s in species_labels:
            header += '\t'
            header += s

        np.savetxt(f'{self.path}/fastchem.atm', np.transpose(data), header=header, fmt="%.18e", comments="")

    def gen_picaso_cloud_file(self, wavelengths, skip_groups=[]):
        carma = self.carma
        beta_exts = []
        beta_scas = []
        g_avgs     = []

        P = carma.results.P
        idx = np.argmin(np.abs(P - 1e4))

        for i in range(len(carma.groups)):
            if i in skip_groups: continue
            beta_ext, beta_sca, g_avg = get_cloud_opacities(carma, i, wavelengths)

            beta_exts.append(beta_ext)
            beta_scas.append(beta_sca)
            g_avgs.append(g_avg)

            print(f"{carma.results.group_names[i]}:\t{np.max(beta_ext)}")
        
        beta_ext = np.sum(np.array(beta_exts), axis=0)
        beta_sca = np.sum(np.array(beta_scas), axis=0)
        g_avg = np.sum(np.array(g_avgs) * np.array(beta_scas), axis=0) / beta_sca
        g_avg = np.where(beta_sca==0, 0, g_avg)

        ssas = beta_sca/beta_ext
        ssas = np.where(beta_ext==0, 0, ssas)

        dz = np.abs(carma.z_levels[1:] - carma.z_levels[:-1])

        d_tau = beta_ext * dz[:, np.newaxis]
        
        with open(os.path.join(carma.name, "clouds.atm"), "w+") as f:
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




def _load_cloud_indices(specie_name):
    data = np.genfromtxt(os.path.join(os.path.dirname(__file__), 
                                      "refractive_indices_txt_files",
                                        opacity_files[specie_name]), 
                        comments="#")
    wavelengths = data[:, 0] * 1e-4 # convert to cm
    n_interp = interp1d(wavelengths, data[:, 1])
    k_interp = interp1d(wavelengths, data[:, 2])
    return n_interp, k_interp

def get_cloud_opacities(carma, i, wavelengths, min_columnden = 1e-25):

    name = carma.results.group_names[i]

    
    data = np.genfromtxt(os.path.join(SRC, "mie_tables", f'{name}.dat'), delimiter='\t', names=True) #TODO

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


    numdens = np.mean(carma.results.numden[:, i, :, -20:], axis=2)
    columndens = np.sum(numdens, axis=0)

    weighted_qext = np.zeros((carma.NZ, carma.NBIN, len(wavelengths)))
    weighted_qsca = np.zeros((carma.NZ, carma.NBIN, len(wavelengths)))
    weighted_g    = np.zeros((carma.NZ, carma.NBIN, len(wavelengths)))


    for ibin in range(carma.NBIN):
        if columndens[ibin] < min_columnden: continue
        for ilambda in range(len(wavelengths)):
            # print(2 * np.pi * carma.results.r[ibin, i]/ wavelengths[ilambda])
            #        

            r = carma.results.r[ibin, i] * 1e-4
            wavelength = wavelengths[ilambda]

            weight_term =  np.pi * r**2 * numdens[:, ibin]



            weighted_qext[:, ibin, ilambda] = weight_term * Qext_interp(r, wavelength)
            weighted_qsca[:, ibin, ilambda] = weight_term * Qsca_interp(r, wavelength)
            weighted_g[:, ibin, ilambda] = weight_term * Qsca_interp(r, wavelength) * g_interp(r, wavelength)

    
    beta_ext = np.sum(weighted_qext, axis=1)
    beta_sca = np.sum(weighted_qsca, axis=1)
    g_avg = np.sum(weighted_g, axis=1) / beta_sca
    g_avg = np.where(beta_sca==0, 0, g_avg)
    
    # for ilambda in range(len(wavelengths)):
    #     m = n_interp(wavelengths[ilambda]) + 1j*k_interp(wavelengths[ilambda])
    #     mie_res = ps.Mie_SD(m, wavelengths[ilambda], carma.results.r[:, i]*2e7, numdens)


    return beta_ext, beta_sca, g_avg


        


    