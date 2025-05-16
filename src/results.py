import numpy as np
from .constants import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button, Slider
import matplotlib as mpl
from itertools import cycle
from scipy.signal import periodogram

# petroff 10 color cycle
petroff10 = ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"]


class Results:
    
    def __init__(self, carma):
        path = carma.name
        path_end = path.split("/")[-1]
        file_path = path+f"/bd {path_end}.txt"
        
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
        
        r = np.zeros(NBIN)
        rmass = np.zeros((NBIN, NGROUP))
        
        for i in range(NGROUP):
            for j in range(NBIN):
                _, _, r[j], rmass[j, i], _, _, _ = np.array(f.readline().split(), dtype=float)
            
        
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
        self.group_names = list(carma.groups.keys())
        self.gas_names = list(carma.gasses.keys())
                        
                        
    def plot_toa_gas(self, skip_gasses = [0], burn_in = 20, **kwargs):
        plt.close()
        fig, ax = plt.subplots()
        ax.set_prop_cycle(mpl.cycler(color=petroff10))
        j = 0
        for i, gas in enumerate(list(self.gas_names)): #TODO get this from header file
            if i not in skip_gasses:
                xs = np.arange(burn_in, len(self.gas_abund[-1, i, :]))
                ax.plot(xs, (self.gas_abund[-1, i, burn_in:]/np.max(self.gas_abund[-1, i, burn_in:])  + j), label=gas, **kwargs)
                j -= 1
        plt.xlabel("Time Step")
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

        cmap = ax.contourf(np.log10(r),
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

        text = plt.text((np.log10(r[0])+np.log10(r[-1]))/2,
                        np.log10(P[-1]) -  (np.log10(P[0]) -np.log10(P[-1]))*.05,  
                        self.group_names[0], 
                        ha="center")

        # adjust the main plot to make room for the sliders
        fig.subplots_adjust(left=0.25, bottom=0.25)

        # Make a horizontal slider to control the frequency.
        ax_t = fig.add_axes([0.25, 0.1, 0.65, 0.03])
        t_slider = Slider(
            ax=ax_t,
            label='Time Step',
            valmin=0,
            valmax=(self.numden.shape[3]-1),
            valinit=0,
            valstep=1
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
            vals = f(t_slider.val, group_slider.val)
            maxv = np.ceil(np.max(vals))
            # maxv = -10
            
            cmap = ax.contourf(np.log10(r),
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
            label='Time Step',
            valmin=0,
            valmax=(self.numden.shape[3]-1),
            valinit=0,
            valstep=1
        )


        def update2(val):
            for i, gas in enumerate(self.gas_names):
                if i != 0:
                    lines[i-1].set_xdata((self.gas_abund[:, i, slider.val]/1e6 * self.P))
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
            valmax=(self.numden.shape[3]-1),
            valinit=0,
            valstep=1
        )


        def update2(val):
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