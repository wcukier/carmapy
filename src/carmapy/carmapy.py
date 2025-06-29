from carmapy.constants import *
from carmapy.results import *
import os
import shutil
import f90nml
import numpy as np
import subprocess
import warnings
import shlex
import contextlib
from numpy.typing import ArrayLike

SRC = os.path.dirname(os.path.dirname(__file__))



@contextlib.contextmanager
def _cd(path):
   old_path = os.getcwd()
   os.chdir(path)
   try:
       yield
   finally:
       os.chdir(old_path)




class Carma:
    def __init__(self, name: str, is_2d=False) -> None:
        
        
        self.is_2d:         bool        = is_2d          # true for 2d carma, false for 1d carma
        self.NZ:            int         = 0                 # number of bins in the vertical (z) direction
        self.NBIN:          int         = 80              # number of bins in particle radius
        self.NLONGITUDE:    int         = 64        # number of longitude bins (ignored if is_2d==False)
        self.P_levels:      ArrayLike   = None        # Pressure levels in barye, starting at bottom of atmosphere
        self.P_centers:     ArrayLike   = None       # Pressure centers in barye, starting at bottom of atmosphere
        self.z_levels:      ArrayLike   = None        # Altitude levels in cm, starting at bottom of atmosphere where z=0
        self.z_centers:     ArrayLike   = None       # Altitude centers in cm, starting at bottom of atmosphere where z=0
        self.T_centers:     ArrayLike   = None       # Temperature centers in K, starting at bottom of atmosphere
        self.T_levels:      ArrayLike   = None        # Temperature levels in K, starting at the bottom of the atmosphere
        self.kzz_levels:    ArrayLike   = None      # Eddy diffusion coefficent k_zz in cm^2/s, starting at bottom of atmosphere
        self.idiag:         int         = 0
        self.iappend:       int         = 0


        self.groups:    dict[str, "Group"]      = {}      # dictionary of carma Group objects
        self.growth:    list["Growth"]          = []      # list of carma Growth objects
        self.elems:     dict[str, "Element"]    = {}      # dictionary of carma Element objects
        self.gasses:    dict[str, "Element"]    = {}      # dictionary of carma Gas objects
        self.nucs:      list["Nuc"]             = []      # list of carma Nuc objecs
        self.coags:     list["Coag"]            = []      # dictionary of carma Coag objects
        
        self.name:      str = name            # Name for the carma object, used to define the directory the object is saved in
        self.surface_grav: float = None    # surface gravity of the planet, in cm/s^2
        self.wt_mol: float = None          # mean molecular weight of the atmosphere, in amu
        self.r_planet: float = 6.991e9     # radius of the planet, in cm
        self.velocity_avg: float = -1      # average longitudinal velocity in cm/s, ignored if is_2d==False
        self.restart: bool = False        # if True will restart the carma run from the saved state
        
        self.dt: int = 1000              # carma timestep in seconds
        self.output_gap: int = 1000      # number of timesteps per output
        self.n_tstep: int = 1_000_000    # total number of timesteps
        
        
        if is_2d:
            self.igridv: int = I_LOGP
        else:
            self.igridv = I_CART

        self.add_gas("H2O")


        
    def set_stepping(self, dt=None, output_gap = None, n_tstep = None) -> None:
        if dt:
            if dt != int(dt):
                raise TypeError("dt must be a integer")
            dt = int(dt)
            if dt < 0:
                raise ValueError("dt must be positive")
        
        if output_gap:
            if output_gap != int(output_gap):
                raise TypeError("dt must be a integer")
            output_gap = int(output_gap)
            if dt <= 0:
                raise ValueError("output_gap must be positive")
            
        if n_tstep:
            if n_tstep != int(n_tstep):
                raise TypeError("n_tstep must be a integer")
            n_tstep = int(n_tstep)
            if n_tstep < 0:
                raise ValueError("n_tstep must be positive")
            
        if dt: self.dt = dt
        if output_gap: self.output_gap = output_gap
        if n_tstep: self.n_tstep = n_tstep
            
        
    def set_physical_params(self, surface_grav = None, wt_mol = None, r_planet = None, velocity_avg = None, use_jovian_radius=False):
        if surface_grav:
            if surface_grav < 0:
                raise ValueError("Surface Gravity must be positive")
       
            
        if wt_mol:
            if wt_mol < 0:
                raise ValueError("Molar Weight must be positive")
            if wt_mol > 3 or wt_mol < 2:
                warnings.warn(f"Typical values of wt_mol are between 2 and 3.  Your value is {wt_mol}.")
           
            
        if r_planet:
            if r_planet < 0:
                raise ValueError("Planet Radius must be positive")
            if r_planet < 20 and not use_jovian_radius:
                warnings.warn("You specified a planetary radius under 20.  Assuming you meant in units of Jovian radius.  \n set use_jovian_radius=True to supress this warning")
                use_jovian_radius = True
            if r_planet > 1e3 and use_jovian_radius:
                raise ValueError(f"The specified planetary radius of {r_planet} jovian radii is too high")
            
        if velocity_avg:
            if velocity_avg < 0:
                raise ValueError("velocity_avg must be positive")
            if not self.is_2d:
                warnings.warn("velocity_avg is ignored in a 1D sim")
        
        
        
        
        
        if surface_grav: self.surface_grav = surface_grav
        if wt_mol: self.wt_mol = wt_mol
        if r_planet:
            if use_jovian_radius:
                self.r_planet = r_planet * JUPITER_RADIUS
            else:
                self.r_planet = r_planet
        if velocity_avg: self.velocity_avg = velocity_avg
                
   
    
    def add_kzz(self, levels: ArrayLike) -> None:
        if self.NZ:
            if len(levels) != self.NZ + 1:
                raise ValueError(f"levels must be {self.NZ+1} long to be compatible with other inputs.\nYour data was {len(levels)} long.")
        else:
            self.NZ = len(levels) - 1
        self.kzz_levels = levels
    
    def add_P(self, levels: ArrayLike) -> None:
        levels = np.array(levels)
        if self.NZ:
          if len(levels) != self.NZ + 1:
                raise ValueError(f"levels must be {self.NZ+1} long to be compatible with other inputs.\nYour data was {len(levels)} long.")
        else:
            self.NZ = len(levels) - 1

        self.P_centers = (levels[:-1] + levels[1:])/2
        self.P_levels = levels
        
    def add_T(self, levels: ArrayLike) -> None:
        levels = np.array(levels)
        if self.NZ:
            if levels.shape[0] != self.NZ + 1:
                raise ValueError(f"levels must be {self.NZ + 1} long to be compatible with other inputs.\nYour data was {len(levels)} long.")
        else:
           self.NZ = levels.shape[0] - 1
        if self.is_2d:
            if len(levels.shape) != 2:
                raise ValueError("Carma is in 2-D mode: T centers must be a 2-D array")
            self.NLONGITUDE = levels.shape[1]

            self.T_centers = (levels[:-1, :] + levels[1:, :])/2

        else:
            if len(levels.shape) != 1:
                raise ValueError("Carma is in 1-D mode: T centers must be a 1-D array")
            self.T_centers = (levels[:-1] + levels[1:])/2
    
        self.T_levels = levels
  


        
    def add_z(self,  levels: ArrayLike) -> None:
        if self.NZ:
             if len(levels) != self.NZ + 1:
                raise ValueError(f"levels must be {self.NZ+1} long to be compatible with other input.\nYour data was {len(levels)} long.")
        else:
            self.NZ = len(levels) - 1
        self.z_centers = (levels[:-1] + levels[1:])/2
        self.z_levels = levels

    
    def add_het_group(self, 
                      gas: str, 
                      seed_group: str, 
                      rmin: float, 
                      mucos: float | None = None, 
                      add_coag: bool = False) -> "Group":
        if type(gas) == type(""):
            gas = self.gasses.get(gas, Gas(gas, len(self.gasses) + 1))
            self.gasses[gas.name] = gas
        
        if type(seed_group) == type(""):
            seed_group = self.groups["Pure "+seed_group.split(" ")[-1]]
        
        name =  gas.name + " on " + seed_group.name.split(" ")[-1]
        group = Group(len(self.groups)+1, name, rmin)
        self.groups[name] = group
   
        if not mucos:
            mucos = mucos_dict[gas.name][seed_group.name.split(" ")[-1]]
        
        self.nucs.append(Nuc(seed_group, group, True, gas, mucos))
        
        mantle_elem = Element(gas.name + " Mantle", len(self.elems)+1, 
                              group, cond_rho[gas.name], "Volatile", 
                              self.gasses[gas.name].igas)
        self.elems[mantle_elem.name] = mantle_elem
        group.mantle = mantle_elem

        core_elem = seed_group.coreify(len(self.elems)+1, group, gas.name)
        self.elems[core_elem.name] = core_elem
        group.core = core_elem
        
        growth = Growth(mantle_elem, gas)
        self.growth.append(growth)
        
        if add_coag:
            self.add_coag(group)
        
        return group
        
    def add_hom_group(self, 
                      gas: str,
                      rmin: float,
                      add_coag: bool = False) -> "Group":
        if type(gas) == type(""):
            gas = self.gasses.get(gas, Gas(gas, len(self.gasses) + 1))
            self.gasses[gas.name] = gas
            
        name = "Pure "+ gas.name
        group = Group(len(self.groups)+1, name, rmin)
        self.groups[name] = group
        

        
        self.nucs.append(Nuc(group, group, False, gas,  0))
        
        elem = Element("Pure "+ gas.name, len(self.elems)+1, 
                            group, cond_rho[gas.name], "Volatile", 
                            self.gasses[gas.name].igas)
        group.core = elem
        self.elems[elem.name] = elem
        self.growth.append(Growth(elem, gas))
        
        if add_coag:
            self.add_coag(group)
        
        return group
    
    def add_gas(self, gas: str, **kwargs) -> "Gas":
        self.gasses[gas] = self.gasses.get(gas, Gas(gas, len(self.gasses)+1, **kwargs))
        return gas
      
    def add_coag(self, group: str):
        if type(group) == type(""):
            g = self.groups.get(group, False)
            if g:
                raise ValueError(f"Group '{group}' not found")
        elif type(group) == type(Group(-1, -1, -1)):
            g = group
        else:
            raise TypeError("Group must be a group object or a string")
            
        self.coags.append(Coag(g))
        
      
    def set_nmr(self, nmr_dict: dict):
        for key in nmr_dict.keys():
            self.gasses[key].nmr = nmr_dict[key]
    
    def calculate_z(self, wt_mol: float | None=None):
        if wt_mol is None:
            wt_mol = self.wt_mol
        if wt_mol is None:
            raise RuntimeError("Carma.wt_mol must be set or a mean molecular weight array must be provided")
        
        if (self.T_levels is None or self.P_levels is None):
            raise RuntimeError("T and P levels must be set")
        if (self.surface_grav is None):
            raise RuntimeError("surface_grav must be set")

        H_levels = k_B * self.T_levels/(wt_mol * PROTON_MASS * self.surface_grav)

        self.z_levels = np.zeros(self.NZ + 1)

        for i in range(1, self.NZ+1):
            dz = H_levels[i] * np.log(self.P_levels[i-1]/self.P_levels[i])
            self.z_levels[i] = self.z_levels[i-1] + dz

        self.z_centers = (self.z_levels[1:] + self.z_levels[:-1])/2



    def extend_atmosphere(self, max_P: float) -> None:
        """
        Extends the atmosphere P, T, z, and k_zz profiles to the specified depth
        adiabatically using the fit from Parmentier et al. (2015) to the equation
        of state described in Saumon (1995).  Requires that the P, T, z, and k_zz
        levels are set

        """
        if (self.P_levels is None or self.T_levels is None or self.kzz_levels is None or self.z_levels is None):
            raise RuntimeError("P_levels, T_levels, z_levels, and/or kzz_levels are not set")
        
        if (self.surface_grav is None or self.wt_mol is None):
            raise RuntimeError("g and/or wt_mol are not set")

        ratio = self.P_levels[0]/self.P_levels[1]
        n = int(np.log(max_P/self.P_levels[0])/np.log(ratio)+1)

        self.NZ = self.NZ + n

        P_new = np.zeros(self.NZ + 1)
        T_new = np.zeros(self.NZ + 1)
        kzz_new = np.zeros(self.NZ + 1)
        z_new = np.zeros(self.NZ + 1)

        P_new[n:] = self.P_levels
        T_new[n:] = self.T_levels
        kzz_new[n:] = self.kzz_levels
        z_new[n:] = self.z_levels

        H0 = k_B * self.T_levels[0]/(self.wt_mol * PROTON_MASS * self.surface_grav)

        def K(P):
            return (self.T_levels[0]/(PARMENTIER_A_COEFF 
                                      - PARMENTIER_B_COEFF * self.T_levels[0])
                        * (P/self.P_levels[0]) ** PARMENTIER_A_COEFF)
    
        def new_T(P):
            return PARMENTIER_A_COEFF * K(P) / (1 + PARMENTIER_B_COEFF * K(P))
        
        for i in range(n-1, -1, -1):
            P_new[i] = self.P_levels[0] * ratio ** (n - i)
            T_new[i] = new_T(P_new[i]) 

            H = k_B * T_new[i]/(self.wt_mol * PROTON_MASS * self.surface_grav)
            kzz_new[i] = self.kzz_levels[0] * (H/H0)**(1/3)

            dz = H * np.log(P_new[i]/P_new[i+1])
            z_new[i] = z_new[i+1] - dz

        z_new -= z_new[0]

        self.z_centers = (z_new[1:] + z_new[:-1])/2
        self.P_centers = (P_new[1:] + P_new[:-1])/2
        self.T_centers = (T_new[1:] + T_new[:-1])/2

        self.z_levels = z_new
        self.P_levels = P_new
        self.T_levels = T_new
        self.kzz_levels = kzz_new


    def calc_H(self, centers=False):
        
        if centers:
            T = self.T_centers
        else:
            T = self.T_levels
       
        return k_B * T/(self.wt_mol * PROTON_MASS * self.surface_grav)

    def run(self, path=None, suppress_output=False):
        if self.is_2d and self.velocity_avg < 0:
            raise RuntimeError("For 2D carma, velocity_avg must be specified")
        
        if (self.wt_mol is None or self.surface_grav is None):
            raise RuntimeError("surface_grav and wt_mol must be set")
        
        if not path: path = self.name
        
        os.makedirs(path, exist_ok=True)
        os.makedirs(os.path.join(path, "inputs"), exist_ok=True)

        shutil.copy(os.path.join(SRC, "carmapy", "carmapy.exe"), path)

        
        path_end = os.path.basename(path) 
        
        nml = {
            "io_files": {
                "filename":             path_end,
                "filename_restart":     path_end+"_restart",
                "fileprefix":           "bd_",
                "gas_input_file":       os.path.join("inputs", "gas_input.txt"),
                "centers_file":         os.path.join("inputs", "centers.txt"),
                "levels_file":          os.path.join("inputs", "levels.txt"),
                "temps_file":           os.path.join("inputs", "temps.txt"),
                "groups_file":          os.path.join("inputs", "groups.txt"),
                "elements_file":        os.path.join("inputs", "elements.txt"),
                "gases_file":           os.path.join("inputs", "gasses.txt"),
                "growth_file":          os.path.join("inputs", "growth.txt"),
                "nuc_file":             os.path.join("inputs", "nucleation.txt"),
                "coag_file":            os.path.join("inputs", "coagulation.txt")
                },
            "physical_params" : {
                "wtmol_air_set":        self.wt_mol,
                "grav_set":             self.surface_grav,
                "rplanet":              self.r_planet,
                "velocity_avg":         self.velocity_avg
                },
            "input_params": {
                "NZ":                   self.NZ,
                "NELEM":                len(self.elems),
                "NGROUP":               len(self.groups),
                "NGAS":                 len(self.gasses),
                "NBIN":                 self.NBIN,
                "NSOLUTE":              1,
                "NWAVE":                0,
                "NLONGITUDE":           self.NLONGITUDE,
                "irestart":             int(self.restart),
                "idiag":                self.idiag,
                "iskip":                self.output_gap,
                "nstep":                self.n_tstep,
                "dtime":                self.dt,
                "NGROWTH":              len(self.growth),
                "NNUC":                 len(self.nucs),
                "NCOAG":                len(self.coags),
                "IS_2D":                int(self.is_2d),
                "igridv":               self.igridv,
                "iappend":              self.iappend
                }        
            }
        nml = f90nml.Namelist(nml)
        nml.write(os.path.join(path, "inputs", "input.nml"), force=True)
        
        io = nml["io_files"]

        with open(os.path.join(path, io["groups_file"]), "w+") as f:
            f.write("name\trmin\n")
            for key in self.groups.keys():
                name = '"'+key + '"'
                f.write(f'{name:24s}{self.groups[key].rmin:.15e}\n')
        
        with open(os.path.join(path, io["gases_file"]), "w+") as f:
            f.write("name\twtmol\tivaprtn\ticomp\twtmol_dif\n")
            for key in self.gasses.keys():
                name = '"'+key + ' Vapor"'
                f.write(f'{name:24s}{self.gasses[key].wtmol:<.6e}\t{self.gasses[key].ivaprtn:2d}\t{self.gasses[key].icomp:2d}\t{self.gasses[key].wtmol_dif:3.4f}\n')
        
        with open(os.path.join(path, io["elements_file"]), "w+") as f:
            f.write("igroup\tname\trho\tprocess\tigas\n")
            for key in self.elems.keys():
                name = '"'+key + '"'
                proc = '"'+self.elems[key].proc + '"'
                f.write(f'{self.elems[key].group.igroup}\t{name:24s}{self.elems[key].rho:2.4f}\t{proc:15s}\t{self.elems[key].igas:2d}\n')
        
        
        
        with open(os.path.join(path, io["nuc_file"]), "w+") as f:
            f.write("ele_from\tele_to\tis_het\tigas\tevap_to\tmucos\n")
            for nuc in self.nucs:
                igas = nuc.gas.igas
                if nuc.is_het:
                    ele_from = nuc.group_from.core.ielem
                    ele_to = nuc.group_to.core.ielem
                    f.write(f'{ele_from:3d}\t{ele_to:3d}\t1\t{igas:3d}\t{ele_from:3d}\t{nuc.mucos:1.8f}\n')

                else:
                    ele_from = nuc.group_from.core.ielem
                    ele_to = ele_from
                    f.write(f'{ele_from:3d}\t{ele_to:3d}\t0\t{igas:3d}\t{0:3d}\t{0:1.8f}\n')
        
        with open(os.path.join(path, io["growth_file"]), "w+") as f:
            f.write("ielem\tigas\n")
            for g in self.growth:
                f.write(f"{g.elem.ielem}\t {g.gas.igas}\n")
        
        with open (os.path.join(path, io["coag_file"]), "w+") as f:
            f.write("igroup\n")
            for c in self.coags:
                f.write(f"{c.group.igroup}\n")

        with open(os.path.join(path, io["centers_file"]), "w+") as f:
            f.write("z_centers\tP_centers\n")
            for i in range(self.NZ):
                f.write(f"{self.z_centers[i]/100}\t{self.P_centers[i]/10}\n")
        
        with open(os.path.join(path, io["levels_file"]), "w+") as f:
            f.write("z_levels\tP_levels\tkzz_levels\n")
            for i in range(self.NZ+1):
                f.write(f"{self.z_levels[i]/100}\t{self.P_levels[i]/10}\t{self.kzz_levels[i]}\n")
        
        np.savetxt(os.path.join(path, io["temps_file"]), self.T_centers, delimiter='\t')
        
        with open(os.path.join(path, io["gas_input_file"]), "w+") as f:
            for key in self.gasses.keys():
                f.write(key+"\t")
            f.write("\n")
            
            for key in self.gasses.keys():
                g = self.gasses[key]
                if type(g) == type(1):
                    if g.nmr < 0:
                        raise AttributeError(f"The nmr for {g.name} was not set.")
                if len(np.shape(g.nmr)) > 0:
                    f.write(f"{g.nmr[0]:10e}\t")
                else:
                    f.write(f"{g.nmr:10e}\t")
            f.write("\n")
            for i in range(1, self.NZ):
                for key in self.gasses.keys():
                    g = self.gasses[key]
                    if len(np.shape(g.nmr)) > 1:
                        if len(g.nmr) != self.NZ:
                            raise ValueError(f"The array for nmr of {g.name} is {len(g.nmr)}.  It should be {self.NZ}.")
                        f.write(f"{g.nmr[i]:10e}\t")
                    else:
                        f.write(f"{0.:10e}\t")
                f.write("\n")
        
        with _cd(path):

            try:
                subprocess.run(["export", "OMP_NUM_THREADS=1"], shell=True,stdout=subprocess.PIPE)
                subprocess.run(["export", "KMP_STACKSIZE=128M"], shell=True,stdout=subprocess.PIPE)
                p = subprocess.Popen(os.path.join(SRC, "carmapy", "carmapy.exe"), shell=False, stdout=subprocess.PIPE)
                
                while p.poll() is None:
                    l = p.stdout.readline() # This blocks until it receives a newline.
                    if not suppress_output: print(l.decode('UTF-8'))
                # When the subprocess terminates there might be unconsumed output 
                # that still needs to be processed.
                if not suppress_output: print(p.stdout.read().decode('UTF-8'))
            except Exception as e:
                print(e)
            
    def read_results(self):
        self.results = Results(self)
        
        
    
    
class Element:
    def __init__(self, 
                 name: str,
                 ielem: int,
                 group: "Group",
                 rho: float,
                 proc: str,
                 igas: int):
        self.name:  str     = name
        self.ielem: int     = ielem
        self.group: "Group" = group
        self.rho:   float   = rho
        self.proc:  str     = proc
        self.igas:  int     = igas
    
class Gas:
    def __init__(self, 
                 gas_name: str, 
                 igas: int, 
                 **kwargs):
        self.name:      str     = gas_name
        self.igas:      int     = igas
        self.wtmol:     float   = kwargs.get("wtmol", wtmol_dict[gas_name])
        self.ivaprtn:   int     = kwargs.get("ivaprtn",vaprtn_dict[gas_name])
        self.icomp:     int     = gcomp_dict[gas_name]
        self.wtmol_dif: float   = kwargs.get("wtmol_dif", wtmol_dif_dict[gas_name])


        self.nmr: float | ArrayLike = kwargs.get("nmr", -1)
        
    
class Nuc:
    def __init__(self, group_from, group_to, is_het, gas, mucos):
        self.group_from = group_from
        self.group_to = group_to
        self.is_het = is_het
        self.gas = gas
        self.ievp2elem = group_from.core
        self.mucos = mucos
    
class Growth:
    def __init__(self, elem, gas):
        self.elem = elem
        self.gas = gas
    
class Group:
    def __init__(self, igroup, name, rmin):
        self.igroup = igroup
        self.name = name
        self.rmin = rmin
        self.core = None
        self.mantle = None
    
    def coreify(self, ielem, group, gas_name="") -> Element:
        core_elem = self.core
        
        name = core_elem.name
        name = name.split(" ")[-1]
        name = name + f" Core ({gas_name})"
        
        elem = Element(name, ielem, group, core_elem.rho, "Core Mass", core_elem.igas)
        return elem


class Coag:
    def __init__(self, group):
        self.group = group




def load_carma(path, restart=1):
    carma = Carma(path)
    carma.restart = restart
    
    nml = f90nml.read(os.path.join(path, "inputs", "inputs.nml"))

    carma.NZ = nml["input_params"]["NZ"] 
    carma.NLONGITUDE = nml["input_params"]["NLONGITUDE"]
    carma.is_2d = nml["input_params"]["IS_2D"]
    carma.igridv - nml["input_params"]["igridv"]
    carma.NBIN = nml["input_params"]["NBIN"]
    carma.output_gap = nml["input_params"]["iskip"]
    carma.n_tstep  = nml["input_params"]["nstep"]
    carma.dt = nml["input_params"]["dtime"]
    
    carma.wt_mol = nml["physical_params"]["wtmol_air_set"]
    carma.surface_grav = nml["physical_params"]["grav_set"]
    carma.r_planet = nml["physical_params"]["rplanet"]
    carma.velocity_avg = nml["physical_params"]["velocity_avg"]
    

    io = nml["io_files"]

    with open(os.path.join(path, io["groups_file"]), "r") as f:
        f.readline()
        for line in f:
            name, rmin = shlex.split(line[:-1])
            carma.groups[name] = Group(len(carma.groups)+1, name, float(rmin))
        
    with open(os.path.join(path, io["gases_file"]), "r") as f:
        f.readline()
        for line in f:
            name, wtmol, ivaprtn, icomp, wtmol_dif = shlex.split(line[:-1])
            name= name[:-len(' Vapor')]
            carma.gasses[name] = Gas(name, len(carma.gasses)+1, wtmol=float(wtmol), ivaprtn=int(ivaprtn), wtmol_dif=float(wtmol_dif))
            
    
    with open(os.path.join(path, io["elements_file"]), "r") as f:
        f.readline()
        for line in f:
            igroup, name, rho, proc, igas = shlex.split(line[:-1])
            group = carma.groups[list(carma.groups.keys())[int(igroup)-1]]
            carma.elems[name] = Element(name, len(carma.elems)+1, group, float(rho), proc, int(igas))
            if "Mantle" in name:
                group.mantle = carma.elems[name]
            else:
                group.core = carma.elems[name]
                
    with open(os.path.join(path, io["nuc_file"])) as f:
        f.readline()
        for line in f:
            ele_from, ele_to, _, igas, _, mucos = shlex.split(line[:-1])
            if ele_to == ele_from:
                is_het = False
            else:
                is_het = True
            group_from = carma.elems[list(carma.elems.keys())[int(ele_from)-1]].group
            group_to = carma.elems[list(carma.elems.keys())[int(ele_to)-1]].group
            gas = carma.gasses[list(carma.gasses.keys())[int(igas)-1]]
            carma.nucs.append(Nuc(group_from, group_to, is_het, gas, float(mucos)))
            
    with open(os.path.join(path, io["growth_file"])) as f:
        f.readline()
        for line in f:
            ielem, igas = shlex.split(line[:-1])
            elem = carma.elems[list(carma.elems.keys())[int(ielem)-1]]
            gas = carma.gasses[list(carma.gasses.keys())[int(igas)-1]]
            carma.growth.append(Growth(elem, gas))
            
    with open(os.path.join(path, io["coag_file"])) as f:
        f.readline()
        for line in f:
            igroup = int(f.readline())
            carma.add_coag(carma.groups[list(carma.elems.keys())[igroup-1]])
            
    centers = np.genfromtxt(os.path.join(path, io["centers_file"]), skip_header=1)
    levels = np.genfromtxt(os.path.join(path, io["levels_file"]), skip_header=1)
    
    carma.z_centers = centers[:,0]*100
    carma.z_levels = levels[:,0]*100

    carma.P_centers = centers[:, 1]*10
    carma.P_levels = levels[:,1]*10

    carma.kzz_levels = levels[:, 2]

    carma.T_centers = np.genfromtxt(os.path.join(path, io["temps_file"]))


    gas_input = np.genfromtxt(os.path.join(path, io["gas_input_file"]))
    for i, key in enumerate(carma.gasses.keys()):
        carma.gasses[key].nmr = gas_input[1:, i]
        
    return carma


def available_species():
    print(list(gas_dict.keys())[1:])


def included_mucos(specie):
    print(gas_dict[specie]["mucos_dict"])