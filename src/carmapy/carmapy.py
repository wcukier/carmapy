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

from typing import Union

# The src/ directory in the package
SRC = os.path.dirname(os.path.dirname(__file__)) 


@contextlib.contextmanager
def _cd(path: str):
   """ Allows for safe cd in python code """
   old_path = os.getcwd()
   os.chdir(path)
   try:
       yield
   finally:
       os.chdir(old_path)

ALLOWED_BCs = ["fixed_conc", "fixed_flux", "zero_grad"]
BC_INTS = [I_FIXED_CONC, I_FLUX_SPEC, I_ZERO_CGRAD]

def _bc2int(bc):
    for i in range(len(ALLOWED_BCs)):
        if bc == ALLOWED_BCs[i]:
            return BC_INTS[i]
    raise ValueError(f"{bc} not found")



class Carma:
    """
    An object representing a CARMA simulation. 

    Parameters
    ----------
    name : str
        A name identifier for the CARMA simulation. This name is used for output 
        directories and file prefixes.
    is_2d : bool, optional
        If True, sets the simulation in 2D mode (longitude x vertical). 
        Defaults to False  (1D vertical column). 
        **Note: `is_2d=True` is not currently supported.**


    Examples
    --------
    >>> carma = Carma("new_carma")
    """

    NZ: int = 0               
    """ Number of vertical layers in the atmospheric grid. """

    NBIN: int = 80             
    """ Number of particle radius bins (default 80). """

    NLONGITUDE: int = -1        
    """ Number of longitude bins (used only in 2D mode). """

    P_levels: ArrayLike = None        
    """ Pressure levels starting from bottom of atmosphere [barye]. """

    P_centers: ArrayLike = None       
    """ Pressure center values (between levels) [barye]. """

    z_levels: ArrayLike = None       
    """ Altitude levels where z=0 is bottom of atmosphere [cm]. """

    z_centers: ArrayLike = None       
    """  Altitude center values (between levels) [cm]. """

    T_levels: ArrayLike = None        
    """  Temperature levels [K]. """

    T_centers: ArrayLike = None      
    """ Temperature center values (between levels) [K]. """

    kzz_levels: ArrayLike = None     
    """ Eddy diffusion coefficient levels [cm²/s]. """

    surface_grav: float = None   
    """ Surface gravity of the planet / brown dwarf [cm/s²] """

    wt_mol: float = None          
    """ Mean molecular weight of the atmosphere [dimentionless] """

    log_metallicity: float = 0
    """ Log Solar Metallicity (Default 0) """

    r_planet: float = 6.991e9    
    """ Planet radius [cm] (Defaults to 1 Jovian radius); ignored in 1D mode. 
    """

    velocity_avg: float = -1    
    """ Mean longitudinal wind speed at the equator [cm/s]; ignored in 1-D mode 
    """

    restart: bool = False        
    """If true, continues the simulation from a previously stored state 
    (Default False) """

    dt: int              
    """ Simulation timestep [s] """

    output_gap: int = 1000      
    """ Number of timesteps between each simulation output """

    n_tstep: int = 1_000_000    
    """ Total number of timesteps """

    is_2d: bool
    """ If true, the simulation is a 2D Carma sim, otherwise 1D"""

    name: str
    """ The name of the directory the carma object is saved / loaded from """

    iappend: int
    """ 1 if and only if a restarted run should add to existing files (default 0)"""
    
    groups: dict[str, "Group"]
    """ Dictionary of group objects indexed by group name used in the sim """

    growth: list["Growth"]
    """ List of the microphysical growth pathways enabled in the sim """

    elems: dict[str, "Element"]
    """ Dictionary of the element objects indexed by elem name used in the sim """

    gases: dict[str, "Gas"]
    """ Dictionary of the gas objects, indexed by gas name, used in the sim """

    nucs: list["Nuc"]
    """ List of the nucleation pathways enabled in the sim """

    coags: list["Coag"]
    """ List of the coagulation pathways enabled in the sim """

    _citation = {}


    def __init__(self, name: str, is_2d=False) -> None:        
        self.is_2d:         bool        = is_2d          # true for 2d carma, false for 1d carma        
        self.name:      str = name            # Name for the carma object, used to define the directory the object is saved in

        self.idiag:         int         = 0
        self.iappend:       int         = 0


        self.groups:    dict[str, "Group"]      = {}      # dictionary of carma Group objects
        self.growth:    list["Growth"]          = []      # list of carma Growth objects
        self.elems:     dict[str, "Element"]    = {}      # dictionary of carma Element objects
        self.gases:    dict[str, "Gas"]    = {}      # dictionary of carma Gas objects
        self.nucs:      list["Nuc"]             = []      # list of carma Nuc objecs
        self.coags:     list["Coag"]            = []      # dictionary of carma Coag objects

        self.winds = None

        self.bot_bound_type_cloud = "fixed_conc"
        self.bot_bound_type_gas = "fixed_conc"
        self.top_bound_type_cloud = "fixed_flux"
        self.top_bound_type_gas = "fixed_flux"

        if is_2d:
            self.igridv: int = I_LOGP
        else:
            self.igridv = I_CART

        self.add_gas("H2O")


        
    def set_stepping(self, 
                     dt: float | None = None, 
                     output_gap: int | None = None, 
                     n_tstep: int | None = None) -> None:
        """Modifies the simulation timestepping behavior.  
        Specifically, sets the simulation timestep, gap between outputs, and 
        total number of timesteps. Will leave as-is any parameter not provided

        Parameters
        ----------
        dt : float, optional
            Simulation timestep [s], left as is if not provided
        output_gap : int, optional
            Gap between simulation outputs, left as is if not provided
        n_tstep : int, optional
            Total number of timesteps to run the simulation for, left as is if 
            not provided
        """
        
  

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
            
        
    def set_physical_params(self, 
                            surface_grav: float = None, 
                            wt_mol: float = None, 
                            log_metallicity: float = None,
                            r_planet: float = None, 
                            velocity_avg: float = None, 
                            use_jovian_radius: bool = False) -> None:
        """Modifies the physical parameters of the simulation.
        Can be used to set the surface gravity, the mean molecular weight, the 
        log metallicity of the planet, the planetary radius, and the average
        longitudinal equatorial wind speed velocity.  Will leave as-is any 
        parameter not provided

        Parameters
        ----------
        surface_grav : float, optional
            Surface gravity [cm/s²], left as is if not provided
        wt_mol : float, optional
            Mean molecular weight [dimensionless], left as is if not provided
        log_metallicity : float, optional
            Log metallicity relative to solar ([Fe/H]), left as is if not 
            provided
        r_planet : float, optional
            The planetary radius [cm (or R_J if use_jovian_radius is true)], 
            left as is if not provided
        velocity_avg : float, optional
            The average longitudinal wind speed at the equator, left as is if 
            not provided
        use_jovian_radius : bool, optional
            Set to true to provide planetary radius in Jovian radii, by default 
            False
        """
        if surface_grav:
            if surface_grav < 0:
                raise ValueError("Surface Gravity must be positive")
       
        if wt_mol:
            if wt_mol < 0:
                raise ValueError("Molar Weight must be positive")
           
        if log_metallicity:
            if ((type(log_metallicity) != type(1)) 
                and (type(log_metallicity) != type(0.5))):
                raise ValueError("log_metallicity must be an int or float")
            
        if r_planet:
            if r_planet < 0:
                raise ValueError("Planet Radius must be positive")
            if r_planet < 20 and not use_jovian_radius:
                warnings.warn("You specified a planetary radius under 20. "
                              "Assuming you meant in units of Jovian radius. "
                        "\n set use_jovian_radius=True to supress this warning")
                use_jovian_radius = True
            if r_planet > 1e3 and use_jovian_radius:
                raise ValueError(f"The specified planetary radius of "
                                 f"{r_planet} jovian radii is too high")
            
        if velocity_avg:
            if velocity_avg < 0:
                raise ValueError("velocity_avg must be positive")
            if not self.is_2d:
                warnings.warn("velocity_avg is ignored in a 1D sim")
        
        
        if surface_grav: self.surface_grav = surface_grav
        if wt_mol: self.wt_mol = wt_mol
        if log_metallicity: self.log_metallicity = log_metallicity
        if r_planet:
            if use_jovian_radius:
                self.r_planet = r_planet * JUPITER_RADIUS
            else:
                self.r_planet = r_planet
        if velocity_avg: self.velocity_avg = velocity_avg
                
   
    
    def add_kzz(self, levels: ArrayLike) -> None:
        """Sets eddy diffusion levels.  Checks to see if the provided levels are
        compatible in array shape with the other inputs already provided to 
        carma. If NZ is not set, sets NZ. The first element
        of the levels array should correspond to the bottom of the atmosphere.

        Parameters
        ----------
        levels : ArrayLike
            Eddy diffusion coefficient values at the boundaries of the 
            atmospheric cells [cm²/s]
        """
        if self.NZ:
            if len(levels) != self.NZ + 1:
                raise ValueError(f"levels must be {self.NZ+1} long to be "
                                 f"compatible with other inputs."
                                 f"\nYour data was {len(levels)} long.")
        else:
            self.NZ = len(levels) - 1
        self.kzz_levels = levels
    
    def add_P(self, levels: ArrayLike) -> None:
        """Sets pressure levels and centers.  Uses provided 
        levels to calculate the pressure centers and checks to see if the
        provided levels are compatible in array shape with the other inputs 
        already provided to carma. If NZ is not set, sets NZ. The first element
        of the levels array should correspond to the bottom of the atmosphere.

        Parameters
        ----------
        levels : ArrayLike
            Pressure values at the boundaries of the atmospheric cells [barye]
        """
        levels = np.array(levels)
        if self.NZ:
          if len(levels) != self.NZ + 1:
                raise ValueError(f"levels must be {self.NZ+1} long to be "
                                 f"compatible with other inputs."
                                 f"\nYour data was {len(levels)} long.")
        else:
            self.NZ = len(levels) - 1

        self.P_centers = (levels[:-1] + levels[1:])/2
        self.P_levels = levels
        
    def add_T(self, levels: ArrayLike) -> None:
        """Sets temperature levels and centers.  Uses provided 
        levels to calculate the temperature centers and checks to see if the
        provided levels are compatible in array shape with the other inputs 
        already provided to carma. If NZ is not set, sets NZ. The first element
        of the levels array should correspond to the bottom of the atmosphere.

        Parameters
        ----------
        levels : ArrayLike
            Temperature values at the boundaries of the atmospheric cells [K].
            For 1-D CARMApy should be of shape (NZ,).  For 2-D CARMApy should
            be of shape (NZ, NLONGITUDEz)
            
        """
        levels = np.array(levels)
        if self.NZ:
            if levels.shape[0] != self.NZ + 1:
                raise ValueError(f"levels must be {self.NZ + 1} long"
                                 + "to be compatible with other inputs."
                                 +"\nYour data was {len(levels)} long.")
        else:
           self.NZ = levels.shape[0] - 1
        if self.is_2d:
            if len(levels.shape) != 2:
                raise ValueError("Carma is in 2-D mode: "
                                 +"T centers must be a 2-D array")
            self.NLONGITUDE = levels.shape[1]

            self.T_centers = (levels[:-1, :] + levels[1:, :])/2

        else:
            if len(levels.shape) != 1:
                raise ValueError("Carma is in 1-D mode: "
                                 +"T centers must be a 1-D array")
            self.T_centers = (levels[:-1] + levels[1:])/2
    
        self.T_levels = levels
  


        
    def add_z(self,  levels: ArrayLike) -> None:
        """Sets altitude levels and centers.  Uses provided 
        levels to calculate the altitude centers and checks to see if the
        provided levels are compatible in array shape with the other inputs 
        already provided to carma. If NZ is not set, sets NZ. z=0 should be the
        first element in the levels array which corresponds to the bottom of the 
        atmosphere.

        Parameters
        ----------
        levels : ArrayLike
            Altitude values at the boundaries of the atmospheric cells [cm]
        """

        if self.NZ:
             if len(levels) != self.NZ + 1:
                raise ValueError(f"levels must be {self.NZ+1} long to be "
                                 f"compatible with other inputs."
                                 f"\nYour data was {len(levels)} long.")
        else:
            self.NZ = len(levels) - 1
        self.z_centers = (levels[:-1] + levels[1:])/2
        self.z_levels = levels

    def add_vertical_winds(self, wind_centers: ArrayLike):
        """Sets the vertical wind speeds as a function of altitude.  Note that
        unlike most other functions, this one requires wind speeds on the centers,
        not the levels, of the atmosphere

        Parameters
        ----------
        wind_centers : ArrayLike
            The vertical wind speed at each altitude center [cm/s]

        """
        if self.NZ:
            if len(wind_centers) != self.NZ:
                raise ValueError(f"wind_centers must be {self.NZ} long to be "
                                f"compatible with other inputs."
                                f"\nYour data was {len(wind_centers)} long.")
        else:
            self.NZ = len(wind_centers)
        
        self.wind_centers = wind_centers

    
    def add_het_group(self, 
                      gas: Union[str, "Gas"], 
                      seed_group: Union[str, "Group"], 
                      rmin: float, 
                      mucos: float = None, 
                      add_coag: bool = False) -> "Group":
        """Adds a heterogeneously nucleating group to the simulation.  Assumes
        the gas nucleates on the seed particle. If a string is passed to `gas`, 
        will use the carmapy default parameters for that gas if that gas does 
        not already exist in the simulation.  Additionally adds the gas and any 
        created elements, nuc objects, growth objects, and coag objects to the 
        simulation.

        Parameters
        ----------
        gas : Union[str, Gas]
            Gas object, name of already created gas, or name of the default 
            carmapy gas that nucleates on the seed particle
        seed_group : Union[str, Group]
            Group object, name of group, or name of gas that formed the group
            which serves as the seed particle for the condensate
        rmin : float
            Minimum radius of the condensate [cm]
        mucos : float, optional
            Cosine of the contact angle between the gas and the seed particle, 
            if not provided defaults to carmapy defaults
        add_coag : bool, optional
            If true, allows coagulation of this particle onto itself,
            by default False

        Returns
        -------
        Group
            The created group consisting of the gas nucleated onto the seed 
            particle
        """
        if type(gas) == type(""):
            gas = self.gases.get(gas, gas)
            if type(gas) == type(""): gas = Gas(gas, len(self.gases) + 1)
        self.gases[gas.name] = gas
        
        if type(seed_group) == type(""):
            seed_group = self.groups["Pure "+seed_group.split(" ")[-1]]
        
        name =  gas.name + " on " + seed_group.name.split(" ")[-1]
        group = Group(len(self.groups)+1, name, rmin)
        self.groups[name] = group
   
        if not mucos:
            mucos = gas_dict[gas.name]["mucos_dict"][seed_group.name.split(" ")[-1]]
        
        self.nucs.append(Nuc(seed_group, group, True, gas, mucos))
        
        mantle_elem = Element(gas.name + " Mantle", len(self.elems)+1, 
                              group, gas.rho_cond, "Volatile", 
                              self.gases[gas.name].igas)
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
                      gas: Union[str, "Gas"],
                      rmin: float,
                      add_coag: bool = False) -> "Group":
        """Adds a heterogeneously nucleating group to the simulation.  Assumes
        the gas homogeneously nucleates. If a string is passed to `gas`, 
        will use the carmapy default parameters for that gas if that gas does 
        not already exist in the simulation.  Additionally adds the gas and any 
        created elements, nuc objects, growth objects, and coag objects to the 
        simulation.

        Parameters
        ----------
        gas : Union[str, Gas]
            Gas object, name of already created gas, or name of the default 
            carmapy gas that homogeneously nucleates
        rmin : float
             Minimum radius of the condensate [cm]
        add_coag : bool, optional
            If true, allows coagulation of this particle onto itself,
            by default False

        Returns
        -------
        Group
            The created group consisting of the homogeneously nucleated gas
        """
        if type(gas) == type(""):
            gas = self.gases.get(gas, Gas(gas, len(self.gases) + 1))
        self.gases[gas.name] = gas
            
        name = "Pure "+ gas.name
        group = Group(len(self.groups)+1, name, rmin)
        self.groups[name] = group
        

        
        self.nucs.append(Nuc(group, group, False, gas,  0))
        
        elem = Element("Pure "+ gas.name, len(self.elems)+1, 
                            group, gas_dict[gas.name]["rho_cond"], "Volatile", 
                            self.gases[gas.name].igas)
        group.core = elem
        self.elems[elem.name] = elem
        self.growth.append(Growth(elem, gas))
        
        if add_coag:
            self.add_coag(group)
        
        return group
    
    def add_gas(self, gas: Union[str, "Gas"], **kwargs) -> "Gas":
        """Adds a gas to the simulation.  
        
        If a string is passed gas then if the gas does not already exist in the 
        simulation, the carmapy default properties of that gas will be used.  
        Passes ``**kwargs`` to the Gas constructor if a new gas is created.

        Parameters
        ----------
        gas : Union[str, Gas]
            Gas object or name of the carmapy default gas to be added to the
            simulation

        Returns
        -------
        Gas
            The gas which was added to the simulation.
        """
        if type(gas) == type(""):
            self.gases[gas] = self.gases.get(gas, 
                                               Gas(gas, len(self.gases)+1, 
                                                **kwargs))
        else:
            self.gases[gas] = Gas
        return gas
      
    def add_coag(self, group: Union[str, "Group"]) -> None:
        """Adds self coagulation of the given group to the simulation. 

        Parameters
        ----------
        group : Union[str, Group]
            Group object or name of group to add self-coagulation to
        """
        if type(group) == type(""):
            g = self.groups.get(group, False)
            if not g:
                raise ValueError(f"Group '{group}' not found in simulation")
        elif type(group) == type(Group(-1, -1, -1)):
            g = group
            self.groups[group.name] = g
        else:
            raise TypeError("Group must be a group object or a string")
            
        self.coags.append(Coag(g))
        
      
    def set_nmr(self, nmr_dict: dict[str, ArrayLike]) -> None:
        """Sets the gas number mixing ratios of the simulation

        Parameters
        ----------
        nmr_dict : dict[str, ArrayLike]
            The number mixing ratio for each species.  The dictionary keys 
            should be the name of gases in the simulation.  The dictionary 
            values should either be a float representing the mixing ratio at the
            bottom of the atmosphere or an array representing the mixing ratio
            at each atmospheric center
        """
        for key in nmr_dict.keys():
            self.gases[key].nmr = nmr_dict[key]
    


        

    def calculate_z(self, wt_mol: float | ArrayLike =None) -> None:
        """Calculate and set the altitude structure in the proper format for the
        version of carma (1D or 2D) used.  Uses P and T levels to 
        calculate the altitude structure using scale heights.  

        Parameters
        ----------
        wt_mol : float | Arraylike, optional
            Mean molecular weight of the atmosphere.  If an array, must be the
            same length as "levels" arrays, if not provided instead uses the
            mean molecular weight stored in the simulation.
        """
        if wt_mol is None:
            wt_mol = self.wt_mol
        if wt_mol is None:
            raise RuntimeError("Carma.wt_mol must be set or a mean molecular "
                               "weight array must be provided")
        
        if (self.T_levels is None or self.P_levels is None):
            raise RuntimeError("T and P levels must be set")
        if (self.surface_grav is None):
            raise RuntimeError("surface_grav must be set")

        if self.is_2d:
            H_levels = k_B * np.mean(self.T_levels[:, :], axis=1)/(wt_mol 
                                            * PROTON_MASS 
                                            * self.surface_grav)
        else:
            H_levels = k_B * self.T_levels/(wt_mol 
                                * PROTON_MASS 
                                * self.surface_grav)

        self.z_levels = np.zeros(self.NZ + 1)
        self.z_centers = np.zeros(self.NZ)

        if self.igridv == I_CART:
            if self.is_2d:
                warnings.warn("carma.igridv is set to I_CART.  For 2d carma it"
                              "usually is set to I_LOGP.")
            for i in range(1, self.NZ+1):
                dz = H_levels[i] * np.log(self.P_levels[i-1]/self.P_levels[i])
                self.z_levels[i] = self.z_levels[i-1] + dz

            self.z_centers = (self.z_levels[1:] + self.z_levels[:-1])/2

        elif self.igridv == I_LOGP:
            for i in range(self.NZ):
                self.z_levels[i+1] = (H_levels[0] 
                        * np.abs(np.log(self.P_levels[i+1]/self.P_levels[0])))
                self.z_centers[i] =  (H_levels[0] 
                        * np.abs(np.log(self.P_centers[i]/self.P_levels[0])))


        else:
            raise ("carma.igridv not set (This should not be reached unless"
                   "you messed with carma's fields manually).")


    def extend_atmosphere(self, max_P: float, #TODO: See if can do non-iteratively b/c now calculating z later
                          wt_mol=None, 
                          method="adiabatic") -> None:
        """Extends the atmospheric structure to deeper pressures.  Modifies the
        pressure, temperature, and eddy diffusion levels and requires
        that they have previously been set.  If ``max_P`` is less than the 
        current maximum pressure, this method does nothing

        Parameters
        ----------
        max_P : float
            The pressure to which the atmosphere is extended
        wt_mol : ArrayLike, optional
            The mean molecular weight of the atmosphere.  If an array, each 
            entry corresponds to one altitude level. Defaults to the mean 
            molecular weight stored in the carma object
        method : string, optional
            The method to extend the atmosphere.  Options are "adiabatic"
            and "isothermal"

        Notes
        -------
        Atmosphere is extended adiabatically using the fit from Parmentier et 
        al. (2015) [1]_ to the equation of state described in Saumon (1995) 
        [2]_.  k_zz is assumed to be proportional to the cube root of the scale
        height


        References
        ----------
        .. [1] Parmentier, V., Guillot, T., Fortney, J. J., & Marley, 
           M. S. 2015, A&A, 574, A35

        .. [2] Saumon, D., Chabrier, G., & van Horn, H. M. 1995, The 
           Astrophysical Journal Supplement Series, 99 (IOP), 713
        """

        if ((self.P_levels is None) or (self.T_levels is None)
             or (self.kzz_levels is None)):

            raise RuntimeError("P_levels, T_levels, z_levels, "
                               "and/or kzz_levels are not set")
        
        if (self.surface_grav is None or self.wt_mol is None):
            raise RuntimeError("g and/or wt_mol are not set")

        if wt_mol is None: wt_mol = self.wt_mol

        ratio = self.P_levels[0]/self.P_levels[1]
        n = int(np.log(max_P/self.P_levels[0])/np.log(ratio)+1)
        if n < 0: return

        self.NZ = self.NZ + n

        P_new = np.zeros(self.NZ + 1)
        kzz_new = np.zeros(self.NZ + 1)

        P_new[n:] = self.P_levels
        kzz_new[n:] = self.kzz_levels

        if self.is_2d:
            T_new = np.zeros((self.NZ + 1, self.NLONGITUDE))
            T_new[n:, :] = self.T_levels
            H0 = k_B * (np.mean(self.T_levels[0, :])
                        / (self.wt_mol * PROTON_MASS * self.surface_grav))
        else:
            T_new = np.zeros(self.NZ + 1)
            T_new[n:] = self.T_levels
            H0 = k_B * (self.T_levels[0]
                        / (self.wt_mol * PROTON_MASS * self.surface_grav))

        def K(P, t0, p0):
            return (t0/(PARMENTIER_A_COEFF 
                                      - PARMENTIER_B_COEFF * t0)
                        * (P/p0) ** PARMENTIER_A_COEFF)
    
        def new_T(P, t0, p0):
            return (PARMENTIER_A_COEFF * K(P, t0, p0) 
                    / (1 + PARMENTIER_B_COEFF * K(P, t0, p0)))

        if method == "adiabatic":
            for i in range(n-1, -1, -1):
                P_new[i] = self.P_levels[0] * ratio ** (n - i)
                
                if not self.is_2d:
                    T_new[i] = new_T(P_new[i], 
                                     self.T_levels[0], 
                                     self.P_levels[0]) 

                    H = k_B * T_new[i]/(self.wt_mol 
                                        * PROTON_MASS 
                                        * self.surface_grav)
                    

                else:
                    for j in range(self.NLONGITUDE):
                        T_new[i, j] = new_T(P_new[i], 
                                            self.T_levels[0, j], 
                                            self.P_levels[0])
                        
                        H = (k_B * np.mean(T_new[i, :])
                            / (self.wt_mol * PROTON_MASS * self.surface_grav))
                        
                kzz_new[i] = self.kzz_levels[0] * (H/H0)**(1/3)
        
        elif method == "isothermal":
            for i in range(n-1, -1, -1):
                P_new[i] = self.P_levels[0] * ratio ** (n - i)
                
                if not self.is_2d:
                    T_new[i] = self.T_levels[0]
                else:
                    for j in range(self.NLONGITUDE):
                        T_new[i, j] = self.T_levels[0, j]
                        
                kzz_new[i] = self.kzz_levels[0] 



        # z_new -= z_new[0]

        # self.z_centers = (z_new[1:] + z_new[:-1])/2
        self.P_centers = (P_new[1:] + P_new[:-1])/2
        
        if not self.is_2d:
            self.T_centers = (T_new[1:] + T_new[:-1])/2
        else:
            self.T_centers = (T_new[1:, :] + T_new[:-1, :])/2

        # self.z_levels = z_new
        self.P_levels = P_new
        self.T_levels = T_new
        self.kzz_levels = kzz_new

        wt_mol_new = np.zeros(self.NZ + 1)
        wt_mol_new[n:] = wt_mol
        wt_mol_new[:n] = wt_mol_new[n]
        self.calculate_z(wt_mol=wt_mol_new)


    def calc_scale_height(self, calc_centers=False) -> np.ndarray:
        """Calculates the atmospheric scale height.

        Parameters
        ----------
        calc_centers : bool, optional
            if True calculates scale height at centers, otherwise calculates it
            at levels, by default False

        Returns
        -------
        np.ndarray
            The scale height of the atmosphere at each center/level point
        """
        if calc_centers:
            T = self.T_centers
        else:
            T = self.T_levels
       
        return k_B * T/(self.wt_mol * PROTON_MASS * self.surface_grav)

    def set_atmospheric_parameters(self, 
                                   rmu_0: float, 
                                   rmu_t0: float, 
                                   rmu_c: float,
                                   thcond_0: float,
                                   thcond_1: float,
                                   thcond_2: float,
                                   cp: float) -> None:
        """Sets the atmospheric viscosity and thermal conductivity

        Notes
        -----
        1. Atmospheric viscosity, `rmu`, is set by the following formula
           (the Sutherland equation): 
           `rmu = rmu_0 * ((rmu_t0 + tmu_c)/(T + rmu_t0)) * (T/rmu_t0) ** 1.5`
           where `T` is the local atmospheric temperature
        
        2. Atmosphgeric thermal conductivity, `thcond` is set by the following
           formula:
           `thcond = thcond_0 + thcond_1*T + thcond_2 * T**2`

        Parameters
        ----------
        rmu_0 : float
            Viscosoty scaling term [Poise]
        rmu_t0 : float
             Viscosity reference temp [K]
        rmu_c : float
            Viscosity Sutherland constant [K]
        thcond_0 : float
            Consant thermal conductivity term [ergs/s/cm/K]
        thcond_1 : float
            Coefficient to linear thermal conductivity term [ergs/s/cm/K^2]
        thcond_2 : float
            Coefficient to quadratic thermal conductivity term [ergs/s/cm/K^3]
        cp : float
            Specific heat capacity of the atmosphere [erg/g/K]
        """
        self.atmo = {
            "rmu_0": rmu_0,
            "rmu_t0": rmu_t0,
            "rmu_c": rmu_c,
            "thcond_0": thcond_0,
            "thcond_1": thcond_1,
            "thcond_2": thcond_2,
            "CP": cp
        }

    def set_atmospheric_parameters_from_defaults(self, default: str) -> None:
        """Sets the atmospheric viscosity and thermal conductivity from a 
        default profile.  Currently available defaults are:

        - "Pure H2"


        Parameters
        ----------
        default : str
            The name of the default profile to use
        """
        profile = atmo_dict.get(default, -1)

        if profile == -1:
            raise KeyError(f"{default} is not a known profile.  Known profiles"
                           f"are {list(atmo_dict.keys())}.")
        
        self.atmo = profile


    def set_cloud_boundary_type(self,
                                top_boundary_type: str,
                                bottom_boundary_type: str,
                                ) -> None:
        """Sets the type of boundary conditions on the cloud particle 
        concentration.  Note that the same boundary condition type must be used 
        for all cloud species (although different types can be chosen for the 
        bottom of the atmosphere vs the top)

        Parameters
        ----------
        top_boundary_type : str
            Which type of boundary condtion to use at the top of the atmosphere.
            Options are "fixed_conc", "fixed_flux", or "zero_grad"
        bottom_boundary_type : str
            Which type of boundary condtion to use at the bottom of the 
            atmosphere. Options are "fixed_conc", "fixed_flux", or "zero_grad"
        """

        if top_boundary_type not in ALLOWED_BCs: 
            raise ValueError(f"top_boundary_type must be one of {ALLOWED_BCs}")

        if bottom_boundary_type not in ALLOWED_BCs: 
            raise ValueError(f"bottom_boundary_type must be one of "
                             f"{ALLOWED_BCs}")

        self.top_bound_type_cloud = top_boundary_type
        self.bot_bound_type_cloud = bottom_boundary_type

    def set_gas_boundary_type(self,
                                    top_boundary_type: str,
                                    bot_boundary_type: str,
                                    ) -> None:
        """Sets type of boundary conditions on the gas concentration.  
        Note that the same boundary condition type must be used for all 
        cloud species (although different types can be chosen for the bottom of
        the atmosphere vs the top)

        Parameters
        ----------
        top_boundary_type : str
            Which type of boundary condtion to use at the top of the atmosphere.
            Options are "fixed_conc", "fixed_flux", or "zero_grad"
        bot_boundary_type : str
            Which type of boundary condtion to use at the bottom of the 
            atmosphere. Options are "fixed_conc", "fixed_flux", or "zero_grad"
        """
        ALLOWED_BCs = ["fixed_conc", "fixed_flux", "zero_grad"]

        if top_boundary_type not in ALLOWED_BCs: 
            raise ValueError(f"top_boundary_type must be one of {ALLOWED_BCs}")

        if bottom_boundary_type not in ALLOWED_BCs: 
            raise ValueError(f"bottom_boundary_type must be one of "
                             f"{ALLOWED_BCs}")

        self.top_bound_type_gas = top_boundary_type
        self.bot_bound_type_gas = bottom_boundary_type
        
    def set_cloud_boundary(self, 
                            group: Union[str, "Group"],
                            bot_conc=0.0, 
                            top_conc=0.0, 
                            bot_flux=0.0, 
                            top_flux=0.0) -> None:
        """Sets the boundary conditions for the group.  Throws an error if used
        on a multi-element group (ex. a heterogeneously nucleated group)

        Parameters
        ----------
        group : Union[str, Group]
            The group to set the boundary conditions for
        bot_conc : ArrayLike, optional
            Either 0 or an array of NBIN elements describing
            the concentration of each bin in the group at the bottom of the 
            atmosphere. Only used if the bottom cloud boundary condition is set 
            to "fixed_conc". [particles/cm^3]. By default 0 for all bins
        top_conc : ArrayLike, optional
            Either 0 or an array of NBIN elements describing
            the concentration of the group at the top of the atmosphere for each  
            bin. Only used if the top cloud boundary condition is set 
            to "fixed_conc". [particles/cm^3]. By default 0 for all bins
        bot_flux : ArrayLike, optional
            Either 0 or an array of NBIN elements describing
            the upwards flux of the group at the bottom of the atmosphere for  
            each bin. Only used if the bottom cloud boundary condition is set to 
            "fixed_flux". [particles/cm^2/s]. By default 0 for all bins
        top_flux : ArrayLike, optional
            Either 0 or an array of NBIN elements describing
            the concentration of the group at the bottom of the atmosphere for  
            each bin. Only used if the cloud boundary condition is set to 
            "fixed_conc". [particles/cm^3]. By default 0 for all bins.
        """

        if type(group) == type(""):
            group = self.groups[group]

        if group.mantle: raise ValueError("Multi-element groups must have boundary"
                                        "conditions of 0.  The is already the"
                                        "default behavior so this function does"
                                        "not need to be called.")

        if np.all(top_conc == 0): top_conc = np.zeros(self.NBIN)
        if np.all(bot_conc == 0): bot_conc = np.zeros(self.NBIN)
        if np.all(top_flux == 0): top_flux = np.zeros(self.NBIN)
        if np.all(bot_flux == 0): bot_flux = np.zeros(self.NBIN)

        if np.shape(top_conc) != (self.NBIN, ): 
            raise ValueError("top_conc must be 0 or an arraylike of length NBIN")
        if np.shape(bot_conc) != (self.NBIN, ): 
            raise ValueError("bot_conc must be 0 or an arraylike of length NBIN")
        if np.shape(top_flux) != (self.NBIN, ): 
            raise ValueError("top_flux must be 0 or an arraylike of length NBIN")
        if np.shape(bot_flux) != (self.NBIN, ): 
            raise ValueError("bot_flux must be 0 or an arraylike of length NBIN")

        group.boundary["top_conc"] = top_conc
        group.boundary["bot_conc"] = bot_conc
        group.boundary["top_flux"] = top_flux
        group.boundary["bot_flux"] = bot_flux


    def set_gas_boundary(self, 
                            gas: Union[str, "Group"],
                            bot_conc=None, 
                            top_conc=None, 
                            bot_flux=None, 
                            top_flux=None) -> None:
        """Sets the boundary conditions for the gas.  Will leave any unspecified
        boundary conditions unchanged, potentially remaining with the default
        behavior for the gas (see Gas).

        Parameters
        ----------
        gas : Union[str, Gas]
            gas group to set the boundary conditions for
        bot_conc : float, optional
            The number mixing ratio of the gas at the bottom of the 
            atmosphere. Only used if the bottom gas boundary condition is set 
            to "fixed_conc". [particles/particles]. 
        top_conc : float, optional
            The number mixing ratio of the gas at the top of the atmosphere 
            Only used if the top gas boundary condition is set 
            to "fixed_conc". [particles/particles]. 
        bot_flux : float, optional
            The upwards flux of the gas at the bottom of the atmosphere for.
            Only used if the bottom gas boundary condition is set to 
            "fixed_flux". [particles/cm^2/s]. 
        top_flux : ArrayLike, optional
            The concentration of the gas at the bottom of the atmosphere.
            Only used if the cloud boundary condition is set to 
            "fixed_conc". [particles/cm^3].
        """

        if top_conc is not None:
            if top_conc < 0: raise ValueError("top_conc must be positive")
            gas.boundary["top_conc"] = top_conc

        if bot_conc is not None:
            if bot_conc < 0: raise ValueError("bot_conc must be positive")
            gas.boundary["bot_conc"] = bot_conc

        if top_flux is not None:
            gas.boundary["top_flux"] = top_flux

        if bot_flux is not None:
            gas.boundary["bot_flux"] = bot_flux

    def run(self, suppress_output=False) -> None:
        """Runs the CARMA Simulation.

        Creates a directory at the path described by the name of the simulation
        and populates it with the input files required to run the CARMA 
        executable.  Runs and print the stdout from the CARMA executable unless
        suppress_output is true.  The carma executable will also write to output
        files in the created directory.

        Parameters
        ----------
        suppress_output : bool, optional
            If true, will not print stdout from the CARMA executable,
            by default False

        """
        if self.is_2d and self.velocity_avg < 0:
            raise RuntimeError("For 2D carma, velocity_avg must be specified")
        
        if (self.wt_mol is None or self.surface_grav is None):
            raise RuntimeError("surface_grav and wt_mol must be set")
        
        if self.is_2d:
            self._citation["2d"] = True

        path = self.name
        
        os.makedirs(path, exist_ok=True)
        os.makedirs(os.path.join(path, "inputs"), exist_ok=True)

        shutil.copy(os.path.join(SRC, "carmapy", "carmapy.exe"), path)

        
        path_end = os.path.basename(path) 
        
        nml = {
            "io_files": {
                "filename":            path_end,
                "filename_restart":    path_end+"_restart",
                "fileprefix":          "",
                "gas_input_file":      os.path.join("inputs", "gas_input.txt"),
                "centers_file":        os.path.join("inputs", "centers.txt"),
                "levels_file":         os.path.join("inputs", "levels.txt"),
                "temps_file":          os.path.join("inputs", "temps.txt"),
                "groups_file":         os.path.join("inputs", "groups.txt"),
                "elements_file":       os.path.join("inputs", "elements.txt"),
                "gases_file":          os.path.join("inputs", "gases.txt"),
                "growth_file":         os.path.join("inputs", "growth.txt"),
                "nuc_file":            os.path.join("inputs", "nucleation.txt"),
                "coag_file":           os.path.join("inputs", "coagulation.txt"),
                "winds_file":          os.path.join("inputs", "winds.txt"),
                "p_boundary_file":     os.path.join("inputs", "pbound.txt"),
                "g_boundary_file":     os.path.join("inputs", "gbound.txt"),
                },
            "physical_params" : {
                "wtmol_air_set":        self.wt_mol,
                "grav_set":             self.surface_grav,
                "rplanet":              self.r_planet,
                "velocity_avg":         self.velocity_avg,
                "met":                  self.log_metallicity,
                "rmu_0":                self.atmo["rmu_0"],
                "rmu_t0":               self.atmo["rmu_t0"],
                "rmu_c":                self.atmo["rmu_c"],
                "thcond_0":             self.atmo["thcond_0"],
                "thcond_1":             self.atmo["thcond_1"],
                "thcond_2":             self.atmo["thcond_2"],
                "CP":                   self.atmo["CP"]
                },
            "input_params": {
                "NZ":                   self.NZ,
                "NELEM":                len(self.elems),
                "NGROUP":               len(self.groups),
                "NGAS":                 len(self.gases),
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
                "iappend":              self.iappend,
                "itbnd_pc":             _bc2int(self.top_bound_type_cloud),
                "ibbnd_pc":             _bc2int(self.bot_bound_type_cloud),
                "itbnd_gc":             _bc2int(self.top_bound_type_gas),
                "ibbnd_gc":             _bc2int(self.bot_bound_type_gas)
                }        
            }
        nml = f90nml.Namelist(nml)
        nml.write(os.path.join(path, "inputs", "input.nml"), force=True)
        
        io = nml["io_files"]

        with open(os.path.join(path, io["groups_file"]), "w+") as f:
            f.write("name\trmin\n")
            for key in self.groups.keys():
                name = '"'+key + '"'
                f.write(f'{name:35s}{self.groups[key].rmin:.15e}\n')
        
        with open(os.path.join(path, io["gases_file"]), "w+") as f:

            f.write("name\twtmol\tivaprtn\ticomp\twtmol_dif\trho_cond\t"
                    "surften_0\tcoldia\tvp_offset\tvp_tcoeff\tis_type3\t"
                    "surften_slope\tvp_metcoeff\tvplogpcoeff\tlat_heat_e\t"
                    "stofact\n")
            
            for key in self.gases.keys():
                gas : "Gas" = self.gases[key]
                name = '"'+key + ' Vapor"'

                if   gas.gcomp == 1: vaprtn = I_VAPRTN_H2O_MURPHY2005
                elif gas.gcomp == 2: vaprtn = I_VAPRTN_H2SO4_AYERS1980
                else: vaprtn = I_VAPRTN_USER

                f.write(f'{name:24s}{gas.wtmol:<.6e}\t'
                        # f'{gas_dict[key]["rtn"]:2d}\t' 
                        f'{vaprtn:2d}\t' 
                        f'{gas.gcomp:2d}\t'
                        f'{gas.wtmol_dif:<.18e}\t'
                        f'{gas.rho_cond:<.18e}\t'
                        f'{gas.surften_0:<.18e}\t'
                        f'{gas.coldia:<.18e}\t'
                        f'{gas.vp_offset:<.18e}\t'
                        f'{gas.vp_tcoeff:<.18e}\t'
                        f'{int(gas.is_typeIII):1d}\t'
                        f'{gas.surften_slope:<.18e}\t'
                        f'{gas.vp_metcoeff:<.18e}\t'
                        f'{gas.vp_logpcoeff:<.18e}\t'
                        f'{gas.lat_heat_e:<.18e}\t'
                        f'{gas.stofact:2d}'
                        '\n')
        
        # TODO: check for names which are too long
        with open(os.path.join(path, io["elements_file"]), "w+") as f:
            f.write("igroup\tname\trho\tprocess\tigas\n")
            for key in self.elems.keys():
                name = '"'+key + '"'
                proc = '"'+self.elems[key].proc + '"'
                f.write(f'{self.elems[key].group.igroup}\t{name:35s}'
                        f'{self.elems[key].rho:2.4f}\t{proc:15s}\t'
                        f'{self.elems[key].igas:2d}\n')
        
        
        
        with open(os.path.join(path, io["nuc_file"]), "w+") as f:
            f.write("ele_from\tele_to\tis_het\tigas\tevap_to\tmucos\n")
            for nuc in self.nucs:
                igas = nuc.gas.igas
                if nuc.is_het:
                    ele_from = nuc.group_core.core.ielem
                    ele_to = nuc.group_mantle.core.ielem
                    f.write(f'{ele_from:3d}\t{ele_to:3d}\t1\t{igas:3d}\t'
                            f'{ele_from:3d}\t{nuc.mucos:1.8f}\n')

                else:
                    ele_from = nuc.group_core.core.ielem
                    ele_to = ele_from
                    f.write(f'{ele_from:3d}\t{ele_to:3d}\t0\t{igas:3d}\t'
                            f'{0:3d}\t{0:1.8f}\n')
        
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
                f.write(f"{self.z_levels[i]/100}\t"
                        f"{self.P_levels[i]/10}\t"
                        f"{self.kzz_levels[i]}\n")
        
        np.savetxt(os.path.join(path, io["temps_file"]),
                    self.T_centers, 
                    delimiter='\t')

        np.savetxt(os.path.join(path, "inputs", "temp_levels.txt"),
                    self.T_levels, 
                    delimiter='\t')
        
        if self.winds is None:
            self.winds = np.zeros(self.NZ, dtype=float)
        
        with open(os.path.join(path, io["winds_file"]), "w+") as f:
            for w in self.winds:
                f.write(f'{w}\n')

        gas_conc_bot_bc = np.zeros(len(self.gases))

        with open(os.path.join(path, io["gas_input_file"]), "w+") as f:
            for key in self.gases.keys():
                f.write(key+"\t")
            f.write("\n")
            
            for i, key in enumerate(self.gases.keys()):
                g = self.gases[key]
                if type(g) == type(1):
                    if g.nmr < 0:
                        raise AttributeError(f"The nmr for {g.name} "
                                             "was not set.")
                if len(np.shape(g.nmr)) > 0:
                    f.write(f"{g.nmr[0]:10e}\t")
                    gas_conc_bot_bc[i] = g.nmr[0] * g.wtmol_dif/self.wt_mol + 1e-50
                else:
                    f.write(f"{g.nmr:10e}\t")
                    gas_conc_bot_bc[i] = g.nmr * g.wtmol_dif/self.wt_mol + 1e-50

            f.write("\n")
            for i in range(1, self.NZ):
                for key in self.gases.keys():
                    g = self.gases[key]
                    if len(np.shape(g.nmr)) > 1:
                        if len(g.nmr) != self.NZ:
                            raise ValueError(f"The array for nmr of {g.name} "
                                             f"is {len(g.nmr)}.  It should be "
                                             f"{self.NZ}.")
                        
                        f.write(f"{g.nmr[i]:10e}\t")
                    else:
                        f.write(f"{0.:10e}\t")
                f.write("\n")
        
        gas_conc_top_bc = np.zeros(len(self.gases))
        gas_flux_bot_bc = np.zeros(len(self.gases))
        gas_flux_top_bc = np.zeros(len(self.gases))

        ## handle gas boundary conditions
        for i, key in enumerate(self.gases.keys()):
            g = self.gases[key] 
            bc = g.boundary.get("bot_conc", -1)
            if bc != -1:
                gas_conc_bot_bc[i] = bc * g.wtmol_dif/self.wt_mol + 1e-50
            
            gas_flux_bot_bc[i] = g.boundary.get("bot_flux", 0)
            
            bc = g.boundary.get("top_conc", -1)
            if (bc == -1) and (self.top_bound_type_gas == "fixed_conc"): 
                raise RuntimeError(f"'top_conc' gas boundary condition not"
                                f"set for {key} even though"
                                "carma.top_bound_type_gas = 'fixed_conc'")
            
            gas_conc_top_bc[i] = bc * g.wtmol_dif/self.wt_mol + 1e-50
            gas_flux_top_bc[i] = g.boundary.get("top_flux", 0)
        
        with open(os.path.join(path, io["g_boundary_file"]), "w+") as f:
            f.write("gas_name\tgctop\tgcbot\tftopg\tfbotg\n")
            for i, key in enumerate(self.gases.keys()):
                f.write(f'"{key}"\t'
                        f"{gas_conc_top_bc[i]}\t"
                        f"{gas_conc_bot_bc[i]}\t"
                        f"{gas_flux_top_bc[i]}\t"
                        f"{gas_flux_bot_bc[i]}\n")


        ## handle cloud bcs
        elem_top_conc_bc = np.zeros((self.NBIN, len(self.elems)))
        elem_bot_conc_bc = np.zeros((self.NBIN, len(self.elems)))
        elem_top_flux_bc = np.zeros((self.NBIN, len(self.elems)))
        elem_bot_flux_bc = np.zeros((self.NBIN, len(self.elems)))
        for i, key in enumerate(self.elems.keys()):
            e = self.elems[key]
            g = e.group
            elem_top_conc_bc[:, i] = g.boundary.get("top_conc", np.zeros(self.NBIN))
            elem_top_flux_bc[:, i] = g.boundary.get("top_flux", np.zeros(self.NBIN))
            elem_bot_conc_bc[:, i] = g.boundary.get("bot_conc", np.zeros(self.NBIN))
            elem_bot_flux_bc[:, i] = g.boundary.get("bot_flux", np.zeros(self.NBIN))


        with open(os.path.join(path, io["p_boundary_file"]), "w+") as f:
            f.write("elem_name\tibin\tpctop\tpcbot\tftopp\tfbotp\n")
            for ibin in range(self.NBIN):
                for j, key in enumerate(self.elems.keys()):
                    f.write(f'"{key}"\t'
                            f"{ibin}\t"
                            f"{elem_top_conc_bc[i, j]}\t"
                            f"{elem_bot_conc_bc[i, j]}\t"
                            f"{elem_top_flux_bc[i, j]}\t"
                            f"{elem_bot_flux_bc[i, j]}\n")

        with _cd(path):

            try:
                subprocess.run(["export", "OMP_NUM_THREADS=1"], 
                               shell=True,
                               stdout=subprocess.PIPE)
                
                subprocess.run(["export", "KMP_STACKSIZE=128M"],
                                shell=True,
                                stdout=subprocess.PIPE)
                p = subprocess.Popen(
                    os.path.join(SRC, "carmapy", "carmapy.exe"), 
                    shell=False, 
                    stdout=subprocess.PIPE)
                
                while p.poll() is None:
                    l = p.stdout.readline() #blocks until it receives a newline.
                    if not suppress_output: print(l.decode('UTF-8'))
                # When the subprocess ends there might be unconsumed output 
                # that still needs to be processed.
                if not suppress_output: print(p.stdout.read().decode('UTF-8'))
            except Exception as e:
                print(e)
            
    def read_results(self, read_diag=False) -> None:
        """Reads in results of the carma simulation.  

        Parameters
        ----------
        If true reads in the microphysical rates and core mass fraction.
        Defaults to False.

        See Also
        --------
        carmapy.Results()
        """
        self.results = Results(self, read_diag=read_diag)
        
        
    def _cite(self) -> str:
        """Generates the required citations based on the modules used

        Returns
        -------
        str
            A reccomended citation and a list of references to cite
        """

        if self.is_2d:
            self._citation["2d"] = True

        cites = ("This work made use of CARMApy (Cukier et al. in prep), "
                "a wrapper of the ExoCARMA (Gao and Benneke 2018, Powell et al."
                " 2018) branch of CARMA (Turco et al. 1979, Toon et al. 1988) "
                "version 3.0 (Bardeen et al. 2008, 2010).  "
        )

        if self._citation.get("2d", False):   
            cites += ("2-D CARMA was implemented by Powell and Zhang 2024.  "
            )

        if self._citation.get("fastchem", False):
            cites += ("Atmospheric chemistry calculations were performed using "
            "fastchem (Stock et al. 2018, 2022).  ")

        if self._citation.get("picaso", False):
            cites += ("Spectra were generated using picaso "
            "(Batalha et al. 2019).  ")
    
        print(cites)

        return cites
    
class Element:
    """An object representing a single species in a Group.  Note that element
    refers not to chemical elements, but to components of the group (ie TiO2 
    Core or MgSiO2 Mantle).

    Parameters
    ----------
    name : str
        The name of the element
    ielem : int
        The index of the element in the Carma simulation
    group : Group
        The group the element is part of
    rho : float
        The density of the element
    proc : str
        Flag that describes the element type.  Currently supported types are:

        - "Volatile": Homogeneously nucleating elements or the outer element in a
          heterogeneously nucleating group

        - "Core Mass": The inner element in a heterogeneously nucleating group
    igas : int
        The index of the gas which the element grows from and evaporates to in t
        he Carma simulation
    """

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
        self.igas:  int     = igas #TODO chance to reference gas directly
    

_DEFAULT_GAS_BC = {"bot_conc": -1,
                    "top_conc": -1,
                    "bot_flux": 0,
                    "top_flux": 0}

class Gas: 
    """An object representing a limiting gas resevoir for a condensate.  Attributes
    not defined in ``**kwargs`` will be populated from carmapy defaults, if available

    Parameters
    ----------
    name : str
        The name of the gas resevoir.  Note currently we name the gas resevoir
        with the name of the condensate, not the gas phase (ie Mg2SiO4 not Mg) 
    igas : int
        The index of the gas in the Carma simulation
    nmr : ArrayLike, optional
        The number mixing ratio.  If a float is provided represents the number
        mixing ratio at the bottom of the atmosphere; if an array then 
        represents the mixing ratio at each "center" location.  Can be provided
        at anytime before simulation run, by default -1 (not initialized)


    Notes
    -----
    1. The vapor pressure of the condensate is calculated as follows:

        ``vp = 1e6 * 10**(offset - vp_tcoeff/T - vp_metcoeff * met - vp_logp_oeff * log10(P*1e-6))``

        with ``vp`` in baryes, ``T`` in K, and ``P`` in baryes.

    
    2. The surface tension of the condensate is calculated as follows:

        ``surften = surften_0 - surften_slope * (T)``

        with ``T`` in K and ``surften`` in dyne / cm

    3. If not directly given, the latent heat of evaportation is calculated
       folowing Charnay et al. 2015 [2]_

       ``lat_heat_e = vp_coeff * log(10) * R/wtmol_dif``

       where R is the ideal gas constant

    4. The ``gas.boundary`` dictionary is structured as follows:
        - ["bot_conc"] describes the number mixing ratio of the gas at the base
          of the atmosphere (only used if the bottom gas boundary condition is
          set to "fixed_conc").  If not set, uses the value of `gas.nmr` at the
          bottom of the atmosphere instead.
        - ["top_conc"] describes the number mixing ratio of the gas at the top
          of the atmosphere (only used if the top gas boundary condition is
          set to "fixed_conc").  Will throw an error if used but not set.
        - ["bot_flux"] describes the upwards flux of the gas to the base
          of the atmosphere (only used if the bottom gas boundary condition is
          set to "fixed_flux") [g/cm^2/s].  Defaults to 0 if not set.
        - ["top_flux"] describes the downwards flux of the gas at the top
          of the atmosphere (only used if the top gas boundary condition is
          set to "fixed_flux"). [g/cm^2/s] Defaults to 0 if not set.


    References
    ----------
    .. [1] Helling, C., & Woitke, P. 2006, Astronomy and Astrophysics, 
       Volume 455, Issue 1, August III 2006, pp325-338, 455, 325

    .. [2] Charnay et al. 2015, ApJL 813, L1
    """
    

    wtmol: float 
    """ Molar mass of the condensate formed by the gas [g/mol] """

    gcomp: int = 0
    """ Integer that indicates composition (1 Water, 2 H2SO4, 0 other) """

    wtmol_dif: float 
    """ Molar mass of the gas phase of the gas [g/mol] """

    rho_cond: float 
    """ Density of the condensate formed by the gas [g/cm³]"""

    surften_0: float 
    """ Surface tension at 0 K assuming linear trend holds (see note 2) [dyne/cm] """

    surften_slope: float = 0
    """ Slope of surface tension with temperature (see note 2) [dyne/cm/K]"""

    coldia: float 
    """ Collisional diameter of the condensate [cm] """

    vp_offset: float 
    """ Constant term in vapor pressure equation (see note 1)"""

    vp_tcoeff: float 
    """ Coeficcient to temperature term in vapor pressure equation (see note 1) [K]"""

    vp_metcoeff: float = 0
    """ Coeficcient to metallicity term in vapor pressure equation (see note 1)"""

    vp_logpcoeff: float = 0
    """ Coeficcent to pressure term in vapor pressure equation (see note 1)"""

    is_typeIII: bool = False
    """ True if condensation reaction is a Type III reaction (see Helling & Woitke 2006) [1]_"""

    lat_heat_e: float = -1
    """ Latent heat of evaporation, if not provided derived from other inputs (see note 3) """

    stofact: int
    """ The stoichiometry factor between the gas phase and the condensate """
    
    hill_formula: str
    """ The chemical formula of the condensate in hill notation"""

    boundary: dict
    """ The boundary conditions for the gas (see note 4) """

    def __init__(self, 
                 name: str, 
                 igas: int, 
                 nmr: ArrayLike = -1,
                 **kwargs):
        self.name: str = name
        self.igas: int = igas        
        self.nmr: ArrayLike = nmr

        defaults = gas_dict.get(name, kwargs)

        self.wtmol          = kwargs.get("wtmol",         defaults["wtmol"])
        self.wtmol_dif      = kwargs.get("wtmol_dif",     defaults["wtmol_dif"])
        self.gcomp          = kwargs.get("gcomp",         defaults.get("gcomp", 0))
        self.rho_cond       = kwargs.get("rho_cond",      defaults["rho_cond"])
        self.surften_0      = kwargs.get("surften_0",     defaults["surften_0"])
        self.surften_slope  = kwargs.get("surften_slope", defaults["surften_slope"])
        self.coldia         = kwargs.get("coldia",        defaults["coldia"])
        self.vp_offset      = kwargs.get("vp_offset",     defaults["vp_offset"])
        self.vp_tcoeff      = kwargs.get("vp_tcoeff",     defaults["vp_tcoeff"])
        self.vp_metcoeff    = kwargs.get("vp_metcoeff",   defaults["vp_metcoeff"])
        self.vp_logpcoeff   = kwargs.get("vp_logpcoeff",  defaults["vp_logpcoeff"])
        self.is_typeIII     = kwargs.get("is_typeIII",    defaults["is_typeIII"])
        self.stofact        = kwargs.get("stofact",       defaults["stofact"])
        self.lat_heat_e     = kwargs.get("lat_heat_e",    defaults.get("lat_heat_e", -1))
        self.hill_formula   = kwargs.get("hill_formula",  defaults["hill_formula"])
        self.boundary       = kwargs.get("boundary",      _DEFAULT_GAS_BC)

        
    
class Nuc:
    """An object representing a nucleation pathway.

    Parameters
    ----------
    group_core : Group
        The group forms the seed particle for the nucleation event.  If 
        homogeneous nucleation, should be the same as ``group_mantle``
    group_mantle : _type_
        The group which forms as a result of the nucleation
    is_het : bool
        if True, the nucleation is heterogeneous.  If False, homogeneous
    gas : Gas
        The gas resevoir which condenses
    mucos : float
        The cosine of the contact angle between the gas and the seed particle. 
        Ignored if ``is_het==False``.
    
    """
    def __init__(self, 
                 group_core: "Group", 
                 group_mantle: "Group",
                 is_het: bool,
                 gas : "Gas",
                 mucos: "float") -> None:
        self.group_core = group_core
        self.group_mantle = group_mantle
        self.is_het = is_het
        self.gas = gas
        self.ievp2elem = group_core.core
        self.mucos = mucos
    
class Growth:
    """An object representing a growth/evaportation pathway.

    Parameters
    ----------
    elem : Element
        The element upon which the gas condenses/evaportates from
    gas : Gas
        The gas resevoir from/to which the gas condenses/evaportates

    """
    def __init__(self, elem: "Element", gas: "Gas") -> None:
        self.elem = elem
        self.gas = gas
    
class Group:
    """An object representing a cloud species
        

    Parameters
    ----------
    igroup : int
        The index of the group in the Carma simulation
    name : str
        A name for the group.  The typical naming convention is "Pure <Species>"
        (ie "Pure TiO2") for homogeneously nucleating groups and "<Mantle> on 
        <Core>" (ie "Mg2SiO4 on TiO2") for heterogeneously nucleating groups
    rmin : float
        The minimum size of the condensate

    Notes
    -----
    1. The ``group.boundary`` describes the boundary conditions of the 
    particulate matter in the atmosphere.  Currently, CARMApy only supports 
    non-zero boundary conditions for single element groups (ex. homogeneously
    nucleated groups).  Noting that each of these entries is either 0 or an 
    array of NBIN elements, the structure of the dictionary is as follows:
    
    - ["bot_conc"] describes the concentration of the group at the base
       of the atmosphere (only used if the bottom cloud boundary condition is
       set to "fixed_conc"). [particles/cm^3]  If not set, defaults to 0.
    - ["top_conc"] describes the concentration of the group at the top
       of the atmosphere (only used if the top cloud boundary condition is
       set to "fixed_conc").  [particles/cm^3] If not set, defaults to 0.
    - ["bot_flux"] describes the upwards flux of the group to the base
       of the atmosphere (only used if the bottom cloud boundary condition is
       set to "fixed_flux") [particles/cm^2/s].  Defaults to 0 if not set.
    - ["top_flux"] describes the downwards flux of the group at the top
       of the atmosphere (only used if the top cloud boundary condition is
       set to "fixed_flux"). [particles/cm^2/s] Defaults to 0 if not set.

    """

    core : "Element"
    """ The element which represents the original seed particle. """

    mantle: "Element"
    """ The element on the surface of the cloud particle. """

    boundary: dict
    """ The boundary conditions for the group (see note 1). """


    def __init__(self, igroup: int, name: str, rmin: float) -> None:
        self.igroup = igroup
        self.name = name
        self.rmin = rmin
        self.core = None
        self.mantle = None
        self.boundary = {}
    
    def coreify(self, ielem: int, group: "Group", gas_name: str) -> Element:
        """Create a core element from the only element of the current group.
        Used to create the core element of a heterogeneously nucleating group
        where this group serves as the seed particles.  Adds the created element
        as the core of the provided group

        Parameters
        ----------
        ielem : int
            The index of the new element in the Carma simulation
        group : Group
            The group to which the new element belongs
        gas_name : str
            Name of the gas resevoir corresponding to the new element

        Returns
        -------
        Element
            The newly created element initialized as the core of a heterogeneous
            group
        """
        core_elem = self.core
        
        name = core_elem.name
        name = name.split(" ")[-1]
        name = name + f" Core ({gas_name})"
        
        elem = Element(name, ielem, group, core_elem.rho, "Core Mass", core_elem.igas)
        group.core = elem
        return elem
    



class Coag:
    """An object which represents a coagulation pathway of a condensate with 
    itself.

    Parameters
    ----------
    group : Group
        The group which represents the condensate which is able to coagulate
    """
    def __init__(self, group: "Group") -> None:
        self.group = group

