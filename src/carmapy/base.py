"""Essential functions for carmapy

Note
----
All of these fuctions are loaded into the top-level carmapy namespace.  In other
words you can access them like ``carmapy.load_carma(...)`` instead of 
``carmapy.base.load_carma(...)``"""

from carmapy.constants import *
from carmapy.results import *
from carmapy.carmapy import *
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

SRC = os.path.dirname(os.path.dirname(__file__))

_int2bc = {
    I_ZERO_CGRAD: "zero_grad",
    I_FIXED_CONC: "fixed_conc",
    I_FLUX_SPEC: "fixed_flux"
}

def load_carma(path: str, restart: int =0) -> Carma:
    """Restores a carma object from files.  Loads in the configuration stored by
    ``Carma.run()``

    Parameters
    ----------
    path : str
        The path to the directory holding the carma files
    restart : int, optional
        If 1, future calls to ``Carma.run()`` for this object will resume the 
        simulation from the most recently saved state.  If 0, instead calls to 
        ``Carma.run()`` will restart the simulation from the beginning.  Either
        way, the current output files will be overwritten by default. By default
        0

    Returns
    -------
    Carma
        The loaded Carma object
    """

    carma = Carma(path)
    carma.restart = restart
    
    nml = f90nml.read(os.path.join(path, "inputs", "input.nml"))

    carma.NZ = nml["input_params"]["NZ"] 
    carma.NLONGITUDE = nml["input_params"]["NLONGITUDE"]
    carma.is_2d = nml["input_params"]["IS_2D"]
    carma.igridv - nml["input_params"]["igridv"]
    carma.NBIN = nml["input_params"]["NBIN"]
    carma.output_gap = nml["input_params"]["iskip"]
    carma.n_tstep  = nml["input_params"]["nstep"]
    carma.dt = nml["input_params"]["dtime"]
    carma.top_bound_type_cloud = _int2bc[nml["input_params"]["itbnd_pc"]]
    carma.bot_bound_type_cloud = _int2bc[nml["input_params"]["ibbnd_pc"]]
    carma.top_bound_type_gas = _int2bc[nml["input_params"]["itbnd_gc"]]
    carma.bot_bound_type_gas = _int2bc[nml["input_params"]["ibbnd_gc"]]

    carma.wt_mol = nml["physical_params"]["wtmol_air_set"]
    carma.surface_grav = nml["physical_params"]["grav_set"]
    carma.r_planet = nml["physical_params"]["rplanet"]
    carma.velocity_avg = nml["physical_params"]["velocity_avg"]
    carma.log_metallicity = nml["physical_params"]["met"]

    carma.atmo = {
        "rmu_0": nml["physical_params"]["rmu_0"],
        "rmu_t0": nml["physical_params"]["rmu_t0"],
        "rmu_c": nml["physical_params"]["rmu_c"],
        "thcond_0": nml["physical_params"]["thcond_0"],
        "thcond_1": nml["physical_params"]["thcond_1"],
        "thcond_2": nml["physical_params"]["thcond_2"],
        "CP" : nml["physical_params"]["CP"]
    }

    io = nml["io_files"]

    with open(os.path.join(path, io["groups_file"]), "r") as f:
        f.readline()
        for line in f:
            name, rmin = shlex.split(line[:-1])
            carma.groups[name] = Group(len(carma.groups)+1, name, float(rmin))
        
    with open(os.path.join(path, io["gases_file"]), "r") as f:
        f.readline()
        for line in f:
            (name, 
             wtmol, 
             _, 
             gcomp, 
             wtmol_dif,
             rho_cond,
             surften_0,
             coldia,
             vp_offset,
             vp_tcoeff,
             is_typeIII,
             surften_slope,
             vp_metcoeff,
             vp_logpcoeff,
             lat_heat_e,
             stofact) = shlex.split(line[:-1])
            
            name= name[:-len(' Vapor')]
            carma.gases[name] = Gas(name, 
                                     len(carma.gases)+1, 
                                     wtmol=float(wtmol), 
                                     wtmol_dif=float(wtmol_dif),
                                     gcomp = int(gcomp),
                                     rho_cond = float(rho_cond),
                                     surften_0 = float(surften_0),
                                     surften_slope = float(surften_slope),
                                     coldia = float(coldia),
                                     vp_offset = float(vp_offset),
                                     vp_tcoeff = float(vp_tcoeff),
                                     vp_metcoeff = float(vp_metcoeff),
                                     vp_logpcoeff = float(vp_logpcoeff),
                                     is_typeIII = bool(is_typeIII),
                                     lat_heat_e=float(lat_heat_e),
                                     stofact=int(stofact)
                                     )
            
    
    with open(os.path.join(path, io["elements_file"]), "r") as f:
        f.readline()
        for line in f:
            igroup, name, rho, proc, igas = shlex.split(line[:-1])
            group = carma.groups[list(carma.groups.keys())[int(igroup)-1]]
            carma.elems[name] = Element(name, 
                                        len(carma.elems)+1, 
                                        group, 
                                        float(rho), 
                                        proc, 
                                        int(igas))
            
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
            group_core = carma.elems[list(carma.elems.keys())[int(ele_from)-1]].group
            group_mantle = carma.elems[list(carma.elems.keys())[int(ele_to)-1]].group
            gas = carma.gases[list(carma.gases.keys())[int(igas)-1]]
            carma.nucs.append(Nuc(group_core, group_mantle, is_het, gas, float(mucos)))
            
    with open(os.path.join(path, io["growth_file"])) as f:
        f.readline()
        for line in f:
            ielem, igas = shlex.split(line[:-1])
            elem = carma.elems[list(carma.elems.keys())[int(ielem)-1]]
            gas = carma.gases[list(carma.gases.keys())[int(igas)-1]]
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
    carma.T_levels = np.genfromtxt(os.path.join(path, "inputs", "temp_levels.txt"))


    gas_input = np.genfromtxt(os.path.join(path, io["gas_input_file"]))
    for i, key in enumerate(carma.gases.keys()):
        carma.gases[key].nmr = gas_input[1:, i]
        
    winds = np.zeros(carma.NZ)
    with open(os.path.join(path, io["winds_file"])) as f:
        for i in range(carma.NZ):
            winds[i] = float(f.readline())

    with open(os.path.join(path, io["g_boundary_file"])) as f:
        f.readline()
        for key in carma.gases.keys():
            g = carma.gases[key]
            _, gctop, gcbot, ftopg, fbotg = shlex.split(f.readline()[:-1])

            g.boundary = {
                "top_conc": float(gctop),
                "bot_conc": float(gcbot),
                "top_flux": float(ftopg),
                "bot_flux": float(fbotg)
            }

    elem_top_conc_bc = np.zeros((carma.NBIN, len(carma.elems)))
    elem_bot_conc_bc = np.zeros((carma.NBIN, len(carma.elems)))
    elem_top_flux_bc = np.zeros((carma.NBIN, len(carma.elems)))
    elem_bot_flux_bc = np.zeros((carma.NBIN, len(carma.elems)))

    with open(os.path.join(path, io["p_boundary_file"])) as f:
        f.readline()
        for ibin in range(carma.NBIN):
            for j in range(len(carma.elems.keys())):
                line = f.readline()
                _, _, pctop, pcbot, ftopp, fbotp = shlex.split(line)

                elem_top_conc_bc[ibin, j] = float(pctop)
                elem_bot_conc_bc[ibin, j] = float(pcbot)
                elem_top_flux_bc[ibin, j] = float(ftopp)
                elem_bot_flux_bc[ibin, j] = float(fbotp)

    for j, key in enumerate(carma.elems.keys()):
        e = carma.elems[key]
        e.group.boundary["top_conc"] = elem_top_conc_bc[:, j]
        e.group.boundary["bot_conc"] = elem_bot_conc_bc[:, j]
        e.group.boundary["top_flux"] = elem_top_flux_bc[:, j]
        e.group.boundary["bot_flux"] = elem_bot_flux_bc[:, j]

    return carma


def available_species() -> None:
    """Prints the gas species available in the carmapy defaults
    """
    print(list(gas_dict.keys())[1:])


def included_mucos(species):
    """Prints the cosines of the contanct angle available in carmapy defaults.
    Prints the dictionary of potential seed particles and their contact angles
    for the requested species from the carmapy defaults.

    Parameters
    ----------
    species : str
        The name of the gas species for which to see the available contact 
        angles
    """
    print(gas_dict[species]["mucos_dict"])