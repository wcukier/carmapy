""" functions that interface with fastchem """

import pyfastchem
from carmapy.constants import *
import numpy as np
from scipy.interpolate import interp1d
import scipy.optimize 
from numpy.typing import ArrayLike
import os
SRC = os.path.dirname(os.path.dirname(__file__))
import io
import sys

def get_fastchem_abundances(T : np.ndarray, 
                            P : np.ndarray, 
                            species : list,
                            metallicity : float = 1):
  """Uses fastchem to calculate the abundances of the provided species.

  Parameters
  ----------
  T : np.ndarray
      Temperature profile of the atmosphere [K]
  P : np.ndarray
      Pressure profile of the atmosphere [barye]
  species : list
      List of 
  metallicity : float, optional
      _description_, by default 1

  Returns
  -------
  _type_
      _description_

  Raises
  ------
  RuntimeError
      _description_
  ValueError
      _description_
  """
    
  temperature = T
  pressure = np.array(P) / BAR_TO_BARYE

  fastchem = pyfastchem.FastChem(
    os.path.join(SRC, "fastchem", "asplund_2009_extended.dat"),
    os.path.join(SRC, "fastchem", "logK.dat"),
    0)

  input_data = pyfastchem.FastChemInput()
  output_data = pyfastchem.FastChemOutput()


  input_data.temperature = temperature
  input_data.pressure = pressure
  fastchem_flag = fastchem.calcDensities(input_data, output_data)

  solar_abundances = np.array(fastchem.getElementAbundances())

  element_abundances = np.copy(solar_abundances)
    
    #scale the element abundances, except those of H and He
  for j in range(0, fastchem.getElementNumber()):
    if fastchem.getElementSymbol(j) != 'H' and fastchem.getElementSymbol(j) != 'He':
      element_abundances[j] *= metallicity
      
  fastchem.setElementAbundances(element_abundances)

  if np.amin(output_data.element_conserved[:]) == 0:
    raise RuntimeError("FastChem - element conservation: fail")

    
  number_densities = np.array(output_data.number_densities)

  nmr = number_densities / np.repeat((P/(k_B * T))[:, np.newaxis], number_densities.shape[1], axis=1)
  
  ret = []
  for s in species:
    index = fastchem.getGasSpeciesIndex(s)
    if index == pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
      raise ValueError(f"{s} is an unknown species")
    ret.append(nmr[:, index])
  
  
  return np.array(ret)

def saturation_vapor_pressure(P: float | ArrayLike, 
                       T: float | ArrayLike, 
                       log_met: float,
                       gas: str) -> float | ArrayLike:

  offset     = gas_dict[gas].get("vp_offset", 0)
  T_coeff    = gas_dict[gas].get("vp_tcoeff", 0)
  met_coeff  = gas_dict[gas].get("vp_metcoeff", 0)
  logp_coeff = gas_dict[gas].get("vp_logpcoeff", 0)

  return 10**(offset
              - T_coeff/T
              - met_coeff * log_met
              - logp_coeff * np.log10(P))

def populate_fastchem_abundances(carma: "Carma", metallicity=1, override = {"H2O": 0}):
  species = []
  
  
  for gas in carma.gasses.keys():
      s = gas_dict[gas].get("fastchem_species", -1) 
      
      if s == -1:
          raise ValueError(f"{gas} is not currently supported by the carmapy fastchem interface")
      species.append(s)
        
  abund = get_fastchem_abundances(carma.T_centers, carma.P_centers, species, metallicity)
  
  
  nmr_dict = {}
  
  for i in range(len(carma.gasses.keys())):
    nmr_dict[list(carma.gasses.keys())[i]] = abund[i, 0]
  
  for key in override.keys():
    nmr_dict[key] = override[key]
  
  carma.set_nmr(nmr_dict)


def find_cloud_base(carma, species, metallicity=1):

  P = carma.P_centers
  T = carma.T_centers


  sat_vp = saturation_vapor_pressure(P, T, np.log10(metallicity), species)
  abund = get_fastchem_abundances(T, P, [gas_dict[species]["fastchem_species"]], metallicity)[0, :]

  i = 0
  while(sat_vp[i]/P[i] > abund[i]): 
    i += 1
  
  guess = T[i]

  p_t = interp1d(T, P)

  def _diff(T):

    sat_vp = saturation_vapor_pressure(p_t(T), T, np.log10(metallicity), species)/p_t(T)
    abund = get_fastchem_abundances(np.array(T), np.array(p_t(T)), [gas_dict[species]["fastchem_species"]], metallicity)[0, :]

    return abund- sat_vp
  

  root = scipy.optimize.root(_diff, guess).x[0]
  return p_t(root), root

def populate_abundances_at_cloud_base(carma, metallicity=1):
  P = carma.P_centers
  T = carma.T_centers

  p_t = interp1d(T, P)

  override= {"H2O": 0}

  for s in list(carma.gasses.keys())[1:]:

    P_int, T_int = find_cloud_base(carma, s, metallicity)

    fast_chem_gas = gas_dict[s].get("fastchem_species", -1) 
    if fast_chem_gas == -1: raise ValueError("{s} is not currently supported by the carmapy fastchem interface")

    override[s] = get_fastchem_abundances(np.array([T_int]), np.array([P_int]), [fast_chem_gas], metallicity)[0]
  
  print(override)
  populate_fastchem_abundances(carma, metallicity, override)
  
  
        
    