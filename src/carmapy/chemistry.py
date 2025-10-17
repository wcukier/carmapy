""" Functions to determine the chemical properties of the atmosphere """

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
from typing import Union


def get_fastchem_abundances(T : np.ndarray, 
                            P : np.ndarray, 
                            species : list[str],
                            metallicity : float = 1) -> np.ndarray:
  """Uses fastchem to calculate the gas phase abundances of the provided 
  species at the provided T-P points

  Parameters
  ----------
  T : np.ndarray
      Temperature profile of the atmosphere [K]
  P : np.ndarray
      Pressure profile of the atmosphere [barye]
  species : list[str]
      List of species in Hill Notation (so TiO2 would be "O2Ti1")
  metallicity : float, optional
      metallicity relative to solar (not log), by default 1

  Returns
  -------
  np.ndarray
      A 1-D (if only 1 species provided) or 2-D (if more than 1 species 
      provided) array of the number mixing ratio of each species at each 
      P-T point

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
    if ((fastchem.getElementSymbol(j) != 'H')
         and (fastchem.getElementSymbol(j) != 'He')):
      element_abundances[j] *= metallicity
      
  fastchem.setElementAbundances(element_abundances)

  if np.amin(output_data.element_conserved[:]) == 0:
    raise RuntimeError("FastChem - element conservation: fail")

    
  number_densities = np.array(output_data.number_densities)

  nmr = number_densities / np.repeat((P/(k_B * T))[:, np.newaxis],
                                      number_densities.shape[1], 
                                      axis=1)
  
  ret = []
  for s in species:
    index = fastchem.getGasSpeciesIndex(s)
    if index == pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
      raise ValueError(f"{s} is an unknown species")
    ret.append(nmr[:, index])
  
  
  return np.array(ret)

def saturation_vapor_pressure(P: ArrayLike, 
                       T:  ArrayLike, 
                       log_met: float,
                       gas: Union[str, "Gas"]) -> ArrayLike:
  """Calculates the saturation vapor pressure for carmapy default gases.

  Parameters
  ----------
  P : ArrayLike
      Atmospheric Pressure [barye]
  T : ArrayLike
      Temperature [K]
  log_met : float
      Log solar metallicity
  gas : Union[str, Gas]
      carma Gas object or name of the carmapy default gas

  Returns
  -------
  ArrayLike
      The saturation vapor pressure of the gas at each requested P-T point
  """

  if type(gas) == type(""):
    offset     = gas_dict[gas].get("vp_offset", 0)
    T_coeff    = gas_dict[gas].get("vp_tcoeff", 0)
    met_coeff  = gas_dict[gas].get("vp_metcoeff", 0)
    logp_coeff = gas_dict[gas].get("vp_logpcoeff", 0)
  else:
    offset     = gas.vp_offset
    T_coeff    = gas.vp_tcoeff
    met_coeff  = gas.vp_metcoeff
    logp_coeff = gas.vp_logpcoeff

  return 1e6 * 10**(offset
              - T_coeff/T
              - met_coeff * log_met
              - logp_coeff * np.log10(1e-6*P))

def _populate_fastchem_abundances(carma: "Carma",
                                  metallicity=1, 
                                  override = {"H2O": 0}) -> None:
 
  species = []
  
  
  for gas in carma.gases.keys():
      s = gas_dict[gas].get("hill_formula", -1) 
      
      if s == -1:
          raise ValueError(f"{gas} is not currently supported by the carmapy "
                           "fastchem interface")
      species.append(s)
        
  abund = get_fastchem_abundances(carma.T_centers, 
                                  carma.P_centers, 
                                  species, 
                                  metallicity)
  
  
  nmr_dict = {}
  
  for i in range(len(carma.gases.keys())):
    nmr_dict[list(carma.gases.keys())[i]] = abund[i, 0]
  
  for key in override.keys():
    nmr_dict[key] = override[key]
  
  carma.set_nmr(nmr_dict)


def find_cloud_base(carma: "Carma", species: str) -> tuple[float, float]:
  """Locates the P-T coordianates of the cloud base.

  Parameters
  ----------
  carma : Carma
      A carma object with initialized gases, P-T structure, and log metallicity
  species : str
      Either the name of the carmapy default gas to find the cloud base of or 
      the hill notation chemical formula of the species

  Returns
  -------
  P : float
    The pressure at the cloud base [barye]
  T : float
    The temperature at the cloud base [K]
  """

  carma._citation["fastchem"] = True

  species = carma.gases[species]
  s = species.hill_formula

  P = carma.P_centers
  T = carma.T_centers

  # if carma.is_2d: T = np.mean(T, axis=1)
  if carma.is_2d: T = T[:, 0] # TODO


  metallicity = 10**carma.log_metallicity

  # print(s)

# TODO make T, P ordering consistent
  sat_vp = saturation_vapor_pressure(P, T, np.log10(metallicity), species) 
  abund = get_fastchem_abundances(T, P, [s], metallicity)[0, :]

  i = 0
  while(sat_vp[i]/P[i] > abund[i]): 
    i += 1
    # print(f"{i:2d}\t{sat_vp[i]/P[i]:.2e}\t{abund[i]:.2e}\t{T[i]:.0f}\t{P[i]:.2e}\t{sat_vp[i]:2e}")
    if (i == len(P)): 
      print(f"{s} does not condense")
      return P[i-1], T[i-1]
  
  guess = T[i]

  p_t = interp1d(T, P)

  def _diff(T):

    sat_vp = saturation_vapor_pressure(p_t(T),
                                       T, 
                                       np.log10(metallicity), 
                                       species) / p_t(T)
    
    abund = get_fastchem_abundances(np.array(T), 
                                    np.array(p_t(T)), 
                                    [s], 
                                    metallicity)[0, :]

    return abund - sat_vp
  

  root = scipy.optimize.root(_diff, guess).x[0]
  return p_t(root), root

def populate_abundances_at_cloud_base(carma: "Carma") -> None:
  """Set the number mixing ratios below the cloud base.  Uses the equilibrium
  concentration at the cloud base to determine what the mixing ratio should be

  Parameters
  ----------
  carma : Carma
      A carma object with initialized gases, P-T structure, and log metallicity

  """
  carma._citation["fastchem"] = True


  P = carma.P_centers
  T = carma.T_centers

  if carma.is_2d: T = np.mean(T, axis=1)

  p_t = interp1d(T, P)

  nmr_dict = {"H2O": 0}

  for s in list(carma.gases.keys())[1:]:

    P_int, T_int = find_cloud_base(carma, s)

    fast_chem_gas = carma.gases[s].hill_formula
    if fast_chem_gas == -1: 
      raise ValueError("{s} is not currently supported by "
                        "the carmapy fastchem interface")

    metallicity = 10 ** carma.log_metallicity

    nmr_dict[s] = get_fastchem_abundances(np.array([T_int]), 
                                          np.array([P_int]), 
                                          [fast_chem_gas], 
                                          metallicity)[0]
  
    carma.set_nmr(nmr_dict)

  
  
        
    