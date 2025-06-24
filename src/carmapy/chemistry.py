import pyfastchem
from carmapy.constants import *
import numpy as np
from scipy.interpolate import interp1d

from numpy.typing import ArrayLike
import os
SRC = os.path.dirname(os.path.dirname(__file__))


def get_fastchem_abundances(T : np.ndarray, 
                            P : np.ndarray, 
                            species : list,
                            metallicity : float = 1):
    
  temperature = T
  pressure = np.array(P) / BAR_TO_BARYE

  fastchem = pyfastchem.FastChem(
    os.path.join(SRC, "fastchem", "asplund_2009_extended.dat"),
    os.path.join(SRC, "fastchem", "logK.dat"),
    1)


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


  print("FastChem reports:")
  print("  -", pyfastchem.FASTCHEM_MSG[fastchem_flag])

  if np.amin(output_data.element_conserved[:]) == 1:
    print("  - element conservation: ok")
  else:
    print("  - element conservation: fail")
    
  number_densities = np.array(output_data.number_densities)

  nmr = number_densities / np.repeat((P/(k_B * T))[:, np.newaxis], number_densities.shape[1], axis=1)
  
  ret = []
  for s in species:
    index = fastchem.getGasSpeciesIndex(s)
    if index == pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
      raise ValueError(f"{s} is an unknown species")
    ret.append(nmr[:, index])
  
  
  return np.array(ret)

def condensation_curve(P: float | ArrayLike, 
                       T: float | ArrayLike, 
                       met: float,
                       gas: str) -> float | ArrayLike:

  offset     = gas_dict[gas]["vaprtn"].get("offset", 0)
  T_coeff    = gas_dict[gas]["vaprtn"].get("T_coeff", 0)
  met_coeff  = gas_dict[gas]["vaprtn"].get("met_coeff", 0)
  logp_coeff = gas_dict[gas]["vaprtn"].get("logp_coeff", 0)

  return 10**(offset
              - T_coeff/T
              - met_coeff * met
              - logp_coeff * np.log10(P))

def populate_fastchem_abundances(carma: "Carma", metalicity = 1.0, override = {"H2O": 0}):
  species = []
  
  
  for gas in carma.gasses.keys():
      s = gas_dict[gas].get("fastchem_specie", -1) 
      
      if s == -1:
          raise ValueError(f"{gas} is not currently supported by the carmapy fastchem interface")
      species.append(s)
        
  abund = get_fastchem_abundances(carma.T_centers, carma.P_centers, species)
  
  
  nmr_dict = {}
  
  for i in range(len(carma.gasses.keys())):
    nmr_dict[list(carma.gasses.keys())[i]] = abund[i, 0]
  
  for key in override.keys():
    nmr_dict[key] = override[key]
  
  carma.set_nmr(nmr_dict)

def populate_abundances_at_cloud_base(carma, species, metalicity):
  P = carma.P_levels
  T = carma.T_levels

  p_t = interp1d(T, P)

  override= {"H2O": 0}

  for s in species:
    cond_curve = interp1d(T, condensation_curve(P, T, np.log10(metalicity), s))
    Ts = np.linspace(np.min(T), np.max(T), 50000)

    intersection = np.argmin(np.abs(p_t(Ts) - cond_curve(Ts)))

    fast_chem_gas = fastchem_species.get(s, -1) 
    if fast_chem_gas == -1: raise ValueError("{s} is not currently supported by the carmapy fastchem interface")

    override[s] = get_fastchem_abundances(np.array([p_t(Ts[intersection])]), np.array([Ts[intersection]]), [fast_chem_gas], metalicity)[0]
  
  populate_fastchem_abundances(carma, metalicity, override)
  pass
  
        
    