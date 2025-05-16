import pyfastchem
from .carmapy import Carma
from .constants import *
import numpy as np


def get_fastchem_abundances(T : np.ndarray, 
                            P : np.ndarray, 
                            species : list,
                            metallicity : float = 1):
    
  temperature = T
  pressure = np.array(P) / BAR_TO_BARYE

  fastchem = pyfastchem.FastChem(
    'inputs/fastchem/asplund_2009_extended.dat',
    'inputs/fastchem/logK.dat',
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
  
  print(number_densities.shape)
  ret = []
  for s in species:
    index = fastchem.getGasSpeciesIndex(s)
    if index == pyfastchem.FASTCHEM_UNKNOWN_SPECIES:
      raise ValueError(f"{s} is an unknown species")
    ret.append(nmr[:, index])
  
  
  return np.array(ret)

def populate_fastchem_abundances(carma: Carma, metalicity = 1.0, override = {"H2O": 0}):
  species = []
  
  
  for gas in carma.gasses.keys():
      s = fastchem_species.get(gas, -1) 
      
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
  
        
    