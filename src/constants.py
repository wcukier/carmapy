## This module is part of the CARMA module and contains enumerations that are part of
## the CARMA and CARMASTATE objects.
##
## @author Chuck Bardeen
## @ version July-2009

#--
# Index values of CARMA's flags.  In a given list, begin with 1
# (instead of 0) so that undefined flags will produce an error. 
#
# For example:
# if( itype(ielem) .eq. I_INVOLATILE )then
#
# If itype(ielem) hasn't been defined (and is still 0), we do not want
# to execute the statements that follow.

#  Define values of flag used for vertical transport
#  boundary conditions (ixxxbnd_pc)
I_FIXED_CONC = 1    ## Fixed Concentration
I_FLUX_SPEC  = 2    ## Flux Specification
I_ZERO_CGRAD = 3   ## Zero Concentration Gradient 

#  Define values of flag used for particle element
#  type specification (itype).
I_INVOLATILE = 1    ## Involatile particle
I_VOLATILE   = 2    ## Volatile particle
I_COREMASS   = 3    ## Core Mass
I_VOLCORE    = 4    ## Voltile Core
I_CORE2MOM   = 5    ## Core Mass - 2 Moments

##  Define values of flag used for nucleation process
##  specification (inucproc).
##
##  NOTE: Some of these can be used in combination, so for aerosol freezing this is treated
##  as a bit mask. When setting for one (or more) of the Aerosol freezing methods, use:
##    IAERFREEZE + I_AF_xxx + I_AF_yyy + ...
I_AF_TABAZADEH_2000 = 1     ## Aerosol Freezing, Tabazadeh[2000]
I_AF_KOOP_2000      = 2     ## Aerosol Freezing, Koop[2000]
I_AF_MOHLER_2010    = 4     ## Aerosol Freezing, Mohler[2010]
I_AF_MURRAY_2010    = 8     ## Glassy Aerosol Freezing, Murray[2010]
I_DROPACT           = 256   ## Droplet Activation
I_AERFREEZE         = 512   ## Aerosol Freezing
I_DROPFREEZE        = 1024  ## Droplet Freezing
I_ICEMELT           = 2048  ## Ice Melting
I_HETNUC            = 4096  ## Heterogeneous Nucleation
I_HOMNUC            = 8192  ## Binary homogeneous gas-to-particle nucleation
I_HOMGEN            = 128   ## General classic homogeneous nucleation
I_HETGEN            = 64    ## General classic heterogeneous nucleation

#  Define values of flag used for collection process (icollec)
I_COLLEC_CONST = 1   ## Constant Collection Efficiency
I_COLLEC_FUCHS = 2   ## Binwise Maxima of Fuchs' and Langmuir's Efficiencies
I_COLLEC_DATA  = 3   ## Input Data

#  Define values of flag used for coagulation operation (icoagop)
I_COAGOP_CONST = 1   ## Constant Coagulation Kernel
I_COAGOP_CALC  = 2   ## Calculate Coagulation Kernel

#  Define values of flag used for particle shape (ishape)
I_SPHERE   = 1   ## spherical
I_HEXAGON  = 2   ## hexagonal prisms or plates
I_CYLINDER = 3   ## circular disks, cylinders, or spheroids  
I_FRACTAL = 4   ## circular disks, cylinders, or spheroids

#  Define values of flag used for particle swelling parameterization (irhswell)
I_NO_SWELLING  = 0   ## No swelling
I_FITZGERALD   = 1   ## Fitzgerald
I_GERBER       = 2   ## Gerber
I_WTPCT_H2SO4  = 3   ## The weight percent method for sulfate aerosol

#  Define vallues of flag used for particle swelling composition (Fiztgerald)
I_SWF_NH42SO4   = 1   ## (NH4)2SO4
I_SWF_NH4NO3    = 2   ## NH4NO3
I_SWF_NANO3     = 3   ## NaNO3
I_SWF_NH4CL     = 4   ## NH4Cl
I_SWF_CACL2     = 5   ## CaCl2
I_SWF_NABR      = 6   ## NaBr
I_SWF_NACL      = 7   ## NaCl
I_SWF_MGCL2     = 8   ## MgCl2
I_SWF_LICL      = 9   ## LiCl

#  Define vallues of flag used for particle swelling composition (Gerber)
I_SWG_NH42SO4   = 11  ## (NH4)2SO4
I_SWG_SEA_SALT  = 12  ## Sea Salt
I_SWG_URBAN     = 13  ## Urban
I_SWG_RURAL     = 14  ## Rural

# Routines to calculate gas vapor pressures
I_VAPRTN_H2O_BUCK1981      = 1   ## H2O, Buck[1981]
I_VAPRTN_H2O_MURPHY2005    = 2   ## H2O, Murphy & Koop [2005]
I_VAPRTN_H2O_GOFF1946      = 3   ## H2O, Goff & Gratch [1946], used in CAM
I_VAPRTN_H2SO4_AYERS1980   = 4   ## H2SO4, Ayers [1980] & Kumala [1990]
I_VAPRTN_S8_FERREIRA2011   = 5   ## S8, Ferreira & Lobo [2011] 
I_VAPRTN_S2_LYONS2008      = 6  ## S2, Lyons [2008]
I_VAPRTN_KCL_MORLEY2012    = 7  ## KCl, Morley et al. [2012]
I_VAPRTN_ZNS_MORLEY2012    = 8  ## ZnS, Morley et al. [2012]
I_VAPRTN_NA2S_MORLEY2012      = 9  ## Na2S, Morley et al. [2012]
I_VAPRTN_MNS_MORLEY2012       = 10  ## MnS, Morley et al. [2012]
I_VAPRTN_CR_MORLEY2012        = 11  ## Cr, Morley et al. [2012]
I_VAPRTN_FE_VISSCHER2010      = 12  ## Fe, Visscher et al. [2010]
I_VAPRTN_MG2SIO4_VISSCHER2010 = 13  ## Mg2SiO4, Visscher et al. [2010]
I_VAPRTN_S8_ZAHNLE2016     = 14   ## S8, Zahnle et al. [2016] 
I_VAPRTN_TIO2_LODDERS1999     = 15   ## TiO2, Lodders [1999] 
I_VAPRTN_TIO2_HELLING2001     = 16   ## TiO2, Helling et al. [2001] 
I_VAPRTN_AL2O3_WAKEFORD2017     = 17   ## Al2O3, Wakeford et al. [2017] 
I_VAPRTN_CO_WYLIE1958     = 18   ## CO, Wylie thesis [1958] 

# Routines to calculate fall velocities
I_FALLRTN_STD              = 1   ## Standard CARMA 2.3 routine (spherical only)
I_FALLRTN_STD_SHAPE        = 2   ## Optional CARMA 2.3 routine (supports shapes)
I_FALLRTN_HEYMSFIELD2010   = 3   ## Heymsfield & Westbrook [2010] (ice only)

# Routines to calculate mie optical properties
I_MIERTN_TOON1981      = 1   ## Shell/Core, Toon & Ackerman [1981]
I_MIERTN_BOHREN1983    = 2   ## Homogeneous Sphere, Bohren & Huffman [1983]

# Gas Composition  
I_GCOMP_H2O             = 1   ## Water Vapor
I_GCOMP_H2SO4           = 2   ## Sulphuric Acid
I_GCOMP_SO2             = 3   ## Sulfer Dioxide
I_GCOMP_S8              = 4   ## Sulfur 8
I_GCOMP_S2              = 5   ## Sulfur 2
I_GCOMP_KCL             = 6   ## KCl
I_GCOMP_ZNS             = 7   ## ZnS
I_GCOMP_NA2S            = 8   ## Na2S
I_GCOMP_MNS             = 9   ## MnS
I_GCOMP_CR              = 10  ## Cr
I_GCOMP_FE              = 11  ## Fe
I_GCOMP_MG2SIO4         = 12  ## Mg2SiO4
I_GCOMP_TIO2            = 13 ## TiO2
I_GCOMP_AL2O3           = 14 ## Al2O3
I_GCOMP_CO              = 15 ## CO

# How is the CARMA group represented in the parent model
I_CNSTTYPE_PROGNOSTIC   = 1   ## Prognostic, advected constituent for each bin
I_CNSTTYPE_DIAGNOSTIC   = 2   ## Diagnostic, bins diagonosed from model state

# Return Codes
#
# NOTE: Also see error handling macros in globaer.h.
RC_OK             = 0   ## Success
RC_ERROR          = -1  ## Failure
RC_WARNING        = 1   ## Warning
RC_WARNING_RETRY  = 2   ## Warning, Retry Suggested


#  Define values of symbols used to specify horizontal & vertical grid type.
#   Grid selection is made by defining each of the variables
#   <igridv> and <igridh> to one of the grid types known to the model.
#
#   Possible values for igridv:
#       I_CART    cartesian
#       I_SIG     sigma
#       I_HYBRID  hybrid
#       I_LOGP    Log-Pressure #DPOW EDIT#
#
#    Possible values for igridh:
#       I_CART   cartesian
#       I_LL     longitude_latitude
#       I_LC     lambert_conformal
#       I_PS     polar_stereographic
#       I_ME     mercator
I_CART       = 1   ## Cartesian
I_SIG        = 2   ## Sigma
I_LL         = 3   ## Longitude & Latitude
I_LC         = 4   ## Lambert Conformal
I_PS         = 5   ## Polar Sterographic
I_ME         = 6   ## Mercator
I_HYBRID     = 7   ## Hybrid
I_LOGP   	   = 8   ## Log-Pressure


WTMOL_H2O = 18.016         # Water vapor
WTMOL_H2SO4 = 98.079       # Sulfuric Acid           
WTMOL_S8 = 256.48          # Sulfur (S8)  
WTMOL_S2 = 256.48 / 4.   # Sulfur (S2)
WTMOL_KCL = 74.5           # KCl
WTMOL_ZNS = 97.474         # ZnS
WTMOL_ZN = 65.38           # Zn
WTMOL_NA2S = 78.0452       # Na2S
WTMOL_NA2 = 2. * 22.9898 # Na2
WTMOL_NA = 22.9898 # Na2
WTMOL_MNS = 87.003         # MnS
WTMOL_MN = 54.938          # Mn
WTMOL_CR = 51.9961         # Cr
WTMOL_FE = 55.845          # Fe
WTMOL_MG2SIO4 = 140.69     # Mg2SiO4
WTMOL_MG2 = 2. * 24.305  # Mg2
WTMOL_MG = 24.305          # Mg
WTMOL_TIO2 = 79.866        # TiO2
WTMOL_TIO = 63.866         # TiO
WTMOL_AL2O3 = 101.961      # Al2O3
WTMOL_AL2 = 2. * 26.98   # Al2
WTMOL_AL = 26.98   # Al
WTMOL_CO = 28.01           # CO

	
	## Define mass density of condensates [ g / cm^3 ]
RHO_W = 1.			# liquid water
RHO_I = 0.93		# water ice
RHO_H2SO4 = 1.84		# sulfuric acid
RHO_SX = 1.96		# sulfur
RHO_KCL = 1.988		# KCl
RHO_ZNS = 4.04		# ZnS
RHO_NA2S = 1.856		# Na2S
RHO_MNS = 4.0		# MnS
RHO_CR = 7.15		# Cr
RHO_FE = 7.87		# Fe
RHO_MG2SIO4 = 3.21		# Mg2SiO4
RHO_TIO2 = 4.25		# TiO2
RHO_AL2O3 = 3.99		# Al2O3
RHO_CO = 1.0288		# CO, Bierhals J; Ullmann's Encyclopedia of Industrial Chemistry.


TCOEFF_KCL = 11382.
TCOEFF_ZNS = 15873.
TCOEFF_NA2S = 13889.
TCOEFF_MNS = 23810.
TCOEFF_CR = 20592.
TCOEFF_MG2SIO4 = 32488.   # Visscher notes
TCOEFF_FE = 20995.    # Visscher et al. 2010, ApJ 716, 1060
TCOEFF_TIO2_LODDERS = 34602.    # Diana email
TCOEFF_TIO2_HELLING = 32456.8678    # Helling et al. 2001
TCOEFF_AL2O3 = 45892.6    # Wakeford et al. 2017

	 # # Value of Collision Diameters [cm]
COLDIA_H2O = 3.11e-8   # Jacobson 2005 (book)
COLDIA_H2SO4 = 4.3e-8   # Michael J. Mills 1996, thesis  (Stratospheric Sulfate Aerosol: A Microphysical Model)
COLDIA_S8 = 6e-8   # Estimated from bond length, from Meyer 1976
COLDIA_S2 = 2e-8   # Estimated from bond length, from Meyer 1976
COLDIA_KCL = 3.08e-8   # Estimated from bond length, from Chemical Bonds and Bonds Energy by R Sanderson, 1976
COLDIA_ZNS =  3.665e-8  # Estimated from bond length, from Zack & Ziurys (2009), Journal of Molecular Spectroscopy 257, 213
COLDIA_NA2S =  4.195e-8  # Estimated from unit cell volume, Glasser & Jenkins 2000
COLDIA_MNS =  3.675e-8  # Estimated from unit cell volume, http://pveducation.org/pvcdrom/materials/MnS
COLDIA_CR =  3.655e-8  # Atomic radius
COLDIA_FE =  3.695e-8  # Atomic radius
COLDIA_MG2SIO4 =  3.845e-8  # Estimated from unit cell volume, Glasser & Jenkins 2000
COLDIA_TIO2 =  3.385e-8  # Estimated from Helling et al. (2001)
COLDIA_AL2O3 =  3.825e-8  # Estimated from Dobrovinskaya et al. 2009
COLDIA_CO =  3.86e-8  # Ramos-Estrada M. et al. Latin Am Appl Res. 34, 41-47, 2004



# p = 10^(offset - t_coeff/T - met_coeff * met - logp_coeff * log10 p)



gas_dict = {
    "H2O": {
        "cond_rho": RHO_W,
        "vaprtn": {
            "rtn": I_VAPRTN_H2O_MURPHY2005,
            "offset": 0,
            "T_coeff": 0,
            "met_coeff": 0,
            "logp_coeff": 0
        },
        "wtmol": WTMOL_H2O,
        "wtmol_dif": WTMOL_H2O,
        "gcomp": I_GCOMP_H2O,
        "mucos_dict": {},
        "fastchem_specie": "H2O1",
        "opacity_files": "WS15/H2O_complex.txt"
    },
        "KCl": {
        "cond_rho": RHO_KCL,
        "vaprtn": {
            "rtn": I_VAPRTN_KCL_MORLEY2012,
            "offset": 13.6106,
            "T_coeff": TCOEFF_KCL,
            "met_coeff": 0,
            "logp_coeff": 0
        },
        "wtmol": WTMOL_KCL,
        "wtmol_dif": WTMOL_KCL,
        "gcomp": I_GCOMP_KCL,
        "mucos_dict": {},
        "fastchem_specie": "Cl1K1",
        "opacity_files": "WS15/KCl_complex.txt"
    },
        "ZnS": {
        "cond_rho": RHO_ZNS,
        "vaprtn": {
            "rtn": I_VAPRTN_ZNS_MORLEY2012,
            "offset": 18.812,
            "T_coeff": TCOEFF_ZNS,
            "met_coeff": 1,
            "logp_coeff": 0
        },
        "wtmol": WTMOL_ZNS,
        "wtmol_dif": WTMOL_ZN,
        "gcomp": I_GCOMP_ZNS,
        "mucos_dict": {"KCl": 0.144356},
        "fastchem_specie": "Zn",
        "opacity_files": "WS15/ZnS_complex.txt"
    },
        "Na2S": {
        "cond_rho": RHO_NA2S,
        "vaprtn": {
            "rtn": I_VAPRTN_NA2S_MORLEY2012,
            "offset": 14.55,
            "T_coeff": TCOEFF_ZNS,
            "met_coeff": 0.5,
            "logp_coeff": 0
        },
        "wtmol": WTMOL_NA2S,
        "wtmol_dif": WTMOL_NA,
        "gcomp": I_GCOMP_NA2S,
        "mucos_dict": {"TiO2": 0.48481},
        "fastchem_specie": "Na",
        "opacity_files": "WS15/Na2S_complex.txt"
    },
        "MnS": {
        "cond_rho": RHO_MNS,
        "vaprtn": {
            "rtn": I_VAPRTN_MNS_MORLEY2012,
            "offset": 17.532,
            "T_coeff": TCOEFF_MNS,
            "met_coeff": 1,
            "logp_coeff": 0
        },
        "wtmol": WTMOL_MNS,
        "wtmol_dif": WTMOL_MN,
        "gcomp": I_GCOMP_MNS,
        "mucos_dict": {"TiO2": 0.214735},
        "fastchem_specie": "Na",
        "opacity_files": "KH18/MnS-KH_complex.txt"
    },
        "Cr": {
        "cond_rho": RHO_CR,
        "vaprtn": {
            "rtn": I_VAPRTN_CR_MORLEY2012,
            "offset": 13.49,
            "T_coeff": TCOEFF_CR,
            "met_coeff": 0,
            "logp_coeff": 0
        },
        "wtmol": WTMOL_CR,
        "wtmol_dif": WTMOL_CR,
        "gcomp": I_GCOMP_CR,
        "mucos_dict": {"TiO2": 0.262189},
        "fastchem_specie": "Cr",
        "opacity_files": "KH18/Cr_complex.txt"
    },
        "Mg2SiO4": {
        "cond_rho": RHO_MG2SIO4,
        "vaprtn": {
            "rtn": I_VAPRTN_MG2SIO4_VISSCHER2010,
            "offset": 22.08,
            "T_coeff": TCOEFF_MG2SIO4,
            "met_coeff": 1.4,
            "logp_coeff": 0.2
        },
        "wtmol": WTMOL_MG2SIO4,
        "wtmol_dif": WTMOL_MG,
        "gcomp": I_GCOMP_MG2SIO4,
        "mucos_dict": {"TiO2": 0.995},
        "fastchem_specie": "Mg",
        "opacity_files": "B21/Mg2SiO4_amorph.txt"
    },
        "Fe": {
        "cond_rho": RHO_FE,
        "vaprtn": {
            "rtn": I_VAPRTN_FE_VISSCHER2010,
            "offset": 13.23,
            "T_coeff": TCOEFF_FE,
            "met_coeff": 0,
            "logp_coeff": 0
        },
        "wtmol": WTMOL_FE,
        "wtmol_dif": WTMOL_FE,
        "gcomp": I_GCOMP_FE,
        "mucos_dict": {"TiO2": 0.221548},
        "fastchem_specie": "Fe",
        "opacity_files": "KH18/Fe_complex.txt"
    },
        "TiO2": {
        "cond_rho": RHO_TIO2,
        "vaprtn": {
            "rtn": I_VAPRTN_TIO2_HELLING2001,
            "offset": 15.5489,
            "T_coeff": TCOEFF_TIO2_HELLING,
            "met_coeff": 0,
            "logp_coeff": 0
        },
        "wtmol": WTMOL_TIO2,
        "wtmol_dif": WTMOL_TIO2,
        "gcomp": I_GCOMP_TIO2,
        "mucos_dict": {},
        "fastchem_specie": "O2Ti1",
        "opacity_files": "gCMCRT/TiO2[s].txt"
    },
        "Al2O3": {
        "cond_rho": RHO_AL2O3,
        "vaprtn": {
            "rtn": I_VAPRTN_AL2O3_WAKEFORD2017,
            "offset": 23.7,
            "T_coeff": TCOEFF_AL2O3,
            "met_coeff": 1.66,
            "logp_coeff": 0
        },
        "wtmol": WTMOL_AL2O3,
        "wtmol_dif": WTMOL_AL,
        "gcomp": I_GCOMP_AL2O3,
        "mucos_dict": {"TiO2": 0.724172},
        "fastchem_specie": "Al",
        "opacity_files": "KH18/Al2O3[s].txt"
    }
}



cond_rho = {
    "KCl": RHO_KCL,
    "ZnS": RHO_ZNS,
    "Na2S": RHO_NA2S,
    "MnS": RHO_MNS,
    "Cr": RHO_CR,
    "Mg2SiO4": RHO_MG2SIO4,
    "Fe": RHO_FE,
    "TiO2": RHO_TIO2,
    "Al2O3": RHO_AL2O3
}




# igas_dict = {
#     "KCl": 1,
#     "ZnS": 2,
#     "Na2S": 3,
#     "MnS": 4,
#     "Cr": 5,
#     "Mg2SiO4": 6,
#     "Fe": 7,
#     "TiO2": 8,
#     "Al2O3": 9
# }

vaprtn_dict = {
    "H2O": I_VAPRTN_H2O_MURPHY2005,
    "KCl": I_VAPRTN_KCL_MORLEY2012,
    "ZnS": I_VAPRTN_ZNS_MORLEY2012,
    "Na2S": I_VAPRTN_NA2S_MORLEY2012,
    "MnS": I_VAPRTN_MNS_MORLEY2012,
    "Cr": I_VAPRTN_CR_MORLEY2012,
    "Mg2SiO4": I_VAPRTN_MG2SIO4_VISSCHER2010,
    "Fe": I_VAPRTN_FE_VISSCHER2010,
    "TiO2": I_VAPRTN_TIO2_HELLING2001,
    "Al2O3": I_VAPRTN_AL2O3_WAKEFORD2017
}

wtmol_dict = {
    "H2O": WTMOL_H2O,
    "KCl": WTMOL_KCL,
    "ZnS": WTMOL_ZNS,
    "Na2S": WTMOL_NA2S,
    "MnS": WTMOL_MNS,
    "Cr": WTMOL_CR,
    "Mg2SiO4": WTMOL_MG2SIO4,
    "Fe": WTMOL_FE,
    "TiO2": WTMOL_TIO2,
    "Al2O3": WTMOL_AL
}

wtmol_dif_dict = {
    "H2O": WTMOL_H2O,
    "KCl": WTMOL_KCL,
    "ZnS": WTMOL_ZN,
    "Na2S": WTMOL_NA,
    "MnS": WTMOL_MN,
    "Cr": WTMOL_CR,
    "Mg2SiO4": WTMOL_MG,
    "Fe": WTMOL_FE,
    "TiO2": WTMOL_TIO2,
    "Al2O3": WTMOL_AL
}

gcomp_dict = { 
    "H2O": I_GCOMP_H2O,           
    "H2SO4": I_GCOMP_H2SO4,          
    "SO2": I_GCOMP_SO2,         
    "S8": I_GCOMP_S8,         
    "S2": I_GCOMP_S2,      
    "KCl": I_GCOMP_KCL,    
    "ZnS": I_GCOMP_ZNS,          
    "Na2S": I_GCOMP_NA2S,        
    "MnS": I_GCOMP_MNS,           
    "Cr": I_GCOMP_CR,          
    "Fe": I_GCOMP_FE,             
    "Mg2SiO4": I_GCOMP_MG2SIO4,      
    "TiO2": I_GCOMP_TIO2,          
    "Al2O3": I_GCOMP_AL2O3,           
    "CO": I_GCOMP_CO             
              }

mucos_dict = {
    "H2O": {},
    "KCl": {},
    "ZnS": {"KCl": 0.144356},
    "Na2S": {"TiO2": 0.48481},
    "MnS": {"TiO2": 0.214735},
    "Cr": {"TiO2": 0.262189},
    "Mg2SiO4": {"TiO2": 0.995},
    "Fe": {"TiO2": 0.221548},
    "TiO2": {},
    "Al2O3": {"TiO2": 0.724172}
}


fastchem_species = {
    "H2O": "H2O",
    "TiO2": "O2Ti1",
    "Fe": "Fe",
    "Mg2SiO4": "Mg",
    "Cr": "Cr",
    "MnS": "Mn",
    "Na2S": "Na",
    "ZnS": "Zn",
    "KCl": "Cl1K1",
    "Al2O3": "Al"
}

opacity_files = {
    "H2O": "WS15/H2O_ice.txt",
    "KCl": "WS15/KCl_complex.txt",
    "ZnS": "WS15/ZnS_complex.txt",
    "Na2S": "WS15/Na2S_complex.txt",
    "MnS": "KH18/MnS-KH_complex.txt",
    "Cr": "KH18/Cr_complex.txt",
    "Mg2SiO4": "B21/Mg2SiO4_amorph.txt",
    "Fe": "KH18/Fe_complex.txt",
    "TiO2": "gCMCRT/TiO2[s].txt", 
    "Al2O3": "KH18/Al2O3[s].txt"
}



JUPITER_RADIUS = 6.995e9 #cm
BAR_TO_BARYE = 1e6
k_B = 1.381e-16 #erg/K
PROTON_MASS =  1.673e-24 # g

PARMENTIER_A_COEFF = 0.32
PARMENTIER_B_COEFF = 1/30000

