module carma_condensate_mod

  use carma_precision_mod
  use carma_constants_mod

  implicit none

	!--
	! Condensate Constants
	
	! Meter-Kilogram-Second (MKS) convention for units
	! This convention is different from CARMA's original 
	!  Centimeter-Gram-Second (CGS) convention.  Be wary of
	!  this conversion to the new convention.
	
	! Use the _f for all literal constants, e.g. 1.2e_f.
	! If you omit the _f in the initialization, a compiler may cast this
	!  number into single precision and then store it as _f precision.
	
	!! Define molecular weight of condensates [ g / mole ]
	real(kind=f), parameter :: WTMOL_H2O = 18.016_f         ! Water vapor
	real(kind=f), parameter :: WTMOL_H2SO4 = 98.079_f       ! Sulfuric Acid           
	real(kind=f), parameter :: WTMOL_S8 = 256.48_f          ! Sulfur (S8)  
	real(kind=f), parameter :: WTMOL_S2 = 256.48_f / 4._f   ! Sulfur (S2)
	real(kind=f), parameter :: WTMOL_KCL = 74.5_f           ! KCl
	real(kind=f), parameter :: WTMOL_ZNS = 97.474_f         ! ZnS
	real(kind=f), parameter :: WTMOL_ZN = 65.38_f           ! Zn
	real(kind=f), parameter :: WTMOL_NA2S = 78.0452_f       ! Na2S
	real(kind=f), parameter :: WTMOL_NA2 = 2._f * 22.9898_f ! Na2
	real(kind=f), parameter :: WTMOL_NA = 22.9898_f ! Na2
	real(kind=f), parameter :: WTMOL_MNS = 87.003_f         ! MnS
	real(kind=f), parameter :: WTMOL_MN = 54.938_f          ! Mn
	real(kind=f), parameter :: WTMOL_CR = 51.9961_f         ! Cr
	real(kind=f), parameter :: WTMOL_FE = 55.845_f          ! Fe
	real(kind=f), parameter :: WTMOL_MG2SIO4 = 140.69_f     ! Mg2SiO4
	real(kind=f), parameter :: WTMOL_MG2 = 2._f * 24.305_f  ! Mg2
	real(kind=f), parameter :: WTMOL_MG = 24.305_f          ! Mg
	real(kind=f), parameter :: WTMOL_TIO2 = 79.866_f        ! TiO2
	real(kind=f), parameter :: WTMOL_TIO = 63.866_f         ! TiO
	real(kind=f), parameter :: WTMOL_AL2O3 = 101.961_f      ! Al2O3
	real(kind=f), parameter :: WTMOL_AL2 = 2._f * 26.98_f   ! Al2
	real(kind=f), parameter :: WTMOL_AL = 26.98_f   ! Al
	real(kind=f), parameter :: WTMOL_CO = 28.01_f           ! CO

	
	!! Define mass density of condensates [ g / cm^3 ]
	real(kind=f), parameter :: RHO_W = 1._f			! liquid water
	real(kind=f), parameter :: RHO_I = 0.93_f		! water ice
	real(kind=f), parameter :: RHO_H2SO4 = 1.84_f		! sulfuric acid
	real(kind=f), parameter :: RHO_SX = 1.96_f		! sulfur
	real(kind=f), parameter :: RHO_KCL = 1.988_f		! KCl
	real(kind=f), parameter :: RHO_ZNS = 4.04_f		! ZnS
	real(kind=f), parameter :: RHO_NA2S = 1.856_f		! Na2S
	real(kind=f), parameter :: RHO_MNS = 4.0_f		! MnS
	real(kind=f), parameter :: RHO_CR = 7.15_f		! Cr
	real(kind=f), parameter :: RHO_FE = 7.87_f		! Fe
	real(kind=f), parameter :: RHO_MG2SIO4 = 3.21_f		! Mg2SiO4
	real(kind=f), parameter :: RHO_TIO2 = 4.25_f		! TiO2
	real(kind=f), parameter :: RHO_AL2O3 = 3.99_f		! Al2O3
	real(kind=f), parameter :: RHO_CO = 1.0288_f		! CO, Bierhals J; Ullmann's Encyclopedia of Industrial Chemistry.


	!! Latent heat of evaporation for water [cm^2/s^2]
	real(kind=f), parameter :: RLHE_CNST = 2.501e10_f	
	
	!! Latent heat of ice melting for water [cm^2/s^2]
	real(kind=f), parameter :: RLHM_CNST = 3.337e9_f

	!! Latent heat of sulfur [cm^2/s^2]
	real(kind=f), parameter :: RLH_CNST_SX = 2.9e9_f

	!! Value of B in log(P) = A - B/T for saturation vapor pressure of exoplanet clouds (Morley et al. 2012) [K]
	real(kind=f), parameter :: TCOEFF_KCL = 11382._f
	real(kind=f), parameter :: TCOEFF_ZNS = 15873._f
	real(kind=f), parameter :: TCOEFF_NA2S = 13889._f
	real(kind=f), parameter :: TCOEFF_MNS = 23810._f
	real(kind=f), parameter :: TCOEFF_CR = 20592._f
	real(kind=f), parameter :: TCOEFF_MG2SIO4 = 32488._f   ! Visscher notes
	real(kind=f), parameter :: TCOEFF_FE = 20995._f   ! Visscher et al. 2010, ApJ 716, 1060
	real(kind=f), parameter :: TCOEFF_TIO2_LODDERS = 34602._f   ! Diana email
	real(kind=f), parameter :: TCOEFF_TIO2_HELLING = 32456.8678_f   ! Helling et al. 2001
	real(kind=f), parameter :: TCOEFF_AL2O3 = 45892.6_f   ! Wakeford et al. 2017

	!! Value of Collision Diameters [cm]
  	real(kind=f), parameter :: COLDIA_H2O = 3.11e-8_f  ! Jacobson 2005 (book)
  	real(kind=f), parameter :: COLDIA_H2SO4 = 4.3e-8_f  ! Michael J. Mills 1996, thesis  (Stratospheric Sulfate Aerosol: A Microphysical Model)
  	real(kind=f), parameter :: COLDIA_S8 = 6e-8_f  ! Estimated from bond length, from Meyer 1976
  	real(kind=f), parameter :: COLDIA_S2 = 2e-8_f  ! Estimated from bond length, from Meyer 1976
  	real(kind=f), parameter :: COLDIA_KCL = 3.08e-8_f  ! Estimated from bond length, from Chemical Bonds and Bonds Energy by R Sanderson, 1976
  	real(kind=f), parameter :: COLDIA_ZNS =  3.665e-8_f ! Estimated from bond length, from Zack & Ziurys (2009), Journal of Molecular Spectroscopy 257, 213
  	real(kind=f), parameter :: COLDIA_NA2S =  4.195e-8_f ! Estimated from unit cell volume, Glasser & Jenkins 2000
  	real(kind=f), parameter :: COLDIA_MNS =  3.675e-8_f ! Estimated from unit cell volume, http://pveducation.org/pvcdrom/materials/MnS
  	real(kind=f), parameter :: COLDIA_CR =  3.655e-8_f ! Atomic radius
  	real(kind=f), parameter :: COLDIA_FE =  3.695e-8_f ! Atomic radius
  	real(kind=f), parameter :: COLDIA_MG2SIO4 =  3.845e-8_f ! Estimated from unit cell volume, Glasser & Jenkins 2000
  	real(kind=f), parameter :: COLDIA_TIO2 =  3.385e-8_f ! Estimated from Helling et al. (2001)
  	real(kind=f), parameter :: COLDIA_AL2O3 =  3.825e-8_f ! Estimated from Dobrovinskaya et al. 2009
  	real(kind=f), parameter :: COLDIA_CO =  3.86e-8_f ! Ramos-Estrada M. et al. Latin Am Appl Res. 34, 41-47, 2004

end module 
