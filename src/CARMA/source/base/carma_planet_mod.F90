module carma_planet_mod

  use carma_precision_mod
  use carma_constants_mod

  implicit none

	!--
	! Planet Constants
	
	! Meter-Kilogram-Second (MKS) convention for units
	! This convention is different from CARMA's original 
	!  Centimeter-Gram-Second (CGS) convention.  Be wary of
	!  this conversion to the new convention.
	
	! Use the _f for all literal constants, e.g. 1.2e_f.
	! If you omit the _f in the initialization, a compiler may cast this
	!  number into single precision and then store it as _f precision.
	
	!! Acceleration of gravity near the planet's surface [ cm/s^2 ]
!	real(kind=f), parameter :: GRAV = 980.6_f			!EARTH
!	real(kind=f), parameter :: GRAV = 887.0_f			!VENUS
!	real(kind=f), parameter :: GRAV = 893.0_f			!GJ1214B
!	real(kind=f), parameter :: GRAV = 3929.0_f                      !GAMMA CEPHEI AB
!	real(kind=f), parameter :: GRAV = 337.0_f                       !HD 192310 C
!	real(kind=f), parameter :: GRAV = 22895._f                      !UPS AND D
!	real(kind=f), parameter :: GRAV = 1780._f                       ! Low-g grid model
!	real(kind=f), parameter :: GRAV = 17800._f                       ! Mid-g grid model
!	real(kind=f), parameter :: GRAV = 178000._f                       ! High-g grid model
	
	!! Define planet equatorial radius [ cm ]
!	real(kind=f), parameter :: RPLANET  = 6.37e+8_f			!EARTH
!	real(kind=f), parameter :: RPLANET  = 6.05e+8_f			!VENUS
!	real(kind=f), parameter :: RPLANET  = 1.708e9_f                 !GJ1214B
!	real(kind=f), parameter :: RPLANET  = 6.991e9_f                 !GAMMA CEPHEI AB, UPS AND D (JUPITER)
!	real(kind=f), parameter :: RPLANET  = 1.036_f * 6.991e9_f       ! WASP-43b
!	real(kind=f), parameter :: RPLANET  = 0.93666_f * 6.991e9_f       ! HAT-P-12b
!	real(kind=f), parameter :: RPLANET  = 2.462e9_f                 !HD 192310 C (NEPTUNE)
!	real(kind=f), parameter :: RPLANET  = 2.58e+8_f                 !TITAN
!	real(kind=f), parameter :: RPLANET  = 1.19e+8_f                 !PLUTO

	!! Define molecular weight of dry air [ g / mole ]
!	real(kind=f), parameter :: WTMOL_AIR = 28.966_f			!EARTH
!	real(kind=f), parameter :: WTMOL_AIR = 43.450_f			!VENUS
!	real(kind=f), parameter :: WTMOL_AIR = 2.3202_f                 !GJ1214B
	!real(kind=f), parameter :: WTMOL_AIR = 2.2_f                 ! GAS GIANTS (H/He-Dominated)
	
	!! Define gas constant for dry air [ erg / deg_K / mole ]
	!real(kind=f), parameter :: R_AIR = RGAS / WTMOL_AIR
	
	!! Define number of seconds per the planet's day [ s / d ]
!	real(kind=f), parameter :: SCDAY = 86400._f			!EARTH
!	real(kind=f), parameter :: SCDAY = 10108800._f			!VENUS
!	real(kind=f), parameter :: SCDAY = 1.58_f * 86400._f            !GJ1214B
!	real(kind=f), parameter :: SCDAY = 36000._f                     !10 hours, generic gas giant
!	real(kind=f), parameter :: SCDAY = 1377648._f                   !TITAN
!	real(kind=f), parameter :: SCDAY = 551815.2_f                   !PLUTO

        !! Define specific heat at constant pres of dry air [ cm^2 / s^2 / deg_K ]
!	real(kind=f), parameter :: CP = 1.004e7_f			!EARTH
!	real(kind=f), parameter :: CP = 7.0e6_f				!VENUS
	real(kind=f), parameter :: CP = 1.3e8_f                         !Hydrogen, for generic gas giant, Kataria et al. 2015, ApJ 801, 86
!	real(kind=f), parameter :: CP = 1.040e7_f                       !TITAN, PLUTO

	!! Define ratio of gas constant for dry air and specific heat
	!real(kind=f), parameter :: RKAPPA = R_AIR / CP
	
end module 
