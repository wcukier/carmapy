! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine defines time-independent parameters used to calculate
!! condensational growth/evaporation.
!!
!! The parameters defined for each gas are 
!1>
!!   gwtmol:   molecular weight [g/mol]
!!   diffus:   diffusivity      [cm^2/s]
!!   rlhe  :   latent heat of evaporation [cm^2/s^2]
!!   rlhm  :   latent heat of melting [cm^2/s^2]
!!<
!! Time-independent parameters that depend on particle radius are
!! defined in setupgkern.f.
!!
!! This routine requires that vertical profiles of temperature <T>,
!! and pressure <p> are defined.
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine setupgrow(carma, cstate, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_planet_mod
  use carma_condensate_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc       !! return code, negative indicates failure
  
  ! Local Variable
  integer                        :: ielem    !! element index
  integer                        :: k        !! z index
  integer                        :: i
  real(kind=f)                   :: rhoa_cgs, aden, wtpctf
  ! Define formats
  1 format(a,':  ',12i6)
  2 format(a,':  ',i6)
  3 format(/' id  gwtmol   gasname',(/,i3,3x,f5.1,3x,a))
  5 format(/,'Particle growth mapping arrays (setupgrow):')


  !-----Check that values are valid------------------------------------------
  do ielem = 1, NELEM
    if( igrowgas(ielem) .gt. NGAS )then
      if (do_print) write(LUNOPRT,*) 'setupgrow::ERROR - component of igrowgas > NGAS'
      rc = -1
      return
    endif
  enddo

  ! Define parameters with weak time-dependence to be used in
  ! growth equation.
  do k = 1, NZ

    ! Diffusivity of water vapor in air from Pruppacher & Klett (eq. 13-3);
    ! units are [cm^2/s].
    if (igash2o .ne. 0) then 
      rhoa_cgs = rhoa(k) / (xmet(k) * ymet(k) * zmet(k))
!      diffus(k, igash2o) = 0.211_f * (1.01325e+6_f / p(k)) * (t(k) / 273.15_f )**1.94_f
      diffus(k, igash2o) = 5._f / (16._f * AVG * COLDIA_H2O**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_H2O + wtmol_air(k)) / (2._f * PI * WTMOL_H2O)) 
  
      ! Latent heat of evaporation for water; units are [cm^2/s^2]
      if (do_cnst_rlh) then
        rlhe(k, igash2o) = RLHE_CNST
      else
        ! from Stull
        rlhe(k, igash2o) = (2.5_f - .00239_f * (t(k) - 273.16_f)) * 1.e10_f      
      end if
  
      ! Latent heat of ice melting; units are [cm^2/s^2]
      if (do_cnst_rlh) then
        rlhm(k, igash2o) = RLHM_CNST
      else
  
        ! from Pruppacher & Klett (eq. 4-85b)
        !
        ! NOTE: This expression yields negative values for rlmh at mesospheric
        ! temperatures.
        rlhm(k, igash2o) = (79.7_f + 0.485_f * (t(k) - 273.16_f) - 2.5e-3_f * &
          ((t(k) - 273.16_f)**2)) * 4.186e7_f
      end if
    end if

    ! Properties for H2SO4
    if (igash2so4 .ne. 0) then
      ! Diffusivity
      rhoa_cgs = rhoa(k) / (xmet(k) * ymet(k) * zmet(k))
      aden     = rhoa_cgs * AVG / wtmol_air(k)
      diffus(k,igash2so4) = 5._f / (16._f * AVG * COLDIA_H2SO4**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_H2SO4 + wtmol_air(k)) / (2._f * PI * WTMOL_H2SO4)) 
  
      wtpctf = wtpct(k)/100._f
      ! HACK: make H2SO4 latent heats same as water
!      rlhe(k,igash2so4) = rlhe(k, igash2o)
!      rlhm(k,igash2so4) = rlhe(k, igash2o)
      ! From Jang et al. 2006, Transactions of the Korean Nuclear Society Autumn Meeting
      rlhe(k,igash2so4) = 1364.93*wtpctf**3._f - 1226.46*wtpctf**2._f + 382.23*wtpctf + 540.52
      rlhe(k,igash2so4) = rlhe(k,igash2so4) * 4.184e7_f
      rlhm(k,igash2so4) = rlhe(k,igash2so4)
    end if

    ! Properties for S8
    if (igass8 .ne. 0) then
      ! Diffusivity
      diffus(k,igass8) = 5._f / (16._f * AVG * COLDIA_S8**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_S8 + wtmol_air(k)) / (2._f * PI * WTMOL_S8)) 
      rlhe(k,igass8) = RLH_CNST_SX ! http://cameochemicals.noaa.gov/chris/SXX.pdf
      rlhm(k,igass8) = rlhe(k, igass8)
    end if

    ! Properties for S2
    if (igass2 .ne. 0) then
      ! Diffusivity
      diffus(k,igass2) = 5._f / (16._f * AVG * COLDIA_S2**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_S2 + wtmol_air(k)) / (2._f * PI * WTMOL_S2)) 
      rlhe(k,igass2) = RLH_CNST_SX ! http://cameochemicals.noaa.gov/chris/SXX.pdf
      rlhm(k,igass2) = rlhe(k, igass2)
    end if

    ! Properties for KCl
    if (igaskcl .ne. 0) then
      ! Diffusivity
      diffus(k,igaskcl) = 5._f / (16._f * AVG * COLDIA_KCL**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_KCL + wtmol_air(k)) / (2._f * PI * WTMOL_KCL)) 
      rlhe(k,igaskcl) = TCOEFF_KCL * log(10._f) * RGAS / WTMOL_KCL ! Charnay et al. 2015, ApJL 813, L1
      rlhm(k,igaskcl) = rlhe(k, igaskcl)
    end if

    ! Properties for ZnS
    if (igaszns .ne. 0) then
      ! Diffusivity
      diffus(k,igaszns) = 5._f / (16._f * AVG * COLDIA_ZNS**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_ZN + wtmol_air(k)) / (2._f * PI * WTMOL_ZN)) 
      rlhe(k,igaszns) = TCOEFF_ZNS * log(10._f) * RGAS / WTMOL_ZN ! Charnay et al. 2015, ApJL 813, L1
      rlhm(k,igaszns) = rlhe(k, igaszns)
    end if

    ! Properties for Na2S
    if (igasna2s .ne. 0) then
      ! Diffusivity
      diffus(k,igasna2s) = 0.5_f * 5._f / (16._f * AVG * COLDIA_NA2S**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_NA + wtmol_air(k)) / (2._f * PI * WTMOL_NA)) 
      rlhe(k,igasna2s) = TCOEFF_NA2S * log(10._f) * RGAS / WTMOL_NA ! Charnay et al. 2015, ApJL 813, L1
      rlhm(k,igasna2s) = rlhe(k, igasna2s)
    end if

    ! Properties for MnS
    if (igasmns .ne. 0) then
      ! Diffusivity
      diffus(k,igasmns) = 5._f / (16._f * AVG * COLDIA_MNS**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_MN + wtmol_air(k)) / (2._f * PI * WTMOL_MN)) 
      rlhe(k,igasmns) = TCOEFF_MNS * log(10._f) * RGAS / WTMOL_MN ! Charnay et al. 2015, ApJL 813, L1
      rlhm(k,igasmns) = rlhe(k, igasmns)
    end if

    ! Properties for Cr
    if (igascr .ne. 0) then
      ! Diffusivity
      diffus(k,igascr) = 5._f / (16._f * AVG * COLDIA_CR**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_CR + wtmol_air(k)) / (2._f * PI * WTMOL_CR)) 
      rlhe(k,igascr) = TCOEFF_CR * log(10._f) * RGAS / WTMOL_CR ! Charnay et al. 2015, ApJL 813, L1
      rlhm(k,igascr) = rlhe(k, igascr)
    end if

    ! Properties for Fe
    if (igasfe .ne. 0) then
      ! Diffusivity
      diffus(k,igasfe) = 5._f / (16._f * AVG * COLDIA_FE**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_FE + wtmol_air(k)) / (2._f * PI * WTMOL_FE)) 
      rlhe(k,igasfe) = TCOEFF_FE * log(10._f) * RGAS / WTMOL_FE ! Charnay et al. 2015, ApJL 813, L1
      rlhm(k,igasfe) = rlhe(k, igasfe)
    end if

    ! Properties for Mg2SiO4
    if (igasmg2sio4 .ne. 0) then
      ! Diffusivity
      diffus(k,igasmg2sio4) = 0.5_f *  5._f / (16._f * AVG * COLDIA_MG2SIO4**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_MG + wtmol_air(k)) / (2._f * PI * WTMOL_MG)) 
      rlhe(k,igasmg2sio4) = TCOEFF_MG2SIO4 * log(10._f) * RGAS / WTMOL_MG ! Charnay et al. 2015, ApJL 813, L1; neglecting PT part
      rlhm(k,igasmg2sio4) = rlhe(k, igasmg2sio4)
    end if

    ! Properties for TiO2
    if (igastio2 .ne. 0) then
      ! Diffusivity
      diffus(k,igastio2) = 5._f / (16._f * AVG * COLDIA_TIO2**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_TIO2 + wtmol_air(k)) / (2._f * PI * WTMOL_TIO2)) 
      rlhe(k,igastio2) = TCOEFF_TIO2_HELLING * log(10._f) * RGAS / WTMOL_TIO2 ! Charnay et al. 2015, ApJL 813, L1
      rlhm(k,igastio2) = rlhe(k, igastio2)
    end if

    ! Properties for Al2O3
    if (igasal2o3 .ne. 0) then
      ! Diffusivity
      diffus(k,igasal2o3) = 0.5_f * 5._f / (16._f * AVG * COLDIA_AL2O3**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_AL + wtmol_air(k)) / (2._f * PI * WTMOL_AL)) 
      rlhe(k,igasal2o3) = TCOEFF_AL2O3 * log(10._f) * RGAS / WTMOL_AL ! Charnay et al. 2015, ApJL 813, L1
      rlhm(k,igasal2o3) = rlhe(k, igasal2o3)
    end if

    ! Properties for CO
    if (igasco .ne. 0) then
      ! Diffusivity
      diffus(k,igasco) = 5._f / (16._f * AVG * COLDIA_CO**2_f * rhoa_cgs * COLINT) * &
  	sqrt(RGAS * t(k) * wtmol_air(k) * (WTMOL_CO + wtmol_air(k)) / (2._f * PI * WTMOL_CO)) 
      rlhe(k,igasco) = 2.63e9_f ! @ triple point; Bierhals et al. Ullmann's Encyclopedia of Industrial Chemistry
      rlhm(k,igasco) = rlhe(k, igasco)
    end if

    
  enddo

#ifdef DEBUG
  ! Report some initialization values
  if (do_print_init) then
    write(LUNOPRT,5)
    write(LUNOPRT,2) 'NGAS    ',NGAS
    write(LUNOPRT,1) 'igrowgas',(igrowgas(i),i=1,NELEM)
    write(LUNOPRT,3) (i,gwtmol(i),gasname(i),i=1,NGAS)
  endif
#endif

  ! Return to caller with particle growth mapping arrays and time-dependent
  ! parameters initialized.
  return
end
