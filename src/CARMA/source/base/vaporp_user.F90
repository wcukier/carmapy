! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  Calculates the vapor pressure for a user defined condensate.
!!
!!  <pvap_liq> and <pvap_ice> are vapor pressures in units of [dyne/cm^2]
!!
!!  Created   Dec-1995  (Ackerman) 
!!  Modified  Sep-1997  (McKie)
!!  Modified Jul-2001 (Mills)
!!  Modified Jul-2025 (Cukier)
!!
!!
!! @author Mike Mills, Tianyi Fan, Wolf Cukier
!! @version Jul-2025
subroutine vaporp_user(carma, cstate, iz, igas, rc, pvap_liq, pvap_ice)
!     types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_planet_mod
  use carma_condensate_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use sulfate_utils
  
  implicit none

  type(carma_type), intent(inout)      :: carma     !! the carma object
  type(carmastate_type), intent(inout) :: cstate    !! the carma state object
  integer, intent(in)                  :: iz        !! z index
  integer, intent(in)                  :: igas      !! gas index
  real(kind=f), intent(out)            :: pvap_liq  !! vapor pressure wrt liquid [dyne/cm2]
  real(kind=f), intent(out)            :: pvap_ice  !! vapor pressure wrt ice [dyne/cm2]
  integer, intent(inout)               :: rc        !! return code, negative indicates failure
  real(kind=f)                         :: offset, tcoeff, metcoeff, logpcoeff

  offset    = carma%f_gas(igas)%f_vp_offset
  tcoeff    = carma%f_gas(igas)%f_vp_tcoeff
  metcoeff  = carma%f_gas(igas)%f_vp_metcoeff
  logpcoeff = carma%f_gas(igas)%f_vp_logpcoeff

  pvap_liq = 1e6_f * 10._f ** (offset - tcoeff/t(iz) - metcoeff*met - logpcoeff*log10(1e-6_f*p(iz)))
  pvap_ice = pvap_liq

  return
end
