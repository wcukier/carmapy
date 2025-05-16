! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  Calculates the vapor pressure for sulfuric acid.
!!
!!  <pvap_liq> and <pvap_ice> are vapor pressures in units of [dyne/cm^2]
!!
!!  Created   Dec-1995  (Ackerman) 
!!  Modified  Sep-1997  (McKie)
!!  Modified Jul-2001 (Mills)
!!
!!  NOTE: To calculate vapor pressure of H2SO4 water vapor pressure (pvapl(iz, igash2o))
!!  should be calculated before this calculation.
!!
!! @author Mike Mills, Tianyi Fan
!! @version Feb-2011
subroutine vaporp_co_wylie1958(carma, cstate, iz, rc, pvap_liq, pvap_ice)
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
  real(kind=f), intent(out)            :: pvap_liq  !! vapor pressure wrt liquid [dyne/cm2]
  real(kind=f), intent(out)            :: pvap_ice  !! vapor pressure wrt ice [dyne[cm2]
  integer, intent(inout)               :: rc        !! return code, negative indicates failure

  ! Local declarations

  real(kind=f)                   :: tt
  real(kind=f)                   :: slope
  real(kind=f)                   :: b

  ! Saturation vapor presure of CO from thesis by Wylie (1958), extrapolated at
  ! higher temperatures. 

  if (t(iz) .lt. 61.57_f) then
    pvap_liq = 1333.2239_f * 10._f ** (2.4482_f - 418.44_f / t(iz) + 4.134_f * log10(t(iz)) - 0.02599_f * t(iz))
  else
    tt = 61.57_f
    slope = tt*log(10._f) * (418.44_f/tt**2._f + 4.134_f/(log(10._f)*tt) - 0.02599_f)
    b = 2.4482_f - 418.44_f/tt + 4.134_f*log10(tt) - 0.02599_f*tt - slope*log10(tt)
    pvap_liq = 1333.2239_f * 10._f ** (slope*log10(t(iz)) + b)
  endif

  pvap_ice = pvap_liq

  return
end
