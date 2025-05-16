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
subroutine vaporp_s2_lyons2008(carma, cstate, iz, rc, pvap_liq, pvap_ice)
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

  ! Saturation vapor presure of (1) liquid, (2) solid S2 from Lyons 2008. 
  if (t(iz) .gt. 392._f) then 
    pvap_liq = 1.01325e6 * 10._f ** (7.0240_f - 6091.2_f/t(iz))
  else
    pvap_liq = 1.01325e6 * 10._f ** (11.5349_f - 7972.9_f/t(iz))
  endif

  ! OLD: Saturation vapor pressure of orthorhomic sulfur (solid) from Lyons 2008

  !pvap_liq = 1.01325e6 * 10._f**(8.5121_f - 5102._f/t(iz)) ! T < 500 K
  !pvap_liq = 1.01325e6 * 10._f**(8.7832_f - 5166._f/t(iz)) 
  pvap_ice = pvap_liq

  return
end
