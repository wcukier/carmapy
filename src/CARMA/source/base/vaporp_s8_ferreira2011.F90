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
subroutine vaporp_s8_ferreira2011(carma, cstate, iz, rc, pvap_liq, pvap_ice)
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

  ! Saturation vapor presure of (1) liquid, (2) (ortho-)rhombic, and (3) monoclinic sulfur from 
  ! Ferreira & Lobo 2011. 
  if (t(iz) .gt. 424.79_f) then 
    pvap_liq = 10._f * exp(-1.4373_f - 4900.7_f/t(iz) + 3.0068_f * log(t(iz)) &
      - 5.229e-7_f*t(iz))  ! Ev(T)*p/T term neglected due to inability to solve such an equation
  else if (t(iz) .lt. 368.39_f) then
    pvap_liq = 10._f * exp(23.641_f - 12695.8_f/t(iz) + 3.124_f * log(t(iz)) &
      - 0.0366*t(iz) + 5.093e-5_f*t(iz)**2._f - 3.574e-8_f*t(iz)**3._f)
  else
    pvap_liq = 10._f * exp(90.604_f - 13332.4_f/t(iz) - 10.955_f * log(t(iz)) &
      + 0.0366*t(iz) - 3.242e-5_f*t(iz)**2._f + 1.079e-8_f*t(iz)**3._f)
  endif

  ! OLD: Saturation vapor pressure of orthorhomic sulfur (solid) from Lyons 2008

  !pvap_liq = 1.01325e6 * 10._f**(8.5121_f - 5102._f/t(iz)) ! T < 500 K
  !pvap_liq = 1.01325e6 * 10._f**(8.7832_f - 5166._f/t(iz)) 
  pvap_ice = pvap_liq

  return
end
