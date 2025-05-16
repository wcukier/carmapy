! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine applies evaporation and nucleation production terms to
!! particle concentrations.
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine downgevapply(carma, cstate, iz, rc)

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
  integer, intent(in)                  :: iz      !! z index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                              :: ibin    !! bin index
  integer                              :: ielem   !! element index


  ! Visit each radius bin for each element to compute particle production 
  ! due to evaporation and element transfer processes for which the source
  ! element number is greater than the target element number
  do ielem = 1,NELEM
    do ibin = 1,NBIN

	!write(*,*) 'downevapply before', iz, ibin, ielem, pc(iz,ibin,ielem), evappe(ibin,ielem), rnucpe(ibin,ielem)

      pc(iz,ibin,ielem) = pc(iz,ibin,ielem) + &
                         dtime * ( evappe(ibin,ielem) + &
                                   rnucpe(ibin,ielem) )

	!write(*,*) 'downevapply after', iz, ibin, ielem, pc(iz,ibin,ielem)

     ! write(*,*) 'downevapply', iz, ibin, ielem, evappe(ibin,ielem), rnucpe(ibin,ielem) 
      !write(*,*) 'downevapply', iz, ibin, ielem, evappe(ibin,ielem)
                                   
      ! Prevent particle concentrations from dropping below SMALL_PC
      call smallconc(carma, cstate, iz, ibin, ielem, rc)

    enddo
  enddo

  ! Return to caller with evaporation and down-grid element transfer
  ! production terms applied to particle concentrations.
  return
end
