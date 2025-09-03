! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates derived mapping arrays and calculates the critical
!! supersaturation <scrit> used to nucleate dry particles (CN) to droplets.
!!
!! This routine requires that array <akelvin> is defined.
!! (i.e., setupgkern.f must be called before this)
!!
!! NOTE: Most of the code from this routine has been moced to CARMA_InitializeGrowth
!! because it does not rely upon the model's state and thus can be called one during
!! CARMA_Initialize rather than being called every timestep if left in this routine.
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine setupnuc(carma, cstate, rc)

  ! types
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

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc       !! return code, negative indicates failure

  ! Local declarations
  integer                        :: igroup   ! group index
  integer                        :: igas     ! gas index
  integer                        :: ielemfrom ! index of nucleating element
  integer                        :: ielem     ! element index
  integer                        :: isol     ! solute index
  integer                        :: ibin     ! bin index
  integer                        :: inuc     ! nucleation index
  integer                        :: k        ! z index
  real(kind=f)                   :: bsol

  integer                        :: i

  !write(*,*) mucos

  ! Define formats
  3 format(a,a)
  6 format(i4,5x,1p2e11.3)
  8 format(/,'Critical supersaturations for ',a,//, '   i        r [cm]     scrit',/)

  ! Define critical supersaturation and target bin for each (dry) particle
  ! size bin that is subject to nucleation.
  ! (only for CN groups subject to nucleation)
  do igroup = 1,NGROUP

    ielemfrom = ienconc(igroup)
    do inuc = 1,nnuc2elem(ielemfrom)
      igas = inucgas(ielemfrom,inuc2elem(inuc, ielemfrom))

      if( igas .ne. 0 .and. (itype( ielemfrom ) .eq. &
        I_INVOLATILE .or. itype( ielemfrom ) .eq. I_VOLATILE) ) then

        isol = isolelem( ienconc( igroup ) )
      
        ! If here is no solute are specified, then no scrit value is defined.
        do k = 1,NZ
          if (isol .ne. 0) then
            ! This is term "B" in Pruppacher and Klett's eqn. 6-28.=
            do ibin = 1,NBIN
                bsol = 3._f*sol_ions(isol)*rmass(ibin,igroup)*gwtmol(igas) &
                  / ( 4._f*PI*solwtmol(isol)*RHO_W )
                scrit(k,ibin,igroup,igas) = sqrt( 4._f * akelvin(k,igas)**3 / ( 27._f * bsol ) )
	    enddo
	  endif


        enddo
      endif
    enddo
  enddo

  ! Return to caller with nucleation mapping arrays and critical
  ! supersaturations defined.
  return
end
