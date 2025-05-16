! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine drives the potentially slower microphysics calculations.
!!
!!  Originally part of microphy.  Now in this separate routine to allow
!!  time splitting of coagulation at a different timestep size from
!!  other microphysical calcs.
!!
!! @author McKie
!! @version Sep-1997
subroutine microslow(carma, cstate, rc)

  ! carma types defs
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
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local Declarations
  integer   :: i		!PETER
  integer   :: iz		!PETER
  integer   :: ibin
  integer   :: ielem
  integer   :: elemmultiple	!PETER

  if (do_printdiag) then										!PETER
    write(lundiag,*) 'PART 4: COAGULATION'								!PETER
    write(lundiag,*) '****************************************'						!PETER
    write(lundiag,*) 'Z: ALTITUDE LEVEL INDEX'								!PETER
    write(lundiag,*) 'BIN: SIZE BIN'									!PETER
    write(lundiag,*) 'COAGPROD X: PRODUCTION RATE OF ELEM X DUE TO COAGULATION [#/cm3/s]'		!PETER
    write(lundiag,*) 'COAGLOSS X: LOSS RATE OF ELEM X DUE TO COAGULATION [#/cm3/s]'			!PETER
    write(lundiag,*) 'NETFLUX X: NET PRODUCTION RATE OF ELEM X DUE TO COAGULATION [#/cm3/s]'		!PETER
    write(lundiag,*) '****************************************'						!PETER 
    elemmultiple = int(NELEM / 3._f)									!PETER
    coagprod(:,:,:) = 0._f										!PETER
    coagloss(:,:,:) = 0._f										!PETER
  end if												!PETER



  ! Calculate (implicit) particle loss rates for coagulation.
  call coagl(carma, cstate, rc)

  ! Calculate particle production terms and solve for particle 
  ! concentrations at end of time step.
  !
  ! NOTE: The order of elements required by CARMA to work with the
  ! element loop first is: if you have a group that is both a source
  ! and product of coagulation, then it needs to come after the
  ! other group that participates in that coagulation in the element
  ! table. For example, icoag(2,1) = 1 will not work, but
  ! icoag(2,1) = 2 should work.
  do ielem = 1,NELEM
    do ibin = 1,NBIN
      call coagp(carma, cstate, ibin, ielem, rc)
      call csolve(carma, cstate, ibin, ielem, rc)
    enddo
  enddo


  if (do_printdiag) then 												!PETER
  do i = 1, elemmultiple												!PETER
    write(lundiag,'(2A5,$)') 'Z', 'BIN'											!PETER
    do ielem = 3*(i-1)+1,3*i-1												!PETER
      write(lundiag,'(A12,i2,A12,i2,A12,i2,$)') 'COAGPROD',ielem,'COAGLOSS',ielem,'NETFLUX',ielem			!PETER
    end do														!PETER
    write(lundiag,'(A12,i2,A12,i2,A12,i2)') 'COAGPROD',3*i,'COAGLOSS',3*i,'NETFLUX',3*i					!PETER
    do iz = 1, NZ													!PETER
      do ibin = 1, NBIN													!PETER
        write(lundiag,'(2i5,$)') iz, ibin										!PETER
        do ielem = 3*(i-1)+1,3*i-1											!PETER
          write(lundiag,'(3e14.3,$)') coagprod(iz,ibin,ielem), coagloss(iz,ibin,ielem), &				!PETER
						      coagprod(iz,ibin,ielem)-coagloss(iz,ibin,ielem)			!PETER
        end do														!PETER
        write(lundiag,'(3e14.3)') coagprod(iz,ibin,3*i), coagloss(iz,ibin,3*i), &					!PETER
						    coagprod(iz,ibin,3*i)-coagloss(iz,ibin,3*i)				!PETER
      end do														!PETER
    end do														!PETER
    write(lundiag,*) ' ' 												!PETER
  end do														!PETER
  if (3*elemmultiple+1 .le. NELEM) then											!PETER
    write(lundiag,'(2A5,$)') 'Z', 'BIN'											!PETER
    do ielem = 3*elemmultiple+1,NELEM-1											!PETER
      write(lundiag,'(A12,i2,A12,i2,A12,i2,$)') 'COAGPROD',ielem,'COAGLOSS',ielem,'NETFLUX',ielem			!PETER
    end do														!PETER
    write(lundiag,'(A12,i2,A12,i2,A12,i2,$)') 'COAGPROD',NELEM,'COAGLOSS',NELEM,'NETFLUX',NELEM
    do iz = 1, NZ													!PETER
      do ibin = 1, NBIN													!PETER
        write(lundiag,'(2i5,$)') iz, ibin										!PETER
        do ielem = 3*elemmultiple+1,NELEM-1										!PETER
          write(lundiag,'(3e14.3,$)') coagprod(iz,ibin,ielem), coagloss(iz,ibin,ielem), &				!PETER
						      coagprod(iz,ibin,ielem)-coagloss(iz,ibin,ielem)			!PETER
        end do														!PETER
        write(lundiag,'(3e14.3)') coagprod(iz,ibin,NELEM), coagloss(iz,ibin,NELEM), &					!PETER
						    coagprod(iz,ibin,NELEM)-coagloss(iz,ibin,NELEM)			!PETER
      end do														!PETER
    end do														!PETER
  end if														!PETER
  end if

  if (do_printdiag) write(lundiag,*) ' '		!PETER
  
  ! Return to caller with new particle concentrations.
  return
end
