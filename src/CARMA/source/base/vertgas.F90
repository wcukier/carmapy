! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!PETER
!! This routine calculates vertrical transport rates for gases.
!! Currently treats diffusion only.
!! Not necessarily generalized for irregular grid.
!!
!! @author Peter Gao
!! @version Oct-2011
subroutine vertgas(carma, cstate, cvert, itbnd, ibbnd, ftop, fbot, cvert_tbnd, cvert_bbnd, vtrans, rc, igas) 

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

  type(carma_type), intent(in)         :: carma           !! the carma object
  type(carmastate_type), intent(inout) :: cstate          !! the carma state object
  real(kind=f), intent(inout)          :: cvert(NZ)       !! quantity being transported
  integer, intent(in)                  :: itbnd           !! top boundary condition
  integer, intent(in)                  :: ibbnd           !! bottom boundary condition
  real(kind=f), intent(in)             :: ftop            !! flux at top boundary
  real(kind=f), intent(in)             :: fbot            !! flux at bottom boundary
  real(kind=f), intent(in)             :: cvert_tbnd      !! quantity at top boundary
  real(kind=f), intent(in)             :: cvert_bbnd      !! quantity at bottom boundary
  real(kind=f), intent(in)             :: vtrans(NZP1)    !! vertical velocity
  integer, intent(inout)               :: rc              !! return code, negative indicates failure
  integer, intent(inout)               :: igas            !! gas index for versub call

  
  ! Local Variables
  integer      :: k
  integer      :: nzm1
  integer      :: itwo
  integer      :: iz
  integer      :: istep
  integer      :: nstep_sed
  real(kind=f) :: dz_avg
  real(kind=f) :: mixratiofactg
  real(kind=f) :: rhofactg
  real(kind=f) :: xexg
  real(kind=f) :: vertdifug(NZP1)
  real(kind=f) :: vertdifdg(NZP1)
  real(kind=f) :: vertadvug(NZP1)
  real(kind=f) :: vertadvdg(NZP1)
  real(kind=f) :: ttheta
  real(kind=f) :: cvert_bbnd_aug
  real(kind=f) :: cvert_tbnd_aug
  real(kind=f) :: fvert(NZ)
  real(kind=f) :: upg(NZP1)
  real(kind=f) :: dng(NZP1)
  real(kind=f) :: cfl_max
  real(kind=f) :: fvert_1
  real(kind=f) :: fvert_nz
  real(kind=f) :: eddydif         !PETER

  ! Set some constants
  nzm1 = max( 1, NZ-1 )
  itwo = min( 2, NZ   )

  vertdifug(:) = 0._f
  vertdifdg(:) = 0._f
  vertadvug(:) = 0._f
  vertadvdg(:) = 0._f

  !  Loop over vertical levels.
  do k = 2, NZ

    dz_avg = dz(k)                            ! layer thickness

    !  Check the vertical coordinate
    if( (igridv .eq. I_CART) .or. ( igridv .eq. I_LOGP)) then !DPOW check logp later
    
    !write(*,*) ekz(k), t(k)/t(1), 1.e9_f/(t(k)/t(1))**2._f
     
     !rhofactg = log(  rhoa(k)/rhoa(k-1) &
     !               * zmet(k-1)/zmet(k) )
     rhofactg = log(  rhoa(k)/rhoa(k-1) )
      !xexg = rhoa(k-1)/rhoa(k) * &
       !     zmet(k)/zmet(k-1)
      xexg = rhoa(k-1)/rhoa(k) 
	  
	  if( abs(ONE - xexg) .lt. ALMOST_ZERO ) xexg = ALMOST_ONE
	  
      vertdifug(k) = ( rhofactg * ekz(k) / dz_avg ) / ( 1._f - xexg )
      vertdifdg(k) = vertdifug(k) * xexg
      
      !write(*,*) rhofactg, xexg, log(  rhoa(k)/rhoa(k-1)), rhoa(k-1)/rhoa(k)
      !write(*,*) xexg, ( rhofactg ) / ( 1._f - xexg ), ekz(k) / dz_avg, ( rhofactg * ekz(k) / dz_avg ) / ( 1._f - xexg )

	!  ...else you're in sigma or hybrid coordinates...
    elseif(( igridv .eq. I_SIG ) .or. ( igridv .eq. I_HYBRID )) then
      vertdifug(k) = ekz(k) / dz_avg
      vertdifdg(k) = ekz(k) / dz_avg

    !  ...else write an error (maybe redundant)...
    else
      if (do_print) write(LUNOPRT,*) 'vertgas::ERROR - Invalid vertical grid type (', igridv, ').'
      rc = -1
      return
    endif
  enddo

  ! Fluxes at boundaries specified by user
  if( ibbnd .eq. I_FLUX_SPEC ) then
    vertdifug(1) = 0._f
    vertdifdg(1) = 0._f
  endif

  if( itbnd .eq. I_FLUX_SPEC ) then
    vertdifug(NZ+1) = 0._f
    vertdifdg(NZ+1) = 0._f
  endif

  if( (ibbnd .eq. I_FIXED_CONC) .or. (ibbnd .eq. I_ZERO_CGRAD) ) then
  
  	if(( igridv .eq. I_SIG ) .or. ( igridv .eq. I_HYBRID )) then
  		dz_avg = dz(NZ) 
  		vertdifug(NZ+1) = ekz(NZ+1) / dz_avg
		vertdifdg(NZ+1) = ekz(NZ+1) / dz_avg
		
		!write(*,*) vertdifug(NZ+1), vertdifdg(NZ+1), 'bottom'
  	else
    	dz_avg = dz(1)                            ! layer thickness
    	!dz_avg = zc(1) - zl(1)
    	!rhofactg = log( p(1)/pl(1) )
    	!rhofactg = log( rhoa(itwo)/rhoa(1) &
         !           * zmet(1)/zmet(itwo) )
        rhofactg = log( rhoa(itwo)/rhoa(1))
                
    	ttheta = rhofactg
    	if( ttheta .ge. 0._f ) then
    	  ttheta = min(ttheta,POWMAX)
    	else
    	  ttheta = max(ttheta,-POWMAX)
    	endif
	
    	xexg = exp(-ttheta)
    	if( abs(ONE - xexg) .lt. ALMOST_ZERO ) xexg = ALMOST_ONE
    	vertdifug(1) = ( rhofactg * ekz(1) / dz_avg ) / ( 1._f - xexg )
    	vertdifdg(1) = vertdifug(1) * xexg
    endif
  endif

  if( (itbnd .eq. I_FIXED_CONC) .or. (itbnd .eq. I_ZERO_CGRAD) ) then
  	if(( igridv .eq. I_SIG ) .or. ( igridv .eq. I_HYBRID )) then
  		dz_avg = dz(1) 
  		
  		vertdifug(1) = ekz(1) / dz_avg
		vertdifdg(1) = ekz(1) / dz_avg
  	else
    	dz_avg = dz(NZ)                            ! layer thickness
    	!rhofactg = log( rhoa(NZ)/rhoa(nzm1) &
        !            * zmet(nzm1)/zmet(NZ) )
        rhofactg = log( rhoa(NZ)/rhoa(nzm1))
    	!dz_avg = zl(NZP1) - zc(NZ)                            ! layer thickness
    	!rhofactg = log( pl(NZP1)/p(NZ) )
    	ttheta = rhofactg
    	if( ttheta .ge. 0._f ) then
    	  ttheta = min(ttheta,POWMAX)
    	else
    	  ttheta = max(ttheta,-POWMAX)
    	endif
	
    	xexg = exp(-ttheta)
    	if( abs(ONE - xexg) .lt. ALMOST_ZERO ) xexg = ALMOST_ONE
    	vertdifug(NZ+1) = ( rhofactg * ekz(NZ+1) / dz_avg ) / ( 1._f - xexg )
    	vertdifdg(NZ+1) = vertdifug(NZ+1) * xexg
    endif
  endif
        !if (igas .eq. 3) then
        !do k = 1,NZ+1
        !write(*,*) 'k,vertdifug,vertdifdg= ', k, vertdifug(k), vertdifdg(k)
        !enddo
        !do k = 1,NZ
        !write(*,*) 'k,cvert,rhoa,mix=',k,cvert(k),rhoa(k),cvert(k)/rhoa(k)*WTMOL_AIR/WTMOL_ZNS
        !enddo
        !endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call vertadv(carma, cstate, vtrans, cvert, itbnd, ibbnd, &
	cvert_tbnd, cvert_bbnd, vertadvug, vertadvdg, rc)
   if (rc < RC_OK) return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (do_explised) then     !PETER

    call versub(carma, cstate, pconmax(:,1)*xmet(:)*ymet(:)*zmet(:), &
		cvert, itbnd, ibbnd, ftop, fbot, cvert_tbnd, cvert_bbnd, &
		vertadvug, vertadvdg, vertdifug, vertdifdg, rc, 0, igas, 1, 1)  								    !PETER
    if (rc < RC_OK) return                                                                          !PETER

  else                                                            !PETER

    call versol(carma, cstate, cvert, itbnd, ibbnd, &             !PETER
                ftop, fbot, cvert_tbnd, cvert_bbnd, &             !PETER
                vertadvug, vertadvdg, vertdifug, vertdifdg, rc, 0, igas, 1, 1)   !PETER
    if (rc < RC_OK) return                                        !PETER

  end if                                                          !PETER

  return
end
