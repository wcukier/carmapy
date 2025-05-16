! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine solves for sedimentation using an explicit substepping approach. It
!!  is faster and handles large cfl and irregular grids better than the normal PPM
!!  solver (versol), but it is more diffusive.
!!
!! @author Andy Ackerman, Chuck Bardeen 
!! version Aug 2010
subroutine versub(carma, cstate, pcmax, cvert, itbnd, ibbnd, ftop, fbot, cvert_tbnd, cvert_bbnd, &
  vertadvu, vertadvd, vertdifu, vertdifd, rc, gas_switch, igas, ibin, ielem)						!PETER

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
  real(kind=f), intent(in)             :: pcmax(NZ)       !! maximum particle concentration (#/x/y/z)
  real(kind=f), intent(inout)          :: cvert(NZ)       !! quantity being transported (#/x/y/z)
  integer, intent(in)                  :: itbnd           !! top boundary condition
  integer, intent(in)                  :: ibbnd           !! bottom boundary condition
  real(kind=f), intent(in)             :: ftop            !! flux at top boundary
  real(kind=f), intent(in)             :: fbot            !! flux at bottom boundary
  real(kind=f), intent(in)             :: cvert_tbnd      !! quantity at top boundary
  real(kind=f), intent(in)             :: cvert_bbnd      !! quantity at bottom boundary
  real(kind=f), intent(in)             :: vertadvu(NZP1)  !! upward vertical transport rate into level k from level k-1 [cm/s]
  real(kind=f), intent(in)             :: vertadvd(NZP1)  !! downward vertical transport rate into level k from level k-1 [cm/s]
  real(kind=f), intent(in)             :: vertdifu(NZP1)  !! upward vertical diffusion rate into level k from level k-1 [cm/s]
  real(kind=f), intent(in)             :: vertdifd(NZP1)  !! downward vertical diffusion rate into level k from level k-1 [cm/s]
  integer, intent(inout)               :: rc              !! return code, negative indicates failure
  integer, optional, intent(in)        :: gas_switch      !! PETER: 0 - being used as part of vertgas, otherwise it's used for particles
  integer, optional, intent(in)        :: igas            !! PETER: gas index for updating gflux
  integer, optional, intent(in)        :: ibin            !! PETER: bin index for updating pflux
  integer, optional, intent(in)        :: ielem           !! PETER: element index for updating pflux

  !  Declare local variables
  integer                        :: iz
  integer                        :: istep
  integer                        :: nstep_sed
  real(kind=f)                   :: fvert(NZ)
  real(kind=f)                   :: fvert_sum(NZP1)     !PETER
  real(kind=f)                   :: boxflux(NZP1)     !PETER
  real(kind=f)                   :: vertupin(NZ)     !PETER
  real(kind=f)                   :: vertupout(NZ)     !PETER
  real(kind=f)                   :: vertdnin(NZ)     !PETER
  real(kind=f)                   :: vertdnout(NZ)     !PETER
!  real(kind=f)                   :: vertin_sum(NZ)     !PETER
!  real(kind=f)                   :: vertout_sum(NZ)     !PETER
  real(kind=f)                   :: up(NZP1)
  real(kind=f)                   :: dn(NZP1)
  real(kind=f)                   :: cfl_max
  real(kind=f)                   :: fvert_1
  real(kind=f)                   :: fvert_nz

  ! Determine the total upward and downward velocities.
  up(:) = vertadvu(:) + vertdifu(:)
  dn(:) = vertadvd(:) + vertdifd(:)

  fvert_sum(:) = 0._f	!PETER
  vertupin(:) = 0._f	!PETER
  vertupout(:) = 0._f	!PETER
  vertdnin(:) = 0._f	!PETER
  vertdnout(:) = 0._f	!PETER
!  vertin_sum(:) = 0._f	!PETER
!  vertout_sum(:) = 0._f	!PETER

  ! Compute the maximum CFL for each bin that has a significant concentration
  ! of particles.
  cfl_max = 0._f

  do iz = 1, NZ
    if (gas_switch .eq. 0) then											!PETER
      cfl_max = max(cfl_max, max(abs(up(iz)), abs(up(iz+1)), abs(dn(iz)), abs(dn(iz+1))) * dtime / dz(iz))	!PETER
    else													!PETER
      if (pcmax(iz) > SMALL_PC) then
        cfl_max = max(cfl_max, max(abs(up(iz)), abs(up(iz+1)), abs(dn(iz)), abs(dn(iz+1))) * dtime / dz(iz))
      end if
    endif													!PETER
  end do

  ! Use the maximum CFL determined above to figure out how much substepping is
  ! needed to sediment explicitly without violating the CFL anywhere in the column.
  if (cfl_max >= 0._f) then
    nstep_sed = int(1._f + cfl_max)
  else
    nstep_sed = 1
  endif
  
  ! If velocities are in both directions, then more steps are needed to make sure
  ! that no more than half of the concentration can be transported in either direction.
  if (maxval(up(:) * dn(:)) > 0._f) then
    nstep_sed = nstep_sed * 2
  end if

  ! Determine the top and bottom boundary fluxes, keeping in mind that
  ! the velocities and grid coordinates are reversed in sigma or hybrid
  ! coordinates
  if ((igridv .eq. I_SIG) .or. (igridv .eq. I_HYBRID)) then
    if (itbnd .eq. I_FLUX_SPEC) then
      fvert_nz = -fbot
    elseif (itbnd .eq. I_ZERO_CGRAD) then !DIANA
      fvert_nz = cvert(1)*up(1)
    else
      fvert_nz = cvert_bbnd*dn(NZ+1)
    end if

    if (ibbnd .eq. I_FLUX_SPEC) then
      fvert_1 = -ftop
    elseif (ibbnd .eq. I_ZERO_CGRAD) then !DIANA
      fvert_1 = cvert(NZ)*dn(NZ+1)
    else
      fvert_1 = cvert_tbnd*up(1)
    end if
  
  else
    if (itbnd .eq. I_FLUX_SPEC) then
      fvert_nz = ftop
    elseif (itbnd .eq. I_ZERO_CGRAD) then
      fvert_nz = cvert(NZ)*dn(NZ+1)
    else 
      fvert_nz = cvert_tbnd*dn(NZ+1)
    end if

    if (ibbnd .eq. I_FLUX_SPEC) then
      fvert_1 = fbot
    elseif (ibbnd .eq. I_ZERO_CGRAD) then
      fvert_1 = cvert(1)*up(1)
    else
      fvert_1 = cvert_bbnd*up(1)
    end if
  endif

  ! Sediment the particles using multiple iterations to satisfy the CFL.
  do istep = 1, nstep_sed

    ! Determine the net particle flux at each gridbox. The first and last levels
    ! need special treatment to handle to bottom and top boundary conditions.
    fvert(1) = (-cvert(1)*dn(1) + fvert_1 + cvert(2)*dn(2) - cvert(1)*up(2))
    boxflux(1) = -cvert(1)*dn(1) + fvert_1 					!PETER
    vertupin(1) = cvert(2)*dn(2)					!PETER
    vertupout(1) = cvert(1)*up(2)				!PETER
    vertdnin(1) = fvert_1 					!PETER
    vertdnout(1) = cvert(1)*dn(1) 				!PETER
    
    do iz = 2, NZ-1
      fvert(iz) = (-cvert(iz)*dn(iz) + cvert(iz-1)*up(iz) + cvert(iz+1)*dn(iz+1) - cvert(iz)*up(iz+1))
      vertupin(iz) = cvert(iz+1)*dn(iz+1)					!PETER
      vertupout(iz) = cvert(iz)*up(iz+1)					!PETER
      vertdnin(iz) = cvert(iz-1)*up(iz)					!PETER
      vertdnout(iz) = cvert(iz)*dn(iz)					!PETER
    end do

    do iz = 2, NZ								!PETER
      boxflux(iz) = -cvert(iz)*dn(iz) + cvert(iz-1)*up(iz)			!PETER
    enddo									!PETER

    fvert(NZ) = (-cvert(NZ)*dn(NZ) + cvert(NZ-1)*up(NZ) + &
	fvert_nz - cvert(NZ)*up(NZ+1))
    boxflux(NZP1) = -fvert_nz + cvert(NZ)*up(NZP1)						!PETER
    vertupin(NZ) = fvert_nz							!PETER
    vertupout(NZ) = cvert(NZ)*up(NZ+1)						!PETER
    vertdnin(NZ) = cvert(NZ-1)*up(NZ)							!PETER
    vertdnout(NZ) = cvert(NZ)*dn(NZ)						!PETER
    
    ! Now update the actual concentrations.
    cvert(:) = cvert(:) + fvert(:) * dtime / nstep_sed / dz(:)
    
    fvert_sum(:) = fvert_sum(:) + boxflux(:) / nstep_sed                !PETER
    if (gas_switch .eq. 0) then									!PETER
      vertupin_sum(:,igas,1) = vertupin_sum(:,igas,1) + vertupin(:) / nstep_sed / dz(:)		!PETER
      vertupout_sum(:,igas,1) = vertupout_sum(:,igas,1) + vertupout(:) / nstep_sed / dz(:)		!PETER
      vertdnin_sum(:,igas,1) = vertdnin_sum(:,igas,1) + vertdnin(:) / nstep_sed / dz(:)		!PETER
      vertdnout_sum(:,igas,1) = vertdnout_sum(:,igas,1) + vertdnout(:) / nstep_sed / dz(:)		!PETER
    else											!PETER
      vertupin_sum(:,ibin,ielem) = vertupin_sum(:,ibin,ielem) + vertupin(:) / nstep_sed / dz(:)	!PETER
      vertupout_sum(:,ibin,ielem) = vertupout_sum(:,ibin,ielem) + vertupout(:) / nstep_sed / dz(:)	!PETER
      vertdnin_sum(:,ibin,ielem) = vertdnin_sum(:,ibin,ielem) + vertdnin(:) / nstep_sed / dz(:)	!PETER
      vertdnout_sum(:,ibin,ielem) = vertdnout_sum(:,ibin,ielem) + vertdnout(:) / nstep_sed / dz(:)	!PETER
    endif 											!PETER
  enddo
   
  if (gas_switch .eq. 0) then		 !PETER
    gflux(:,igas) = fvert_sum(:)         !PETER
  else			        	 !PETER
    pflux(:,ibin,ielem) = fvert_sum(:)   !PETER 
  endif					 !PETER

!  if (gas_switch .eq. 0) then		 											!PETER
!    if (do_printdiag) then													!PETER
!      do iz = 1, NZ														!PETER
!        write(lundiag,'(2I5,I15,3E11.3)') iz, 0, igas, vertin_sum(iz), vertout_sum(iz), vertin_sum(iz)-vertout_sum(iz)          !PETER
!      end do															!PETER
!    endif															!PETER
!  else			        	 											!PETER
!    if (do_printdiag) then 								        				!PETER
!      do iz = 1, NZ														!PETER
!        write(lundiag,'(2I5,I15,3E11.3)') iz, ibin, ielem, vertin_sum(iz), vertout_sum(iz), vertin_sum(iz)-vertout_sum(iz)	!PETER
!      end do															!PETER
!    endif															!PETER
!  endif					 											!PETER
  
  return
end subroutine versub
