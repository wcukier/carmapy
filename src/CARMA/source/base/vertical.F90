! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine drives the vertical transport calculations.
!!
!!  NOTE: Since this is only for sedimentation and brownian diffusion of a column within
!! a parent model, the advection of air density, gases and potential temperature have
!! been removed. Also, the divergence corrections (divcor) for 1D transport are not
!! applied, since these columns exist within a parent model that is responsible for the
!! advection.
!!
!! @author Eric Jensen
!! version Mar-1995
!! 
!! WC -- Note, I removed a lot of diagnostic print statements from this file
subroutine vertical(carma, cstate, rc)

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
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Declare local variables
  integer        :: ielem
  integer        :: i			!PETER
  integer        :: ibin
  integer        :: ig
  integer        :: igas                   !PETER
  integer        :: iz                     !PETER
  integer        :: nzm1                   !PETER
  integer        :: elemmultiple           !PETER
  integer        :: gasmultiple           !PETER
  integer        :: gas_switch             !PETER
  real(kind=f)   :: vertadvu(NZP1)
  real(kind=f)   :: vertadvd(NZP1)
  real(kind=f)   :: vertdifu(NZP1)
  real(kind=f)   :: vertdifd(NZP1)
  real(kind=f)   :: vtrans(NZP1)
  real(kind=f)   :: vtransg(NZP1)          !PETER
  real(kind=f)   :: netflux(NZ)          !PETER
  real(kind=f)   :: vertin_sum(NZ,NBIN,NELEM)          !PETER
  real(kind=f)   :: vertout_sum(NZ,NBIN,NELEM)          !PETER
  real(kind=f)   :: old_pc(NZ)

  rc = RC_OK

  vertadvu(:) = 0._f
  vertadvd(:) = 0._f
  vertdifu(:) = 0._f
  vertdifd(:) = 0._f
 
  do ielem = 1,NELEM          ! Loop over particle elements
    ig = igelem(ielem)        ! particle group
    
    ! Should this group participate in sedimentation?
    if (grp_do_vtran(ig)) then

        do ibin = 1,NBIN          ! Loop over particle mass bins          
          vtrans(:) = -vf(:,ibin,ig)
    
          ! If dry deposition is enabled for this group, then set
          ! the deposition velocity at the surface.
          if (grp_do_drydep(ig)) then
            if ((igridv .eq. I_CART) .or. (igridv .eq. I_LOGP)) then
              vtrans(1) = -vd(ibin, ig)
            else
              vtrans(NZP1) = -vd(ibin, ig)
            end if
          end if

          !  Calculate particle transport rates due to vertical advection
          !  and vertical diffusion, and solve for concentrations at end of time step.
          call vertadv(carma, cstate, vtrans, pc(:,ibin,ielem), itbnd_pc, ibbnd_pc, &
            pc_topbnd(ibin,ielem), pc_botbnd(ibin,ielem), vertadvu, vertadvd, rc)
          if (rc < RC_OK) return

          call vertdif(carma, cstate, ig, ibin, itbnd_pc, ibbnd_pc, vertdifu, vertdifd, rc)
          if (rc < RC_OK) return

          old_pc(:) = pc(:,ibin,ielem)
          
          ! There are 2 different solvers, versol with uses a PPM scheme and versub
          ! which using an explicit substepping approach.
          if (do_explised) then
            call versub(carma, cstate, pconmax(:,ig)*xmet(:)*ymet(:)*zmet(:), &
		          pc(:,ibin,ielem), itbnd_pc, ibbnd_pc, &
              ftoppart(ibin,ielem), fbotpart(ibin,ielem), &
              pc_topbnd(ibin,ielem), pc_botbnd(ibin,ielem), &
              vertadvu, vertadvd, vertdifu, vertdifd, rc, 1, 1, ibin, ielem)                        !PETER
            if (rc < RC_OK) return
          else
            call versol(carma, cstate, pc(:,ibin,ielem), itbnd_pc, ibbnd_pc, &
              ftoppart(ibin,ielem), fbotpart(ibin,ielem), &
              pc_topbnd(ibin,ielem), pc_botbnd(ibin,ielem), &
              vertadvu, vertadvd, vertdifu, vertdifd, rc, 1, 1, ibin, ielem)
            if (rc < RC_OK) return
          end if
          
          ! A clunky way to get the mass flux to the surface and to conserve mass
          ! is to determine the total before and after. Anything lost went to the
          ! surface.
          !
          ! NOTE: This only works if you assume nothing is lost out the top. It would be
          ! better to figure out how to get this directly from versol.
          pc_surf(ibin,ielem) = pc_surf(ibin, ielem) + &
		        sum(old_pc(:) * dz(:) / xmet(:) / ymet(:)) - &
            sum(pc(:,ibin,ielem) * dz(:) / xmet(:) / ymet(:))
          sedimentationflux(ibin,ielem) = ( sum(old_pc(:) * dz(:) / xmet(:) / ymet(:)) - &
            sum(pc(:,ibin,ielem) * dz(:) / xmet(:) / ymet(:)) ) / dtime
        enddo  ! ibin
      !endif
    endif
  enddo  ! ielem

  do iz = 1,NZ                                                                  !PETER
    vtransg(iz) = winds(iz)							                                        !PETER
  enddo   								                                                      !PETER

  vtransg(NZP1) = vtransg(NZ) 						                                      !PETER

  nzm1 = max(1, NZ-1)                                                           !PETER

  if (NZ .gt. 1) then                                                           !PETER
    vtransg(NZ) = (vtransg(nzm1) + vtransg(NZ)) / 2._f                          !PETER
 
    if (NZ .gt. 2) then                                                         !PETER
      do iz = NZ-1, 2, -1                                                       !PETER
        vtransg(iz) = (vtransg(iz-1) + vtransg(iz)) / 2._f                      !PETER
      enddo                                                                     !PETER
    endif                                                                       !PETER
  endif                                                                         !PETER

  do igas = 1,NGAS                                                              !PETER    
    call vertgas(carma, cstate, gc(:,igas), itbnd_gc, ibbnd_gc, &               !PETER
              ftopgas(igas), fbotgas(igas), &                                   !PETER
              gc_topbnd(igas), gc_botbnd(igas), vtransg, rc, igas)              !PETER
    if (rc < RC_OK) return                                                      !PETER
  enddo                                                                         !PETER 

  ! Return to caller with new particle concentrations.
  return
end
