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
  real(kind=f)   :: aa					!PETER
  real(kind=f)   :: bb					!PETER
  real(kind=f)   :: cc					!PETER
  real(kind=f)   :: dd					!PETER
  real(kind=f)   :: ee					!PETER
  real(kind=f)   :: ff					!PETER
  
  rc = RC_OK

  vertadvu(:) = 0._f
  vertadvd(:) = 0._f
  vertdifu(:) = 0._f
  vertdifd(:) = 0._f

  if (do_printdiag) then							!PETER
    elemmultiple = int(NELEM / 3._f)						!PETER
    write(lundiag,*) 'PART 3: ADVECTION AND DIFFUSION'							!PETER
    write(lundiag,*) '****************************************'						!PETER
    write(lundiag,*) 'Z: ALTITUDE LEVEL INDEX'								!PETER
    write(lundiag,*) 'BIN: SIZE BIN'									!PETER
    write(lundiag,*) 'VERTUPIN: TOTAL PARTICLE MASS FLUX INTO Z FROM ABOVE [g/cm3/s]'		!PETER
    write(lundiag,*) 'VERTUPOUT: TOTAL PARTICLE MASS FLUX OUT OF Z FROM ABOVE [g/cm3/s]'		!PETER
    write(lundiag,*) 'VERTDNIN: TOTAL PARTICLE MASS FLUX INTO Z FROM BELOW [g/cm3/s]'		!PETER
    write(lundiag,*) 'VERTDNOUT: TOTAL PARTICLE MASS FLUX OUT OF Z FROM BELOW [g/cm3/s]'		!PETER
    write(lundiag,*) 'VERTINMASS: TOTAL PARTICLE MASS FLUX INTO Z [g/cm3/s]'		!PETER
    write(lundiag,*) 'VERTOUTMASS: TOTAL PARTICLE MASS FLUX OUT OF Z [g/cm3/s]'		!PETER
    write(lundiag,*) 'NETFLUXMASS: NET PARTICLE MASS FLUX INTO Z [g/cm3/s]'		!PETER
    write(lundiag,*) 'VERTIN X: FLUX OF ELEMENT X INTO BIN AT Z [#/cm3/s]'		!PETER
    write(lundiag,*) 'VERTOUT X: FLUX OF ELEMENT X OUT OF BIN AT Z [#/cm3/s]'		!PETER
    write(lundiag,*) 'NETFLUX X: NETFLUX OF ELEMENT X INTO BIN AT Z [#/cm3/s]'		!PETER
    write(lundiag,*) 'VERTUPING X: FLUX OF GAS X INTO Z FROM ABOVE [g/cm3/s]'			!PETER
    write(lundiag,*) 'VERTUPOUTG X: FLUX OF GAS X OUT OF Z FROM ABOVE [g/cm3/s]'		!PETER
    write(lundiag,*) 'VERTDNING X: FLUX OF GAS X INTO Z FROM BELOW [g/cm3/s]'			!PETER
    write(lundiag,*) 'VERTDNOUTG X: FLUX OF GAS X OUT OF Z FROM BELOW [g/cm3/s]'		!PETER
    write(lundiag,*) 'VERTING X: FLUX OF GAS X INTO Z [g/cm3/s]'			!PETER
    write(lundiag,*) 'VERTOUTG X: FLUX OF GAS X OUT OF Z [g/cm3/s]'		!PETER
    write(lundiag,*) 'NETFLUXG X: NETFLUX OF GAS X INTO Z [g/cm3/s]'		!PETER
    write(lundiag,*) 'VERTUPING: TOTAL FLUX OF GAS INTO Z FROM ABOVE [g/cm3/s]'			!PETER
    write(lundiag,*) 'VERTUPOUTG: TOTAL FLUX OF GAS OUT OF Z FROM ABOVE [g/cm3/s]'		!PETER
    write(lundiag,*) 'VERTDNING: TOTAL FLUX OF GAS INTO Z FROM BELOW [g/cm3/s]'			!PETER
    write(lundiag,*) 'VERTDNOUTG: TOTAL FLUX OF GAS OUT OF Z FROM BELOW [g/cm3/s]'		!PETER
    write(lundiag,*) 'VERTING: TOTAL TOTAL FLUX OF GAS INTO Z [g/cm3/s]'			!PETER
    write(lundiag,*) 'VERTOUTG: TOTAL FLUX OF GAS OUT OF Z [g/cm3/s]'		!PETER
    write(lundiag,*) 'NETFLUXG: NETFLUX OF GAS INTO Z [g/cm3/s]'		!PETER
    write(lundiag,*) '****************************************'			!PETER  
    vertupin_sum(:,:,:) = 0._f		!PETER
    vertupout_sum(:,:,:) = 0._f		!PETER
    vertdnin_sum(:,:,:) = 0._f		!PETER
    vertdnout_sum(:,:,:) = 0._f		!PETER
    aa = 0._f	!PETER
    bb = 0._f	!PETER
    cc = 0._f	!PETER
    dd = 0._f	!PETER
    ee = 0._f	!PETER
    ff = 0._f	!PETER 
  end if									!PETER


  
 
  do ielem = 1,NELEM          ! Loop over particle elements
    ig = igelem(ielem)        ! particle group
    
    ! Should this group participate in sedimentation?
    if (grp_do_vtran(ig)) then
    
      ! Are there enough particles in the column to bother?
    !  if (maxval(pconmax(:,ig)) .gt. FEW_PC) then

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
           ! write(*,*) 'vertadv', pc(:,ibin,ielem)
          call vertdif(carma, cstate, ig, ibin, itbnd_pc, ibbnd_pc, vertdifu, vertdifd, rc)
         ! call vertdif(carma, cstate, ig, ibin, itbnd_pc, ibbnd_pc, &
	!	pc(:,ibin,ielem), pc_topbnd(ibin,ielem), pc_botbnd(ibin,ielem), vertdifu, vertdifd, rc)      !PETER
          if (rc < RC_OK) return
            !write(*,*) 'vertdif', pc(:,ibin,ielem)
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
	  !  if (ielem .eq. 1) write(*,*) ibin, ielem, ftoppart(ibin,ielem)
            call versol(carma, cstate, pc(:,ibin,ielem), itbnd_pc, ibbnd_pc, &
              ftoppart(ibin,ielem), fbotpart(ibin,ielem), &
              pc_topbnd(ibin,ielem), pc_botbnd(ibin,ielem), &
              vertadvu, vertadvd, vertdifu, vertdifd, rc, 1, 1, ibin, ielem)
            if (rc < RC_OK) return
          end if
           ! write(*,*) 'versol', pc(:,ibin,ielem)
          
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

  if (do_printdiag) then									!PETER
  vertin_sum = vertupin_sum + vertdnin_sum							!PETER
  vertout_sum = vertupout_sum + vertdnout_sum							!PETER
  do iz = 1, NZ											!PETER
    netflux(iz) = sum(vertupin_sum(iz,:,1)*rmass(:,1)+&
			vertupin_sum(iz,:,2)*rmass(:,2)) + &			!PETER
		  sum(vertdnin_sum(iz,:,1)*rmass(:,1)+vertdnin_sum(iz,:,2)*rmass(:,2)) - &		!PETER
	 	  sum(vertupout_sum(iz,:,1)*rmass(:,1)+vertupout_sum(iz,:,2)*rmass(:,2)) - &		!PETER
		  sum(vertdnout_sum(iz,:,1)*rmass(:,1)+vertdnout_sum(iz,:,2)*rmass(:,2))		!PETER
  end do											!PETER
  write(lundiag,'(A5,7A15)') 'Z', 'VERTUPIN', 'VERTUPOUT', 'VERTDNIN', 'VERTDNOUT', &		!PETER
                             'VERTINMASS', 'VERTOUTMASS', 'NETFLUXMASS'				!PETER
  do iz = 1, NZ														!PETER
    write(lundiag,'(i5, 7e15.3)') iz, sum(vertupin_sum(iz,:,1)*rmass(:,1)+&
					vertupin_sum(iz,:,2)*rmass(:,2)), &		!PETER
				  sum(vertupout_sum(iz,:,1)*rmass(:,1)+vertupout_sum(iz,:,2)*rmass(:,2)), &		!PETER
	 			  sum(vertdnin_sum(iz,:,1)*rmass(:,1)+vertdnin_sum(iz,:,2)*rmass(:,2)), &		!PETER
				  sum(vertdnout_sum(iz,:,1)*rmass(:,1)+vertdnout_sum(iz,:,2)*rmass(:,2)), &		!PETER
				  sum(vertupin_sum(iz,:,1)*rmass(:,1)+vertupin_sum(iz,:,2)*rmass(:,2)) + &		!PETER
				  sum(vertdnin_sum(iz,:,1)*rmass(:,1)+vertdnin_sum(iz,:,2)*rmass(:,2)), &		!PETER
				  sum(vertupout_sum(iz,:,1)*rmass(:,1)+vertupout_sum(iz,:,2)*rmass(:,2)) + &		!PETER
				  sum(vertdnout_sum(iz,:,1)*rmass(:,1)+&
					vertdnout_sum(iz,:,2)*rmass(:,2)), netflux(iz)	!PETER
  end do														!PETER
  write(lundiag,*) ''													!PETER
  do iz = 1, NZ														!PETER
    aa = aa + sum(vertupin_sum(iz,:,1)*rmass(:,1)+&
	vertupin_sum(iz,:,2)*rmass(:,2))					!PETER
    bb = bb + sum(vertupout_sum(iz,:,1)*rmass(:,1)+&
	vertupout_sum(iz,:,2)*rmass(:,2))					!PETER
    cc = cc + sum(vertdnin_sum(iz,:,1)*rmass(:,1)+&
	vertdnin_sum(iz,:,2)*rmass(:,2))					!PETER
    dd = dd + sum(vertdnout_sum(iz,:,1)*rmass(:,1)+&
	vertdnout_sum(iz,:,2)*rmass(:,2))					!PETER
    ee = ee + sum(vertupin_sum(iz,:,1)*rmass(:,1)+&
	vertupin_sum(iz,:,2)*rmass(:,2)) + &					!PETER
	sum(vertdnin_sum(iz,:,1)*rmass(:,1)+&
	vertdnin_sum(iz,:,2)*rmass(:,2))						!PETER
    ff = ff + sum(vertupout_sum(iz,:,1)*rmass(:,1)+&
	vertupout_sum(iz,:,2)*rmass(:,2)) + &					!PETER
	sum(vertdnout_sum(iz,:,1)*rmass(:,1)+&
	vertdnout_sum(iz,:,2)*rmass(:,2))						!PETER
  enddo  														!PETER
  write(lundiag,'(A5,7e15.3)') 'TOT:', aa, bb, cc, dd, ee, ff, sum(netflux)							!PETER
  write(lundiag,*) ' '										!PETER
  write(lundiag,*) 'COLUMN PART FLUX = TOT * NZ * deltaz (in meters) * (100 cm / m)'		!PETER
  vertpartflux = sum(netflux)*NZ*200._f*100._f							!PETER
  write(lundiag,'(A28, e22.10)') 'COLUMN PART FLUX [g/cm2/s]: ', vertpartflux			!PETER
  write(lundiag,*) ' '										!PETER
  do i = 1, elemmultiple												!PETER
    write(lundiag,'(2A5,$)') 'Z', 'BIN'								!PETER
    do ielem = 3*(i-1)+1,3*i-1												!PETER
      write(lundiag,'(A12,i2,A12,i2,A12,i2,$)') 'VERTIN',ielem,'VERTOUT',ielem,'NETFLUX',ielem	!PETER
    end do														!PETER
    write(lundiag,'(A12,i2,A12,i2,A12,i2)') 'VERTIN',3*i,'VERTOUT',3*i,'NETFLUX',3*i			!PETER
    do iz = 1, NZ													!PETER
      do ibin = 1, NBIN													!PETER
        write(lundiag,'(2i5,$)') iz, ibin								!PETER
        do ielem = 3*(i-1)+1,3*i-1											!PETER
          write(lundiag,'(3e14.3,$)') vertin_sum(iz,ibin,ielem), vertout_sum(iz,ibin,ielem), &					!PETER
						      vertin_sum(iz,ibin,ielem)-vertout_sum(iz,ibin,ielem)		!PETER
        end do														!PETER
        write(lundiag,'(3e14.3)') vertin_sum(iz,ibin,3*i), vertout_sum(iz,ibin,3*i), &						!PETER
						    vertin_sum(iz,ibin,3*i)-vertout_sum(iz,ibin,3*i)			!PETER
      end do														!PETER
    end do														!PETER
    write(lundiag,*) ' ' 										!PETER
  end do														!PETER
  if (3*elemmultiple+1 .le. NELEM) then											!PETER
    write(lundiag,'(2A5,$)') 'Z', 'BIN'								!PETER
    do ielem = 3*elemmultiple+1,NELEM-1											!PETER
      write(lundiag,'(A12,i2,A12,i2,A12,i2,$)') 'VERTIN',ielem,'VERTOUT',ielem,'NETFLUX',ielem	!PETER
    end do														!PETER
    write(lundiag,'(A12,i2,A12,i2,A12,i2,$)') 'VERTIN',NELEM,'VERTOUT',NELEM,'NETFLUX',NELEM  		!PETER
    do iz = 1, NZ													!PETER
      do ibin = 1, NBIN													!PETER
        write(lundiag,'(2i5,$)') iz, ibin								!PETER
        do ielem = 3*elemmultiple+1,NELEM-1										!PETER
           write(lundiag,'(3e14.3,$)') vertin_sum(iz,ibin,ielem), vertout_sum(iz,ibin,ielem), &					!PETER
						      vertin_sum(iz,ibin,ielem)-vertout_sum(iz,ibin,ielem)		!PETER
        end do														!PETER
        write(lundiag,'(3e14.3)') vertin_sum(iz,ibin,NELEM), vertout_sum(iz,ibin,NELEM), &					!PETER
						    vertin_sum(iz,ibin,NELEM)-vertout_sum(iz,ibin,NELEM)		!PETER
      end do														!PETER
    end do														!PETER
  end if														!PETER

  end if														!PETER
  

  do iz = 1,NZ                                                                 !PETER
!    vtransg(iz) = 8.0e-6_f / (rhoa(iz) / (xmet(iz)*ymet(iz)*zmet(iz)))         !PETER
!    vtransg(iz) = 0._f						               !PETER
!    vtransg(iz) = -8.0e-6_f / (rhoa(iz) / (xmet(iz)*ymet(iz)*zmet(iz)))         !PETER
    vtransg(iz) = winds(iz)							!PETER
  enddo   								       !PETER

  vtransg(NZP1) = vtransg(NZ) 						       !PETER

  nzm1 = max(1, NZ-1)                                                          !PETER

  if (NZ .gt. 1) then                                                          !PETER
    vtransg(NZ) = (vtransg(nzm1) + vtransg(NZ)) / 2._f                         !PETER
 
    if (NZ .gt. 2) then                                                        !PETER
      do iz = NZ-1, 2, -1                                                      !PETER
        vtransg(iz) = (vtransg(iz-1) + vtransg(iz)) / 2._f                     !PETER
      enddo                                                                    !PETER
    endif                                                                      !PETER
  endif                                                                        !PETER

  do igas = 1,NGAS                                                             !PETER
 !       write(*,*) gc(:,igas)
    
    call vertgas(carma, cstate, gc(:,igas), itbnd_gc, ibbnd_gc, &              !PETER
              ftopgas(igas), fbotgas(igas), &                                  !PETER
              gc_topbnd(igas), gc_botbnd(igas), vtransg, rc, igas) !PETER
    if (rc < RC_OK) return                                                     !PETER
  !      write(*,*) 'after', gc(:,igas)
  enddo                                                                        !PETER 


!  gasmultiple = int(NGAS / 3._f)

  if (do_printdiag) then									!PETER
    vertupin_sum(:,:,:) = 0._f		!PETER
    vertupout_sum(:,:,:) = 0._f		!PETER
    vertdnin_sum(:,:,:) = 0._f		!PETER
    vertdnout_sum(:,:,:) = 0._f		!PETER
    aa = 0._f	!PETER
    bb = 0._f	!PETER
    cc = 0._f	!PETER
    dd = 0._f	!PETER
    ee = 0._f	!PETER
    ff = 0._f	!PETER
    vertin_sum = vertupin_sum + vertdnin_sum							!PETER 
    vertout_sum = vertupout_sum + vertdnout_sum							!PETER
    do igas = 1, NGAS										!PETER
      write(lundiag,'(A5,A15,i2,A15,i2,A15,i2,A15,i2,A15,i2,A15,i2,A15,i2)') &			!PETER
	    'Z', 'VERTUPING', igas, 'VERTUPOUTG', igas, 'VERTDNING', igas, &			!PETER
	    'VERTDNOUTG', igas, 'VERTINMASSG', igas, 'VERTOUTMASSG', igas, 'NETFLUXMASSG', igas	!PETER
      do iz = 1, NZ										!PETER
        write(lundiag,'(i5, 7e15.3)') iz, vertupin_sum(iz,igas,1), vertupout_sum(iz,igas,1), &	!PETER
	  			      vertdnin_sum(iz,igas,1), vertdnout_sum(iz,igas,1), &	!PETER
	  			      vertupin_sum(iz,igas,1) + vertdnin_sum(iz,igas,1), &	!PETER
				      vertupout_sum(iz,igas,1) + vertdnout_sum(iz,igas,1), &	!PETER
				      vertupin_sum(iz,igas,1) + vertdnin_sum(iz,igas,1) - &	!PETER
				      vertupout_sum(iz,igas,1) - vertdnout_sum(iz,igas,1)	!PETER
      end do											!PETER
      write(lundiag,*) ' '									!PETER
    end do											!PETER
    write(lundiag,'(A5,7A15)') 'Z', 'VERTUPING', 'VERTUPOUTG', 'VERTDNING', 'VERTDNOUTG', &		!PETER
			       'VERTINMASSG', 'VERTOUTMASSG', 'NETFLUXMASSG'			!PETER
    do iz = 1, NZ										!PETER
      write(lundiag,'(i5, 7e15.3)') iz, sum(vertupin_sum(iz,:,1)), sum(vertupout_sum(iz,:,1)), &!PETER
	  			    sum(vertdnin_sum(iz,:,1)), sum(vertdnout_sum(iz,:,1)), &	!PETER
	  			    sum(vertupin_sum(iz,:,1)) + sum(vertdnin_sum(iz,:,1)), &	!PETER
				    sum(vertupout_sum(iz,:,1)) + sum(vertdnout_sum(iz,:,1)), &	!PETER
				    sum(vertupin_sum(iz,:,1)) + sum(vertdnin_sum(iz,:,1)) - &	!PETER
				    sum(vertupout_sum(iz,:,1)) - sum(vertdnout_sum(iz,:,1))	!PETER
    end do											!PETER
    write(lundiag,*) ' '									!PETER
    do iz = 1, NZ										!PETER
      netflux(iz) = sum(vertupin_sum(iz,:,1)) + sum(vertdnin_sum(iz,:,1)) - &			!PETER
		    sum(vertupout_sum(iz,:,1)) - sum(vertdnout_sum(iz,:,1))			!PETER
    end do											!PETER
    write(lundiag,*) ''													!PETER
    do iz = 1, NZ														!PETER
      aa = aa + sum(vertupin_sum(iz,:,1))					!PETER
      bb = bb + sum(vertupout_sum(iz,:,1))					!PETER
      cc = cc + sum(vertdnin_sum(iz,:,1))					!PETER
      dd = dd + sum(vertdnout_sum(iz,:,1))					!PETER
      ee = ee + sum(vertupin_sum(iz,:,1)) + sum(vertdnin_sum(iz,:,1))	!PETER
      ff = ff + sum(vertupout_sum(iz,:,1)) + sum(vertdnout_sum(iz,:,1))	!PETER
    enddo  														!PETER
    write(lundiag,'(A5,7e15.3)') 'TOT:', aa, bb, cc, dd, ee, ff, sum(netflux)								!PETER
    write(lundiag,*) ''
    write(lundiag,*) 'COLUMN GAS FLUX = TOT * NZ * deltaz (in meters) * (100 cm / m)'			!PETER
    vertgasflux = sum(netflux)*NZ*200._f*100._f
    write(lundiag,'(A28, e22.10)') 'COLUMN GAS FLUX [g/cm2/s]: ', vertgasflux	!PETER
    write(lundiag,*) ' '									!PETER
  end if											!PETER


!  if (do_printdiag) then										!PETER

!  do i = 1, gasmultiple													!PETER
!    write(lundiag,'(A5,$)') 'Z'									!PETER
!    do ielem = 3*(i-1)+1,3*i-1												!PETER
!      write(lundiag,'(A12,i2,A12,i2,A12,i2,$)') 'VERTING',ielem,'VERTOUTG',ielem,'NETFLUXG',ielem	!PETER
!    end do														!PETER
!    write(lundiag,'(A12,i2,A12,i2,A12,i2)') 'VERTING',3*i,'VERTOUTG',3*i,'NETFLUXG',3*i		!PETER
!    do iz = 1, NZ													!PETER
!      write(lundiag,'(i5,$)') iz									!PETER
!      do igas = 3*(i-1)+1,3*i-1												!PETER
!        write(lundiag,'(3e14.3,$)') vertin_sum(iz,igas,1), vertout_sum(iz,igas,1), &						!PETER
!						      vertin_sum(iz,igas,1)-vertout_sum(iz,igas,1)			!PETER
!      end do														!PETER
!      write(lundiag,'(3e14.3)') vertin_sum(iz,3*i,1), vertout_sum(iz,3*i,1), &						!PETER
!						  vertin_sum(iz,3*i,1)-vertout_sum(iz,3*i,1)				!PETER
!    end do														!PETER
!    write(lundiag,*) ' ' 										!PETER
!  end do														!PETER
!  if (3*gasmultiple+1 .le. NGAS) then											!PETER
!    write(lundiag,'(A5,$)') 'Z'								 	!PETER
!    do ielem = 3*gasmultiple+1,NGAS-1											!PETER
!      write(lundiag,'(A12,i2,A12,i2,A12,i2,$)') 'VERTING',ielem,'VERTOUTG',ielem,'NETFLUXG',ielem	!PETER
!    end do														!PETER
!    write(lundiag,'(A12,i2,A12,i2,A12,i2)') &	 		 					!PETER
!  							'VERTING',NGAS,'VERTGOUT',NGAS,'NETFLUXG',NGAS  		!PETER
 !   do iz = 1, NZ													!PETER
!        write(lundiag,'(i5,$)') iz								 	!PETER
 !       do igas = 3*gasmultiple+1,NGAS-1										!PETER
!          write(lundiag,'(3e14.3,$)') vertin_sum(iz,igas,1), vertout_sum(iz,igas,1), &						!PETER
!						      vertin_sum(iz,igas,1)-vertout_sum(iz,igas,1)			!PETER
!        end do														!PETER
!        write(lundiag,'(3e14.3)') vertin_sum(iz,NGAS,1), vertout_sum(iz,NGAS,1), &						!PETER
!						    vertin_sum(iz,NGAS,1)-vertout_sum(iz,NGAS,1)			!PETER
!    end do														!PETER
!  end if														!PETER

!  end if										!PETER

!  if (do_printdiag) write(lundiag,*) ' '					!PETER

  ! Return to caller with new particle concentrations.
  return
end
