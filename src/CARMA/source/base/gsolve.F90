! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine calculates new gas concentrations.
!!
!! @author Andy Ackerman, Bill McKie, Chuck Bardeen
!! @version Dec-1995, Sep-1997, Nov-2009
subroutine gsolve(carma, cstate, iz, previous_ice, previous_liquid, rc)

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
  real(kind=f), intent(in)             :: previous_ice(NGAS)      !! total ice at the start of substep
  real(kind=f), intent(in)             :: previous_liquid(NGAS)   !! total liquid at the start of substep
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local Variables
  integer                              :: igas    !! gas index
  real(kind=f)                         :: gc_cgs
  real(kind=f)                         :: gc_old
  real(kind=f)                         :: rvap
  real(kind=f)                         :: stofact
  real(kind=f)                         :: total_ice(NGAS)      ! total ice
  real(kind=f)                         :: total_liquid(NGAS)   ! total liquid
  real(kind=f)                         :: gc_threshold         ! threshold for changes to gc
  
  
  1 format(/,'gsolve::ERROR - negtive gas concentration for ',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2,',gc=',e10.3,',gasprod=',e10.3,',supsati=',e10.3, &
              ',supsatl=',e10.3,',t=',f6.2)
  2 format('gsolve::ERROR - conditions at beginning of the step : gc=',e10.3,',supsati=',e17.10, &
              ',supsatl=',e17.10,',t=',f6.2,',d_gc=',e10.3,',d_t=',f6.2)
  3 format(/,'microfast::WARNING - gas concentration change exceeds threshold: ',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2, ', (gc-gcl)/gcl=', e10.3)
  4 format(/,'gsolve::WARNING: negtive gas concentration reached for ',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2,',gc=',e10.3,',gasprod=',e10.3,',supsati=',e10.3, &
              ',supsatl=',e10.3,',t=',f6.2, ', setting gc to 1e-50*rhoa(iz)')
  

  ! Determine the total amount of condensate for each gas.
  call totalcondensate(carma, cstate, iz, total_ice, total_liquid, rc)
  
  do igas = 1,NGAS
    stofact = carma%f_gas(igas)%f_wtmol_dif/carma%f_gas(igas)%f_wtmol !WC
  !  if (igas .eq. igash2o) then
  !    stofact = 1._f
  !  else if (igas .eq. igash2so4) then
  !    stofact = 1._f
  !  else if ((igas .eq. igass8) .or. (igas .eq. igass2)) then
  !    stofact = 1._f
  !  else if (igas .eq. igaskcl) then
  !    stofact = 1._f
  !  else if (igas .eq. igaszns) then
  !    stofact = WTMOL_ZN/WTMOL_ZNS
    if (igas .eq. igasna2s) then
     stofact = WTMOL_NA2/WTMOL_NA2S
  !  else if (igas .eq. igasmns) then
  !    stofact = WTMOL_MN/WTMOL_MNS
  !  else if (igas .eq. igascr) then
  !    stofact = 1._f
  !  else if (igas .eq. igasfe) then
  !    stofact = 1._f
    else if (igas .eq. igasmg2sio4) then !TODO: Comment out
     stofact = WTMOL_MG2/WTMOL_MG2SIO4
    ! endif
  !  else if (igas .eq. igastio2) then
  !    stofact = 1._f
   else if (igas .eq. igasal2o3) then ! TODO check this
     stofact = WTMOL_AL2/WTMOL_AL2O3
  !  else if (igas .eq. igasco) then
  !    stofact = 1._f
   endif


    ! We do not seem to be conserving mass and energy, so rather than relying upon gasprod
    ! and rlheat, recalculate the total change in condensate to determine the change
    ! in gas and energy.
    !
    ! This is because in the old scheme, the particles were solved for implicitly, but the
    ! gas and latent heat were solved for explicitly using the same rates.
    gasprod(igas) = ((previous_ice(igas) - total_ice(igas)) + &
		(previous_liquid(igas) - total_liquid(igas))) / dtime
    rlprod        = rlprod - ((previous_ice(igas) - total_ice(igas)) * &
		(rlhe(iz,igas) + rlhm(iz,igas)) + &
                       (previous_liquid(igas) - total_liquid(igas)) * &
		(rlhe(iz,igas))) / (CP * rhoa(iz) * dtime) 

    gc_old = gc(iz,igas)

    ! Don't let the gas concentration go negative.
    gc(iz,igas) = gc(iz,igas) + dtime * (gasprod(igas)*stofact + phochemprod_gas(iz,igas))
    !if (igas .eq. 2) then

    !write(*,*) 'iz, igas, gc_old, gc(iz,igas), phochemprod_gas(iz,igas)    '
    !write(*,*) iz, igas, gc_old, gc(iz,igas), phochemprod_gas(iz,igas)    
    !if (igas .eq. 2) then
    !  write(*,*) iz, igas, gc(iz,igas)
    !endif

   ! if (igas .eq. igaszns) then
!	write(*,*) 'gsolve', iz, igas, gc_old, dtime * gasprod(igas), gc(iz,igas)
!		write(*,*) ' '
!	endif

!    if (gc(iz,igas) < 0.0_f) then
      !redugrow(igas) = abs(gc_old / 2._f / (dtime * gasprod(igas)))
      !write(*,*) 'NOPE', iz, igas, redugrow(igas), gc_old, dtime * gasprod(igas)
!      if (do_substep) then
!        if (nretries == maxretries) then 
!          if (do_print) write(LUNOPRT,1) trim(gasname(igas)), iz, &
!                lat, lon, gc(iz,igas), gasprod(igas), &
!                supsati(iz,igas), supsatl(iz,igas), t(iz)
!          if (do_print) write(LUNOPRT,2) gcl(iz,igas), supsatiold(iz,igas), &
!                supsatlold(iz,igas), told(iz), d_gc(iz,igas), d_t(iz)
!        end if
!      else
!        if (do_print) write(LUNOPRT,1) trim(gasname(igas)), iz, &
!                lat, lon, gc(iz,igas), gasprod(igas), &
!                supsati(iz, igas), supsatl(iz,igas), t(iz)
!      end if
!
!      rc = RC_WARNING_RETRY
!    end if

    if (gc(iz,igas) < 0.0_f) then
      if (do_substep) then
        if (nretries == maxretries) then 
          !if (do_print) write(LUNOPRT,1) trim(gasname(igas)), iz, &
          !      lat, lon, gc(iz,igas), gasprod(igas), &
          !      supsati(iz,igas), supsatl(iz,igas), t(iz)
          if (do_print) write(LUNOPRT,4) trim(gasname(igas)), iz, &
                lat, lon, gc(iz,igas), gasprod(igas), &
                supsati(iz,igas), supsatl(iz,igas), t(iz)
          if (do_print) write(LUNOPRT,2) gcl(iz,igas), supsatiold(iz,igas), &
                supsatlold(iz,igas), told(iz), d_gc(iz,igas), d_t(iz)
          gc(iz,igas) = 1e-50_f * rhoa(iz)
        else
          rc = RC_WARNING_RETRY
        endif
      else
        if (do_print) write(LUNOPRT,1) trim(gasname(igas)), iz, &
                lat, lon, gc(iz,igas), gasprod(igas), &
                supsati(iz, igas), supsatl(iz,igas), t(iz)
        rc = RC_WARNING_RETRY
      end if
    end if
    
    ! If gas changes by too much, then retry the calculation.
    gc_threshold = dgc_threshold(igas)
    if (do_incloud) gc_threshold = gc_threshold / cldfrc(iz)
    !
    ! NOTE: If doing incloud calculations, then the threshold needs to be scaled by the
    ! cloud fraction since the cloud mass has been scaled by the cloud fraction.
    
    if (gc_threshold /= 0._f) then
      if ((dtime * gasprod(igas) / gc(iz,igas)) > gc_threshold) then
        if (do_substep) then
          if (nretries == maxretries) then 
            if (do_print) write(LUNOPRT,3) trim(gasname(igas)), iz, &
		lat, lon, dtime * gasprod(igas) / gc(iz,igas)
            if (do_print) write(LUNOPRT,2) gcl(iz,igas), supsatiold(iz,igas), &
		supsatlold(iz,igas), told(iz), d_gc(iz,igas), d_t(iz)
          end if
        else
          if (do_print) write(LUNOPRT,3) trim(gasname(igas)), iz, &
		lat, lon, dtime * gasprod(igas) / gc(iz,igas)
        end if
  
        rc = RC_WARNING_RETRY
      end if
    end if
  end do

  ! Return to caller with new gas concentrations.
  return
end
