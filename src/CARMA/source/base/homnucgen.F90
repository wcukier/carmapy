! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  Calculates particle production rates due to nucleation <rhompe>:
!!  classic homogeneous nucleation. Numerical method follows Pandis,
!!  2005, Fundamentals of Atmospheric Modeling, 2nd Edition, Cambridge
!!  University Press, pp. 486.
!!
!!  @author Peter Gao
!!  @version Apr-2016
subroutine homnucgen(carma,cstate, iz, rc) 
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
  
  type(carma_type), intent(inout)      :: carma       !! the carma object
  type(carmastate_type), intent(inout) :: cstate      !! the carma state object
  integer, intent(in)                  :: iz          !! level index
  integer, intent(inout)               :: rc          !! return code, negative indicates failure

  !  Local declarations     
  integer           :: igroup     ! group index
  integer           :: ielem      ! concentration element index
  integer           :: igas      ! gas index
  integer           :: nucbin     ! bin in which nucleation takes place    
  real(kind=f)      :: nucrate    ! nucleation rate (#/x/y/z/s)
  real(kind=f)      :: cmass      ! critical cluster mass (g)
  real(kind=f)      :: rho_cond   ! condensate mass density (g/cm3)
  real(kind=f)      :: molgerm      ! Zeldovich factor (for homogeneous nucleation) 
  real(kind=f)      :: gc_cgs
  real(kind=f)      :: surftens
  real(kind=f)      :: rvap
  real(kind=f)      :: agnuc
  real(kind=f)      :: deltafg
  real(kind=f)      :: fluxmol
 
  !write(*,*) igroup, ielem, igas, nucbin, nucrate, cmass, rho_cond, molgerm, gc_cgs, surftens, rvap, agnuc, deltafg, fluxmol    

  ! Cycle through each group, only proceed if BHN
  do igroup = 1 , NGROUP
    
    ielem = ienconc(igroup)      
    igas = inucgas(ielem,ielem)
    ! write(*,*) iz,ielem,igas
    
    if ((igas .ne. 0) .and. (inucproc(ielem,ielem) .eq. I_HOMGEN)) then

      if (igas .eq. igash2so4) then
        rho_cond = sulfate_density(carma, wtpct(iz), t(iz), rc)
      else ! WC 
        rho_cond = carma%f_gas(igas)%f_rho_cond
      endif
      ! else if (igas .eq. igash2o) then
      !   rho_cond = RHO_W
      ! else if (igas .eq. igash2so4) then
      !   rho_cond = sulfate_density(carma, wtpct(iz), t(iz), rc)
      ! else if ((igas .eq. igass8) .or. (igas .eq. igass2)) then
      !   rho_cond = RHO_SX
      ! else if (igas .eq. igaskcl) then
	    !   rho_cond = RHO_KCL
      ! else if (igas .eq. igaszns) then
	    !   rho_cond = RHO_ZNS
      ! else if (igas .eq. igasna2s) then
      !   rho_cond = RHO_NA2S
      ! else if (igas .eq. igasmns) then
      !   rho_cond = RHO_MNS
      ! else if (igas .eq. igascr) then
      !   rho_cond = RHO_CR
      ! else if (igas .eq. igasfe) then
      !   rho_cond = RHO_FE
      ! else if (igas .eq. igasmg2sio4) then
      !   rho_cond = RHO_MG2SIO4
      ! else if (igas .eq. igastio2) then
      !   rho_cond = RHO_TIO2
      ! else if (igas .eq. igasal2o3) then
      !   rho_cond = RHO_AL2O3
      ! else if (igas .eq. igasco) then
      !   rho_cond = RHO_CO
      ! endif

      surftens = akelvin(iz,igas) * t(iz) * rho_cond * RGAS / (2._f * gwtmol(igas))
      rvap = RGAS / gwtmol_dif(igas)
      gc_cgs = gc(iz,igas) / (zmet(iz)*xmet(iz)*ymet(iz))

      ! write(*, *) 'akelvin:', akelvin(iz,igas), 'supsatl:', supsatl(iz,igas), 'iz:', iz, 'igas:', igas


      agnuc = max(0._f,akelvin(iz,igas)/log(supsatl(iz,igas) + 1._f))

      if (agnuc .eq. 0._f) then
        nucbin  = 0
        nucrate = 0._f
      else
        deltafg = 4._f / 3._f * PI * surftens * agnuc**2._f 		! For user's choice condensate only  
	        fluxmol = gc_cgs * rvap * t(iz) / sqrt(2._f * PI * &
	        (gwtmol_dif(igas) / AVG) * BK * t(iz))
        	molgerm = 4._f / 3._f * PI * rho_cond * agnuc**3._f &
	        / gwtmol(igas) * AVG

        
        cmass = 4._f / 3._f * PI * agnuc**3._f * rho_cond

        nucrate = 4._f * PI * agnuc**2._f * fluxmol &
	      * sqrt(deltafg / (3._f * PI * BK * t(iz) * &
        molgerm**2._f)) * (gc_cgs / gwtmol_dif(igas) * AVG) * &
  	    exp( -1._f * deltafg / (BK * t(iz)))
  	  
!        if (igas .eq. 2) then
!        write(*,*) 'iz, igas, dtime, nucrate, surftens, agnuc, fluxmol, &
!	  deltafg, molgerm, gc_cgs,  -1._f * deltafg / (BK * t(iz)),supsatl(iz,igas) + 1._f'
        ! write(*,*) iz, igas, dtime, nucrate, surftens, agnuc, fluxmol, &
        !  deltafg, molgerm, gc_cgs,  -1._f * deltafg / (BK * t(iz)),supsatl(iz,igas) + 1._f
!       endif

       !   Calc bin # of crit nucleus
        if (cmass .lt. rmassup(1,igroup)) then
          nucbin = 1
        else
          nucbin = 2 + int(log(cmass / rmassup(1,igroup)) / log(rmrat(igroup)))
        endif
 
 		    ! write(*,*) 'hello', nucrate, cmass, akelvin(iz,igas), supsatl(iz,igas), log(supsatl(iz,igas) + 1._f), log(1._f)
 		
        ! If none of the bins are large enough for the critical radius, then
        ! no nucleation will occur.
        if ((nucbin > NBIN).OR.(nucbin < 0)) then
          nucbin  = 0 
          nucrate = 0._f
        else
          ! write(*, *) nucbin
  	      nucrate = nucrate*cmass/rmass(nucbin,igroup)
        endif
      endif

      ! Scale to #/x/y/z/s
      nucrate = nucrate * zmet(iz) * xmet(iz) * ymet(iz)

      !if (nucbin .lt. 2 .and. nucbin .gt. 0) then
      !write(*,*) iz, nucrate, nucbin, cmass, agnuc(iz,igas), fluxmol(iz,igas), deltafg(iz, igas), molgerm, 4._f * PI * agnuc(iz,igas)**2._f * fluxmol(iz,igas) &
!	  * sqrt(deltafg(iz,igas) / (3._f * PI * BK * t(iz) * &
 !         molgerm**2._f)) * (gc_cgs / gwtmol(igas) * AVG), exp( -1._f * deltafg(iz,igas) / (BK * t(iz)))
  !    endif
               
      ! Do further calculations only if nucleation occurred
      if (nucrate .gt. 0._f) then

        rhompe(nucbin, ielem) = rhompe(nucbin, ielem) + nucrate !* redugrow(igas)

	!write(*,*) "hom", iz, nucbin, rhompe(nucbin, ielem), agnuc(iz,igas)**2._f, fluxmol(iz,igas), &
	!  sqrt(deltafg(iz,igas) / (3._f * PI * BK * t(iz) * &
        !  molgerm**2._f)), (gc_cgs / gwtmol(igas) * AVG), &
  	!  exp( -1._f * deltafg(iz,igas) / (BK * t(iz)))
        
        ! Since homogeneous nucleation doesn't go through upgxfer or downgxfer, then
        ! then the effects of latent heat need to be accounted for here.
!        rlprod = rlprod + rhompe(nucbin, ielem) * rmass(nucbin,igroup) * rlh_nuc(ielem,ielem) / (CP * rhoa(iz))
      end if
    end if
  end do

  return
end
