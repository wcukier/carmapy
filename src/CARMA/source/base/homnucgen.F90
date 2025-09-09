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
    
    if ((igas .ne. 0) .and. (inucproc(ielem,ielem) .eq. I_HOMGEN)) then

      if (igas .eq. igash2so4) then
        rho_cond = sulfate_density(carma, wtpct(iz), t(iz), rc)
      else ! WC 
        rho_cond = carma%f_gas(igas)%f_rho_cond
      endif

      ! invert formula from setupgkern.f90
      surftens = akelvin(iz,igas) * t(iz) * rho_cond * RGAS / (2._f * gwtmol(igas))
      rvap = RGAS / gwtmol_dif(igas)
      gc_cgs = gc(iz,igas) / (zmet(iz)*xmet(iz)*ymet(iz))


      ! a_c (Gao 2017 Eqn A.2)
      agnuc = max(0._f,akelvin(iz,igas)/log(supsatl(iz,igas) + 1._f))

      if (agnuc .eq. 0._f) then
        nucbin  = 0
        nucrate = 0._f
      else
        ! F (Gao 2017 Eqn A.2)
        deltafg = 4._f / 3._f * PI * surftens * agnuc**2._f 		! For user's choice condensate only  
        
        ! Gao 2017 Eqn A.4 (note typo in paper: p is *partial* pressure)
        fluxmol = gc_cgs * rvap * t(iz) / sqrt(2._f * PI * &
        (gwtmol_dif(igas) / AVG) * BK * t(iz))

        ! g (number of molecules in a particle with radius a_c)
        molgerm = 4._f / 3._f * PI * rho_cond * agnuc**3._f &
        / gwtmol(igas) * AVG
        
        cmass = 4._f / 3._f * PI * agnuc**3._f * rho_cond

        ! Gao 2017 Eqn A.1
        nucrate = 4._f * PI * agnuc**2._f * fluxmol &
	      * sqrt(deltafg / (3._f * PI * BK * t(iz) * molgerm**2._f)) & ! Z (Eqn A.5)
         * (gc_cgs / gwtmol_dif(igas) * AVG) * & ! n
  	    exp( -1._f * deltafg / (BK * t(iz)))

        !   Calc bin # of crit nucleus
        if (cmass .lt. rmassup(1,igroup)) then
          nucbin = 1
        else
          nucbin = 2 + int(log(cmass / rmassup(1,igroup)) / log(rmrat(igroup)))
        endif
  		
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
               
      ! Do further calculations only if nucleation occurred
      if (nucrate .gt. 0._f) then

        rhompe(nucbin, ielem) = rhompe(nucbin, ielem) + nucrate !* redugrow(igas)
        
        ! Since homogeneous nucleation doesn't go through upgxfer or downgxfer, then
        ! then the effects of latent heat need to be accounted for here.
        ! rlprod = rlprod + rhompe(nucbin, ielem) * rmass(nucbin,igroup) * rlh_nuc(ielem,ielem) / (CP * rhoa(iz))
      end if
    end if
  end do

  return
end
