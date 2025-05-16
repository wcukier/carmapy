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
!  logical			 :: is_hetgen
  real(kind=f)                   :: bsol
!  real(kind=f)                   :: rvap
!  real(kind=f)                   :: gc_cgs
!  real(kind=f)                   :: rho_cond
!  real(kind=f)                   :: surftens      
!  real(kind=f)                   :: radratio
!  real(kind=f)                   :: fo
!  real(kind=f)                   :: phi
!  real(kind=f)                   :: molgerm
!  real(kind=f)                   :: temp1
!  real(kind=f)                   :: imucos
!  real(kind=f), parameter        :: osfreq = 1.0e13_f
!  real(kind=f), parameter        :: dsrpnrgy = 0.18_f ! for water only! Units are eV
!  real(kind=f), parameter        :: mucos = 0.95_f ! guess ~20 degrees for now. 
!  real(kind=f), parameter        :: mucos = 0._f ! guess ~90 degrees for now. 
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
     ! write(*,*) inuc, ielemfrom, inuc2elem(inuc, ielemfrom), igas

      if( igas .ne. 0 .and. (itype( ielemfrom ) .eq. &
        I_INVOLATILE .or. itype( ielemfrom ) .eq. I_VOLATILE) ) then

        isol = isolelem( ienconc( igroup ) )
      
       ! is_hetgen = .FALSE.
       ! do ielem = 1,NELEM
       !   if (inucproc(ielemfrom,ielem) .eq. I_HETGEN) is_hetgen = .TRUE.
       ! enddo
      
        ! If here is no solute are specified, then no scrit value is defined.
      
 !       do ibin = 1,NBIN
        do k = 1,NZ
        
          if (isol .ne. 0) then
            ! This is term "B" in Pruppacher and Klett's eqn. 6-28.=
	    do ibin = 1,NBIN
              bsol = 3._f*sol_ions(isol)*rmass(ibin,igroup)*gwtmol(igas) &
                / ( 4._f*PI*solwtmol(isol)*RHO_W )
              scrit(k,ibin,igroup,igas) = sqrt( 4._f * akelvin(k,igas)**3 / ( 27._f * bsol ) )
	     ! write(*,*) k, ibin, igroup, bsol, scrit(k,ibin,igroup)
	    enddo
	  endif

!          if (igas .eq. igash2o) then
!            rho_cond = RHO_W
!          else if (igas .eq. igash2so4) then
!            rho_cond = sulfate_density(carma, wtpct(k), t(k), rc)
!          else if ((igas .eq. igass8) .or. (igas .eq. igass2)) then
!            rho_cond = RHO_SX
!          else if (igas .eq. igaskcl) then
!	    rho_cond = RHO_KCL
!          else if (igas .eq. igaszns) then
!  	    rho_cond = RHO_ZNS
!          endif
  
!          surftens = akelvin(k,igas) * t(k) * rho_cond * RGAS / (2._f * gwtmol(igas))
!	  rvap = RGAS / gwtmol(igas)
!          gc_cgs = gc(k,igas) / (zmet(k)*xmet(k)*ymet(k))

!	  agnuc(k,igas) = max(0._f,akelvin(k,igas)/log(supsatl(k,igas) + 1._f))
!          deltafg(k,igas) = 4._f / 3._f * PI * surftens * agnuc(k,igas)**2._f 		! For user's choice condensate only  
!	  fluxmol(k,igas) = gc_cgs * rvap * t(k) / sqrt(2._f * PI * &
!	    (gwtmol(igas) / AVG) * BK * t(k))
!	  surfcond(k,igas) = fluxmol(k,igas) / osfreq * exp(dsrpnrgy * &
!	    1.6e-12_f / (BK * t(k)))
!	  molgerm = 4._f / 3._f * PI * rho_cond * agnuc(k,igas)**3._f &
!	    / gwtmol(igas) * AVG

	!  write(*,*) k, igas, agnuc(k,igas),deltafg(k,igas),fluxmol(k,igas), &
	!	surfcond(k,igas), molgerm
	  
          ! Loop over vertical grid layers because of temperature dependence
          ! in solute term.
  !        do k = 1,NZ
!	  if (is_hetgen) then
!            imucos = mucos(igas,igroup)	
!	    !write(*,*) igas, igroup, mucos(igas,igroup)	
!            do ibin = 1,NBIN
!  	      if (agnuc(k,igas) .eq. 0._f) then
!	        curvfact(k,ibin,igroup,igas) = 0._f
!	        zv(k,ibin,igroup,igas) = 0._f
! 	      else
!	        radratio = r(ibin,igroup)/agnuc(k,igas)
!  	        phi = sqrt(1._f - 2._f * imucos * radratio + radratio**2._f)
!	        fo = (radratio - imucos) / phi
!	        temp1 = 2._f - 3._f*fo + fo**3._f
!	        if (log10(abs(temp1)) < -15._f) temp1 = 0._f
!	        curvfact(k,ibin,igroup,igas) = 0.5_f *(1._f + ((1._f - &
!	 	    imucos*radratio)/phi)**3._f + radratio**3._f * &
!		    temp1 + 3._f * imucos * &
!	            radratio**2._f * (fo - 1._f))
!	        zv(k,ibin,igroup, igas) = sqrt(deltafg(k,igas) &
!	  	    / (3._f * PI * BK * t(k) * molgerm**2._f)) * &
!	            sqrt(4._f/(2._f + ((1._f - imucos*radratio) * &
!	 	    (2._f - 4._f*imucos*radratio - (imucos**2._f - &
!		    3._f)*radratio**2._f)) / phi**3._f)) ! H. Vehkamaki, et al. 2007, Technical Note: The heterogeneous Zeldovich factor, Atmos. Chem. Phys., 7, 309-313, 2007		    
!	      endif
!	  write(LUNOPRT,*) k, igas, ibin, igroup, radratio, phi, fo, &
!		0.5_f *(1._f + ((1._f - &
!	    	  imucos*radratio)/phi)**3._f + radratio**3._f * &
!		  (2._f - 3._f*fo + fo**3._f) + 3._f * imucos * &
!	          radratio**2._f * (fo - 1._f)), zv(k,ibin,igroup,igas), imucos

!	  write(LUNOPRT,*) k, igas, ibin, igroup, temp1, &
!		2._f, -3._f*fo,fo**3._f, zv(k,ibin,igroup,igas)


!             scrit(k,ibin,igroup) = exp(akelvin(k,igas)/r(ibin,igroup)) - 1._f !!DP Change for getting rid of solute effect
!		write(LUNOPRT,*) akelvin(k,igas), r(ibin,igroup), k, ibin, igroup
             
!            enddo
!	  endif
        enddo
      endif
    enddo
  enddo

!#ifdef DEBUG
!  if (do_print_init) then
!    do isol = 1,NSOLUTE
  
!      write(LUNOPRT,3) 'solute name:    ',solname(isol)
  
!      do igroup = 1,NGROUP
!        if( isol .eq. isolelem(ienconc(igroup)) )then
!          write(LUNOPRT,8) groupname(igroup)
!          write(LUNOPRT,6) (i,r(i,igroup),scrit(1,i,igroup),i=1,NBIN)
!        endif
!      enddo
!	  enddo
!	endif
!#endif

  ! Return to caller with nucleation mapping arrays and critical
  ! supersaturations defined.
  return
end
