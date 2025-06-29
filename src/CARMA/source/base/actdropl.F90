! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle loss rates due to nucleation <rnuclg>:
!! droplet activation only.
!!  
!! The loss rates for all particle elements in a particle group are equal.
!!
!! To avoid nucleation into an evaporating bin, this subroutine must
!! be called after growp, which evaluates evaporation loss rates <evaplg>.
!!
!! @author Andy Ackerman
!! @version Dec-1995
!subroutine actdropl(carma, cstate, iz, rc)
subroutine actdropl(carma, cstate, iz, rc, maxrate)           !PETER

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
  integer, intent(in)                  :: iz      !! z index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure
  real(kind=f), optional, intent(in)   :: maxrate !! PETER

  !  Local declarations
  integer                              :: igas    !! gas index
  integer                              :: igroup  !! group index
  integer                              :: ibin    !! bin index
  integer                              :: iepart  !! element for condensing group index
  integer                              :: inuc    !! nucleating element index
  integer                              :: ienucto !! index of target nucleation element
  integer                              :: ignucto !! index of target nucleation group
  integer                              :: inucto  !! index of target nucleation bin
  integer                              :: ccn_gas  !! gas that makes up seed
  logical                              :: evapfrom_nucto !! .true. when target droplets are evaporating
  real(kind=f)      :: rho_cond   ! condensate mass density (g/cm3)
  real(kind=f)      :: molgerm      ! Zeldovich factor (for homogeneous nucleation) 
  real(kind=f)      :: gc_cgs
  real(kind=f)      :: surftens
  real(kind=f)      :: rvap
  real(kind=f)      :: agnuc
  real(kind=f)      :: deltafg
  real(kind=f)      :: surfcond
  real(kind=f)      :: fluxmol
  real(kind=f)      :: curvfact
  real(kind=f)      :: zv
  real(kind=f)      :: radratio
  real(kind=f)                   :: fo
  real(kind=f)                   :: phi
  real(kind=f)                   :: temp1
  real(kind=f)                   :: imucos
  real(kind=f)                   :: Wsl !adhesion, used to calculate contact angle! use early Berthelot relation (approximate as entirely polar or nonpolar)
  real(kind=f)                   :: interfacial !interfacial tension
  !real(kind=f), parameter        :: osfreq = 1.0e13_f
  !real(kind=f), parameter        :: dsrpnrgy = 0.18_f ! for water only! Units are eV
  !real(kind=f), parameter        :: osfreq = 3.0e12_f
  !real(kind=f), parameter        :: dsrpnrgy = 2._f ! for KCl on Ag! Units are eV
  real(kind=f)        :: osfreq

  ! This calculation is only necessary for temperatures greater
  ! than -40C.

  ! PETER: By passing temperature restrictions for H2SO4 calculations
!  if( t(iz) .ge. (T0 - 40._f) ) then
  if( t(iz) .ge. 0._f ) then

    ! Loop over particle groups.
    do igroup = 1,NGROUP
    
      ! Bypass calculation if few particles are present
      if( pconmax(iz,igroup) .gt. FEW_PC )then
  
        iepart = ienconc( igroup )            ! particle number density element
  
        ! Calculate nucleation loss rates.  Do not allow nucleation into
        ! an evaporating bin.
        do inuc = 1,nnuc2elem(iepart)

          igas = inucgas(iepart,inuc2elem(inuc,iepart))                ! condensing gas
	 ! write(*,*) iz, inuc, iepart, inuc2elem(inuc,iepart), igas

	  if( igas .ne. 0 )then
          
            ienucto = inuc2elem(inuc,iepart)
            if( ienucto .ne. 0 )then
              ignucto = igelem( ienucto )
            else
              ignucto = 0
            endif

	  !  write(*,*) iz,igas,igroup,iepart,inuc,ienucto,ignucto
      
	!write(*,*) iepart,ienucto,	inucproc(iepart,ienucto)

            ! Only compute nucleation rate for droplet activation
            if( (inucproc(iepart,ienucto) .eq. I_DROPACT) .or.  &
	      (inucproc(iepart,ienucto) .eq. I_HETGEN) ) then

 	!write(*,*) "inucproc"
      
              ! Loop over particle bins.  Loop from largest to smallest for 
              ! evaluation of index of smallest bin nucleated during time step <inucstep>.
              do ibin = NBIN, 1, -1
      
                if( ignucto .ne. 0 )then
                  inucto = inuc2bin(ibin,igroup,ignucto)
                else
                  inucto = 0
                endif
      
                ! Set <evapfrom_nucto> to .true. when target droplets are evaporating
                if( inucto .ne. 0 )then
                 evapfrom_nucto = evaplg(inucto,ignucto) .gt. 0._f
                else
                 evapfrom_nucto = .false.
                endif
                  
!                if( (supsatl(iz,igas) .gt. scrit(iz,ibin,igroup)) .and. &
!                    (.not. evapfrom_nucto) .and. &
!                    (pc(iz,ibin,iepart) .gt. SMALL_PC) )then
                if( (.not. evapfrom_nucto) .and. &
                    (pc(iz,ibin,iepart) .gt. SMALL_PC) )then

                  if( inucproc(iepart,ienucto) .eq. I_DROPACT ) then
		    if (supsatl(iz,igas) .gt. scrit(iz,ibin,igroup,igas)) then
 		      rnuclg(ibin,igroup,ignucto) = 1.0e3_f
		    endif
		  else  
!TODO - WC
                    if (igas .eq. igash2o) then
                      rho_cond = RHO_W
                    else if (igas .eq. igash2so4) then
                      rho_cond = sulfate_density(carma, wtpct(iz), t(iz), rc)
                    else if ((igas .eq. igass8) .or. (igas .eq. igass2)) then
                      rho_cond = RHO_SX
                    else if (igas .eq. igaskcl) then
	              rho_cond = RHO_KCL
                    else if (igas .eq. igaszns) then
  	              rho_cond = RHO_ZNS
                    else if (igas .eq. igasna2s) then
  	              rho_cond = RHO_NA2S
                    else if (igas .eq. igasmns) then
  	              rho_cond = RHO_MNS
                    else if (igas .eq. igascr) then
  	              rho_cond = RHO_CR
                    else if (igas .eq. igasfe) then
  	              rho_cond = RHO_FE
                    else if (igas .eq. igasmg2sio4) then
  	              rho_cond = RHO_MG2SIO4
                    else if (igas .eq. igastio2) then
  	              rho_cond = RHO_TIO2
                    else if (igas .eq. igasal2o3) then
  	              rho_cond = RHO_AL2O3
                    else if (igas .eq. igasco) then
  	              rho_cond = RHO_CO
                    endif
	  
                    surftens = akelvin(iz,igas) * t(iz) * rho_cond * RGAS / (2._f * gwtmol(igas))
	            rvap = RGAS / gwtmol_dif(igas)
                    gc_cgs = gc(iz,igas) / (zmet(iz)*xmet(iz)*ymet(iz))

	            agnuc = max(0._f,akelvin(iz,igas)/log(supsatl(iz,igas) + 1._f))
	            
	            ccn_gas = igrowgas(iepart)

  	            if (agnuc .eq. 0._f) then
                      rnuclg(ibin,igroup,ignucto) = 0._f
                    else
                    ! DK POWELL, check if there is an updated seed composition
                     !if (mucos(igas,igroup) .eq. 0._f) then
                     !
                     !  Wsl = 2._f * sqrt(surfacetens(iz,ccn_gas)*surftens)
	                 !  imucos =  Wsl/surftens - 1._f
	                 !  
	                 !  if (imucos .gt. 0.995_f) then
	                 !     imucos = 0.995_f
	                 !  endif
	                 !  
	                 !  if (imucos .lt. -1._f) then
	                 !     imucos = -1._f
	                 !  endif
	                 !  
	                 !else
	                 !  imucos = mucos(igas,igroup)
	                 !endif
	                 
	                if (mucos(igas,igroup) .eq. 0._f) then
                     ! Set minimum contact angle to 0.1 degrees
                     imucos = min(surfacetens(iz,igrowgas(iepart))/surftens,0.9999984769_f)
                     !print*,igrowgas(iepart),surf_tens(iz,igrowgas(iepart)),surftens,agnuc
                   else
                     imucos = mucos(igas,igroup)
                   endif

                      osfreq = 1.6e11_f * sqrt(desorption(igas)*EV2ERG/BK/gwtmol_dif(igas))
                      deltafg = 4._f / 3._f * PI * surftens * agnuc**2._f 		! For user's choice condensate only  
	              fluxmol = gc_cgs * rvap * t(iz) / sqrt(2._f * PI * &
	                (gwtmol_dif(igas) / AVG) * BK * t(iz))
	              !surfcond = fluxmol / osfreq * exp(dsrpnrgy * &
	              surfcond = fluxmol / osfreq * exp(desorption(igas) * &
	                EV2ERG / (BK * t(iz)))
	              molgerm = 4._f / 3._f * PI * rho_cond * agnuc**3._f &
	                / gwtmol(igas) * AVG
	              radratio = r(ibin,igroup)/agnuc
  	              phi = sqrt(1._f - 2._f * imucos * radratio + radratio**2._f)
	              fo = (radratio - imucos) / phi
	              temp1 = 2._f - 3._f*fo + fo**3._f
	              if (log10(abs(temp1)) < -15._f) temp1 = 0._f
	              curvfact = 0.5_f *(1._f + ((1._f - &
	 	        imucos*radratio)/phi)**3._f + radratio**3._f * &
		        temp1 + 3._f * imucos * &
	                radratio**2._f * (fo - 1._f))
                      zv = sqrt(deltafg &
	  	        / (3._f * PI * BK * t(iz) * molgerm**2._f)) * &
	                sqrt(4._f/(2._f + ((1._f - imucos*radratio) * &
	 	        (2._f - 4._f*imucos*radratio - (imucos**2._f - &
		        3._f)*radratio**2._f)) / phi**3._f)) ! H. Vehkamaki, et al. 2007, Technical Note: The heterogeneous Zeldovich factor, Atmos. Chem. Phys., 7, 309-313, 2007		    
                      rnuclg(ibin,igroup,ignucto) = (4._f * PI**2._f * r(ibin,igroup)**2._f * &
		        agnuc**2._f * fluxmol * surfcond * &
		        zv * exp( -1._f * deltafg * &
	                curvfact / (BK * t(iz))))! * redugrow(igas)

				!DPOW Check --> change units of rnuclg to CARMA units if not in I_CART, possibly need to divide instead of multiply
                
                rnuclg(ibin,igroup,ignucto) = rnuclg(ibin,igroup,ignucto)
                        !write(*,*) "iz, igas, igroup, ibin, dtime, ignucto, rnuclg, agnuc, fluxmol, surfcond, zv, deltafg, molgerm, curvfact, imucos, radratio, phi, fo, temp1"
                        !write(*,*) iz, igas, igroup, ibin, dtime, ignucto, rnuclg(ibin,igroup,ignucto), agnuc, fluxmol, surfcond, zv, deltafg, molgerm, curvfact, imucos, radratio, phi, fo, temp1


	             endif

		  endif
                endif
              enddo   ! ibin = 1,NBIN
            endif    ! inucproc(iepart,ienucto) .eq. I_DROPACT
          endif     ! inuc = 1,nnuc2elem(iepart)
        enddo      ! (igas = inucgas(igroup)) .ne. 0 
      endif       ! pconmax(iz,igroup) .gt. FEW_PC
    enddo        ! igroup = 1,NGROUP
  endif         ! t(iz) .ge. T0-40.

  ! Return to caller with particle loss rates due to nucleation evaluated.
  return
end
