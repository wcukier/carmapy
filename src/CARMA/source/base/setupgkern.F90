! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine defines radius-dependent but time-independent parameters
!! used to calculate condensational growth of particles.  Growth rates
!! are calculated at bin boundaries: the parameters calculated here 
!! ( <gro>, <gro1>, <gro2>, and <akelvin> ) 
!! are defined at lower bin boundaries through the growth rate expression
!! (for one particle) used in growevapl.f:
!!>
!!   dm = gro*pvap*( S + 1 - Ak*As - gro1*gro2*qrad )
!!   --   -------------------------------------------
!!   dt               1 + gro*gro1*pvap
!!
!!  where 
!!
!!  S    = supersaturation
!!  Ak   = exp(akelvin/r)
!!  As   = exp(-sol_ions * solute_mass/solwtmol * gwtmol/condensate_mass)
!!  pvap = saturation vapor pressure [dyne cm**-2]
!!  qrad = radiative energy absorbed
!!<
!! This routine requires that vertical profiles of temperature <T>,
!! and pressure <p> are defined.
!!
!! This routine also requires that particle Reynolds' numbers are
!! defined (setupvfall.f must be called before this).
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine setupgkern(carma, cstate, rc)

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
  integer                        :: igas     !! gas index
  integer                        :: ielem    !! element index
  integer                        :: inuc     ! nucleation index
  integer                        :: k        !! z index
  integer                        :: igroup   !! group index
  integer                        :: i
  real(kind=f)                   :: gstick
  real(kind=f)                   :: cor
  real(kind=f)                   :: phish
  real(kind=f)                   :: esh1
  real(kind=f)                   :: a1
  real(kind=f)                   :: br
  real(kind=f)                   :: rknudn
  real(kind=f)                   :: rknudnt
  real(kind=f)                   :: rlam
  real(kind=f)                   :: rlamt
  real(kind=f)                   :: rhoa_cgs(NZ, NGAS)
  real(kind=f)                   :: freep(NZ, NGAS)
  real(kind=f)                   :: freept(NZ, NGAS)
  real(kind=f)                   :: rlh
  real(kind=f)                   :: diffus1
  real(kind=f)                   :: thcond1
  real(kind=f)                   :: reyn_shape
  real(kind=f)                   :: schn
  real(kind=f)                   :: prnum
  real(kind=f)                   :: x1
  real(kind=f)                   :: x2
  real(kind=f)                   :: fv
  real(kind=f)                   :: surf_tens(NZ)  ! surface tension of non-H2O particles
  real(kind=f)                   :: rho_part  ! density of non-H2O particle
  

  ! Calculate gas properties for all of the gases. Better to do them all once, than to
  ! repeat this for multiple groups.
  do igas = 1, NGAS
 
    ! Radius-independent parameters for condensing gas
    !
    ! This is <rhoa> in cgs units.
    !
    rhoa_cgs(:, igas) = rhoa(:) / (xmet(:)*ymet(:)*zmet(:))

    if (igas .eq. igash2o) then
        
      ! Condensing gas is water vapor
      !
      ! <surfctwa> is surface tension of water-air interface (valid from 0 to 40 C)
      ! from Pruppacher and Klett (eq. 5-12).
      surfctwa(:) = 76.10_f - 0.155_f*( t(:) - 273.16_f )

      !PETER: This is the surfctwa for liquid H2SO4 vs. air, from Myhre et al., 1998 (data).
 !     surfctwa(:) = 70.03_f - 0.039_f*( t(:) - T0 )

      ! <surfctiw> is surface tension of water-ice interface
      ! from Pruppacher and Klett (eq. 5-48).!
      surfctiw(:) = 28.5_f + 0.25_f*( t(:) - 273.16_f )

      ! <surfctiw> is surface tension of water-ice interface
      ! from Hale and Plummer [J. Chem. Phys., 61, 1974].
      surfctia(:) = 141._f - 0.15_f * t(:)
	  
	  surfacetens(:,igas) = surfctwa(:) 
	  
      ! <akelvin> is argument of exponential in kelvin curvature term.
      akelvin(:,igas) = 2._f*gwtmol(igas)*surfctwa(:) &
                        / ( t(:)*RHO_W*RGAS )

      akelvini(:,igas) = 2._f*gwtmol(igas)*surfctia(:) & 
                        / ( t(:)*RHO_W*RGAS )
                        
    ! condensing gas is H2SO4                    
    else if (igas .eq. igash2so4) then
      ! Calculate Kelvin curvature factor for H2SO4 interactively with temperature:
      do k = 1, NZ  
        surf_tens(k) = sulfate_surf_tens(carma, wtpct(k), t(k), rc)
        rho_part = sulfate_density(carma, wtpct(k), t(k), rc)
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens(k) / (t(k) * rho_part * RGAS)
        surfacetens(k, igas) = surf_tens(k)
        
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvini(k, igas)
      end do   
    else if ((igas .eq. igass8) .or. (igas .eq. igass2)) then
      ! Calculate Kelvin curvature factor for polysulfur interactively with temperature:
      do k = 1, NZ  
        surf_tens(k) = 60.8_f ! http://cameochemicals.noaa.gov/chris/SXX.pdf, Fanelli 1950
        rho_part = 1.96_f ! http://periodictable.com/Elements/016/data.html
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens(k) / (t(k) * RHO_SX * RGAS)
        surfacetens(k, igas) = surf_tens(k)
        
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvin(k, igas)
      end do 
    else if (igas .eq. igaskcl) then
      ! Calculate Kelvin curvature factor for polysulfur interactively with temperature:
      do k = 1, NZ  
        surf_tens(k) = 160.4_f - 0.07_f*(t(k)-273.15_f) ! http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_5.html 
!        surf_tens(k) = 4._f*(160.4_f - 0.07_f*t(k)) ! http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_5.html 
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens(k) / (t(k) * RHO_KCL * RGAS)
        surfacetens(k, igas) = surf_tens(k)
        
        desorption(igas) = 0.5_f
        
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvin(k, igas)
      end do   
    else if (igas .eq. igaszns) then
      ! Calculate Kelvin curvature factor for polysulfur interactively with temperature:
      do k = 1, NZ  
        !surf_tens(k) = 1672._f ! Celikkaya & Akinc 1990, Journal of the American Ceramic Society 73, 2360
        surf_tens(k) = 860._f ! Zhang et al. 2003, J. Phys. Chem. B 2003, 107, 13051-13060
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens(k) / (t(k) * RHO_ZNS * RGAS)
        surfacetens(k, igas) = surf_tens(k)
        
        desorption(igas) = 0.5_f
        
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvin(k, igas)
      end do  
    else if (igas .eq. igasna2s) then
      ! Calculate Kelvin curvature factor for polysulfur interactively with temperature:
      do k = 1, NZ  
        !surf_tens(k) = 160.4_f - 0.07_f*t(k) ! Assume same as KCl 
        surf_tens(k) = 1033._f ! Graham Lee's calculations, estimated from enthalpy data  
        !surf_tens(k) = 50._f ! DON'T KNOW
        surfacetens(k, igas) = surf_tens(k)
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens(k) / (t(k) * RHO_NA2S * RGAS)
        
        desorption(igas) = 0.5_f
        
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvin(k, igas)
      end do  
    else if (igas .eq. igasmns) then
      ! Calculate Kelvin curvature factor for polysulfur interactively with temperature:
      do k = 1, NZ  
        !surf_tens(k) = 50._f ! DON'T KNOW
        !surf_tens(k) = 160.4_f - 0.07_f*t(k) ! Assume same as KCl  
        surf_tens(k) = 2326._f ! Graham Lee's calculations, estimated from enthalpy data  
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens(k) / (t(k) * RHO_MNS * RGAS)
        surfacetens(k, igas) = surf_tens(k)
        
        desorption(igas) = 0.5_f
        
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvin(k, igas)
      end do  
    else if (igas .eq. igascr) then
      ! Calculate Kelvin curvature factor for polysulfur interactively with temperature:
      do k = 1, NZ  
        surf_tens(k) = 1642._f - 0.2_f * ( t(k) - 2133.15_f ) ! http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_5.html 
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens(k) / (t(k) * RHO_CR * RGAS)
        surfacetens(k, igas) = surf_tens(k)
        
        desorption(igas) = 0.5_f
        
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvin(k, igas)
      end do  
    else if (igas .eq. igasfe) then
      ! Calculate Kelvin curvature factor for polysulfur interactively with temperature:
      do k = 1, NZ  
        surf_tens(k) = 1862._f - 0.39_f * ( t(k) - 1803.15_f ) ! http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_5.html 
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens(k) / (t(k) * RHO_FE * RGAS)
        surfacetens(k, igas) = surf_tens(k)
        
        desorption(igas) = 0.5_f
        
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvin(k, igas)
      end do  
    else if (igas .eq. igasmg2sio4) then
      ! Calculate Kelvin curvature factor for polysulfur interactively with temperature:
      do k = 1, NZ  
        surf_tens(k) = 436._f !Kazasa et al. 1989!1280._f ! Miura et al. (2010)
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens(k) / (t(k) * RHO_MG2SIO4 * RGAS)
        surfacetens(k, igas) = surf_tens(k)
        
        desorption(igas) = 0.5_f
             
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvin(k, igas)
      end do  
    else if (igas .eq. igastio2) then
      ! Calculate Kelvin curvature factor for polysulfur interactively with temperature:
      do k = 1, NZ  
        !surf_tens(k) = 480._f ! Lee et al. (2014)
        surf_tens(k) = 535.124_f - 0.04396_f * t(k) ! Lee et al. (2015)
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens(k) / (t(k) * RHO_TIO2 * RGAS)
        surfacetens(k, igas) = surf_tens(k)
        
        desorption(igas) = 0.5_f
        
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvin(k, igas)
      end do  
    else if (igas .eq. igasal2o3) then
      ! Calculate Kelvin curvature factor for polysulfur interactively with temperature:
      do k = 1, NZ  
        surf_tens(k) = 690._f !Kazasa et al. 1989!900._f ! Dobrovinskaya et al. (2009)
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens(k) / (t(k) * RHO_AL2O3 * RGAS)
        surfacetens(k, igas) = surf_tens(k)
        
        desorption(igas) = 0.5_f
        
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvin(k, igas)
      end do  
    else if (igas .eq. igasco) then
      ! Calculate Kelvin curvature factor for polysulfur interactively with temperature:
      do k = 1, NZ  
        !surf_tens(k) = 9.8_f ! Bierhals J; Ullmann's Encyclopedia of Industrial Chemistry. (2002).
        surf_tens(k) = max(27.77_f*(1._f - t(k)/132.92_f)**1.126_f,1e-5_f) ! Sprow & Prausnitz (1965)
        akelvin(k, igas) = 2._f * gwtmol(igas) * surf_tens(k) / (t(k) * RHO_CO * RGAS)
        surfacetens(k, igas) = surf_tens(k)
        
        ! Not doing condensation of h2So4 on ice, so just set it to the value
        ! for water vapor.
        akelvini(k, igas) = akelvin(k, igas)
      end do  
    else

      ! Condensing gas is not yet configured.
      if (do_print) write(LUNOPRT,*) 'setupgkern::ERROR - invalid igas'
      rc = RC_ERROR
      return
    endif

    ! Molecular free path of condensing gas 
    freep(:,igas)  = 3._f*diffus(:,igas) &
             * sqrt( ( PI*gwtmol_dif(igas) ) / ( 8._f*RGAS*t(:) ) )

    ! Thermal free path of condensing gas
    freept(:,igas) = freep(:,igas)*thcond(:) / &
               ( diffus(:,igas) * rhoa_cgs(:, igas) &
             * ( CP - RGAS/( 2._f*wtmol_air(:) ) ) )
  end do
  
  !write(LUNOPRT,*) 'rmu ', rmu(:)
  !write(LUNOPRT,*) 'thcond ', thcond(:)
  !write(LUNOPRT,*) 'cp ', CP
  !write(LUNOPRT,*) 'pvapl ', pvapl(:,igash2so4)
  !write(LUNOPRT,*) 'surf tens ', surf_tens(:)
  !write(LUNOPRT,*) 'diffus ', diffus(:,igash2so4)
  !write(LUNOPRT,*) 'rlhe ', rlhe(:,igash2so4)


  ! Loop over aerosol groups only (no radius, gas, or spatial dependence).
  do igroup = 1, NGROUP

    ! Use gstickl or gsticki, depending on whether group is ice or not
    if( is_grp_ice(igroup) ) then
      gstick = gsticki
    else
      gstick = gstickl
    endif

    ! Non-spherical corrections (need a reference for these)
    if( (ishape(igroup) .eq. I_SPHERE) .or. (ishape(igroup) .eq. I_FRACTAL))then

      !   Spheres
      cor = 1._f
      phish = 1._f
    else 

      if( ishape(igroup) .eq. I_HEXAGON )then
        
        ! Hexagons
        phish = 6._f/PI*tan(PI/6._f)*( eshape(igroup) + 0.5_f ) &
               * ( PI / ( 9._f*eshape(igroup)*tan(PI/6._f) ) )**(2._f/3._f)

      else if( ishape(igroup) .eq. I_CYLINDER )then

        ! Spheroids
        phish = ( eshape(igroup) + 0.5_f ) &
               * ( 2._f / ( 3._f*eshape(igroup) ) )**(2._f/3._f)
      endif

      if( eshape(igroup) .lt. 1._f )then

        ! Oblate spheroids
        esh1 = 1._f / eshape(igroup)
        a1 = sqrt(esh1**2 - 1._f)
        cor = a1 / asin( a1 / esh1 ) / esh1**(2._f/3._f)
      else 

        ! Prolate spheroids
        a1 = sqrt( eshape(igroup)**2 - 1._f )
            cor = a1 / log( eshape(igroup) + a1 ) &
                / eshape(igroup)**(ONE/3._f)
      endif
    endif

    ! Evaluate growth terms only for particle elements that grow.
    ! particle number concentration element
    ielem = ienconc(igroup)
    
    ! condensing gas is <igas>
    igas = igrowgas(ielem)
    
    ! If the group doesn't grow, but is involved in aerosol
    ! freezing, then the gas properties still need to be calculated.
    if( igas .eq. 0 ) then
      do inuc = 1,nnuc2elem(ielem)
        igas = inucgas(ielem,inuc2elem(inuc, ielem))
	!write(*,*) 'IGAS', inuc, ielem,inuc2elem(inuc,ielem), igas

        if( igas .ne. 0 )then

          do k = 1, NZ

            ! Latent heat of condensing gas 
            if( is_grp_ice(igroup) )then
              rlh = rlhe(k,igas) + rlhm(k,igas)
            else
              rlh = rlhe(k,igas)
            endif

            ! Radius-dependent parameters 
            do i = 1, NBIN

              br = rlow_wet(k,i,igroup)     ! particle bin Boundary Radius

              ! These are Knudsen numbers
              rknudn  = freep(k,igas) / br
              rknudnt = freept(k,igas) / br

              ! These are "lambdas" used in correction for gas kinetic effects.
              rlam  = ( 1.33_f*rknudn  + 0.71_f ) / ( rknudn  + 1._f ) &
                    + ( 4._f*( 1._f - gstick ) ) / ( 3._f*gstick )

              rlamt = ( 1.33_f*rknudnt + 0.71_f ) / ( rknudnt + 1._f ) &
                    + ( 4._f*( 1._f - tstick ) ) / ( 3._f*tstick )

              ! Diffusion coefficient and thermal conductivity modified for
              ! free molecular limit and for particle shape.
              diffus1 = diffus(k,igas)*cor / ( 1._f + rlam*rknudn*cor/phish )
              thcond1 = thcond(k)*cor / ( 1._f + rlamt*rknudnt*cor/phish )

              ! Save the modified thermal conductivity off so it can be used in pheat.
              thcondnc(k,i,igroup,igas) = thcond1
          
              ! Reynolds' number based on particle shape <reyn_shape>
              if( (ishape(igroup) .eq. I_SPHERE) .or. (ishape(igroup) .eq. I_FRACTAL) )then
                reyn_shape = re(k,i,igroup)

              else if( eshape(igroup) .lt. 1._f )then
                reyn_shape = re(k,i,igroup) * ( 1._f + 2._f*eshape(igroup) )

              else
                reyn_shape = re(k,i,igroup) * PI*( 1._f+2._f*eshape(igroup) ) &
                       / ( 2._f*( 1._f + eshape(igroup) ) )
              endif

              ! Particle Schmidt number
              schn = rmu(k) / ( rhoa_cgs(k,igas) * diffus1 )

              ! Prandtl number
              prnum = rmu(k)*CP/thcond1

              ! Ventilation factors <fv> and <ft> from Pruppacher and Klett
              x1 = schn **(ONE/3._f) * sqrt( reyn_shape )
              x2 = prnum**(ONE/3._f) * sqrt( reyn_shape )

              if( is_grp_ice(igroup) )then

                ! Ice crystals
                if( x1 .le. 1._f )then
                  fv = 1._f   + 0.14_f*x1**2
                else
                  fv = 0.86_f + 0.28_f*x1
                endif

                if( x2 .le. 1._f )then
                  ft(k,i,igroup,igas) = 1._f   + 0.14_f*x2**2
                else
                  ft(k,i,igroup,igas) = 0.86_f + 0.28_f*x2
                endif
              else
          
                ! Liquid water drops
                if( x1 .le. 1.4_f  )then
                  fv = 1._f   + 0.108_f*x1**2
                else
                  fv = 0.78_f + 0.308_f*x1
                endif

                if( x2 .le. 1.4_f )then
                  ft(k,i,igroup,igas) = 1._f   + 0.108_f*x2**2
                else
                  ft(k,i,igroup,igas) = 0.78_f + 0.308_f*x2
                endif
              endif
		!write(*,*) 'IZ', k, i, igroup, igas, ft(k,i,igroup,igas)

              ! Growth kernel for particle without radiation or heat conduction at
              ! radius lower boundary [g cm^3 / erg / s]
              gro(k,i,igroup,igas) = 4._f*PI*br &
                            * diffus1*fv*gwtmol_dif(igas) &
                            / ( BK*t(k)*AVG )
      
              ! Coefficient for conduction term in growth kernel [s/g]
              gro1(k,i,igroup,igas) = gwtmol_dif(igas)*rlh**2 &
                    / ( RGAS*t(k)**2*ft(k,i,igroup,igas)*thcond1 ) &
                    / ( 4._f*PI*br )
  
              ! Coefficient for radiation term in growth kernel [g/erg]
              ! (note: no radial dependence).
              if( i .eq. 1 )then
                gro2(k,igroup,igas) = 1._f / rlh
              endif
 
            enddo   ! i=1,NBIN
          enddo    ! k=1,NZ
        endif     ! igas ne 0
      enddo
    else
      do k = 1, NZ

        ! Latent heat of condensing gas 
        if( is_grp_ice(igroup) )then
          rlh = rlhe(k,igas) + rlhm(k,igas)
        else
          rlh = rlhe(k,igas)
        endif

        ! Radius-dependent parameters 
        do i = 1, NBIN

          br = rlow_wet(k,i,igroup)     ! particle bin Boundary Radius

          ! These are Knudsen numbers
          rknudn  = freep(k,igas) / br
          rknudnt = freept(k,igas) / br

          ! These are "lambdas" used in correction for gas kinetic effects.
          rlam  = ( 1.33_f*rknudn  + 0.71_f ) / ( rknudn  + 1._f ) &
                + ( 4._f*( 1._f - gstick ) ) / ( 3._f*gstick )

          rlamt = ( 1.33_f*rknudnt + 0.71_f ) / ( rknudnt + 1._f ) &
                + ( 4._f*( 1._f - tstick ) ) / ( 3._f*tstick )

          ! Diffusion coefficient and thermal conductivity modified for
          ! free molecular limit and for particle shape.
          diffus1 = diffus(k,igas)*cor / ( 1._f + rlam*rknudn*cor/phish )
          thcond1 = thcond(k)*cor / ( 1._f + rlamt*rknudnt*cor/phish )

          ! Save the modified thermal conductivity off so it can be used in pheat.
          thcondnc(k,i,igroup,igas) = thcond1
          
          ! Reynolds' number based on particle shape <reyn_shape>
          if( (ishape(igroup) .eq. I_SPHERE) .or. (ishape(igroup) .eq. I_FRACTAL)  )then
            reyn_shape = re(k,i,igroup)

          else if( eshape(igroup) .lt. 1._f )then
            reyn_shape = re(k,i,igroup) * ( 1._f + 2._f*eshape(igroup) )

          else
            reyn_shape = re(k,i,igroup) * PI*( 1._f+2._f*eshape(igroup) ) &
                   / ( 2._f*( 1._f + eshape(igroup) ) )
          endif

          ! Particle Schmidt number
          schn = rmu(k) / ( rhoa_cgs(k,igas) * diffus1 )

          ! Prandtl number
          prnum = rmu(k)*CP/thcond1

          ! Ventilation factors <fv> and <ft> from Pruppacher and Klett
          x1 = schn **(ONE/3._f) * sqrt( reyn_shape )
          x2 = prnum**(ONE/3._f) * sqrt( reyn_shape )

          if( is_grp_ice(igroup) )then

            ! Ice crystals
            if( x1 .le. 1._f )then
              fv = 1._f   + 0.14_f*x1**2
            else
              fv = 0.86_f + 0.28_f*x1
            endif

            if( x2 .le. 1._f )then
              ft(k,i,igroup,igas) = 1._f   + 0.14_f*x2**2
            else
              ft(k,i,igroup,igas) = 0.86_f + 0.28_f*x2
            endif
          else
          
            ! Liquid water drops
            if( x1 .le. 1.4_f  )then
              fv = 1._f   + 0.108_f*x1**2
            else
              fv = 0.78_f + 0.308_f*x1
            endif

            if( x2 .le. 1.4_f )then
              ft(k,i,igroup,igas) = 1._f   + 0.108_f*x2**2
            else
              ft(k,i,igroup,igas) = 0.78_f + 0.308_f*x2
            endif
          endif
		!write(*,*) 'IZ', k, i, igroup, igas, ft(k,i,igroup,igas)

              ! Growth kernel for particle without radiation or heat conduction at
              ! radius lower boundary [g cm^3 / erg / s]
          gro(k,i,igroup,igas) = 4._f*PI*br &
                        * diffus1*fv*gwtmol_dif(igas) &
                        / ( BK*t(k)*AVG )
      
          ! Coefficient for conduction term in growth kernel [s/g]
          gro1(k,i,igroup,igas) = gwtmol_dif(igas)*rlh**2 &
                / ( RGAS*t(k)**2*ft(k,i,igroup,igas)*thcond1 ) &
                / ( 4._f*PI*br )
  
          ! Coefficient for radiation term in growth kernel [g/erg]
          ! (note: no radial dependence).
          if( i .eq. 1 )then
            gro2(k,igroup,igas) = 1._f / rlh
          endif
 
        enddo   ! i=1,NBIN
      enddo    ! k=1,NZ
    endif
  enddo       ! igroup=1,NGROUP

  ! Return to caller with time-independent particle growth 
  ! parameters initialized.
  return
end
