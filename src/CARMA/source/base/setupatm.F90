! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine setups up parameters related to the atmospheric state. It assumes that the
!! pressure, temperature, and dimensional fields (xc, dx, yc, dy, zc, zl) have already been
!! specified and all state arrays allocated via CARMASTATE_Create().
!! 
!! @author Chuck Bardeen
!! @ version Feb-1995
!! @see CARMASTATE_Create
subroutine setupatm(carma, cstate, rescale, rplanet, rc)

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

  type(carma_type), intent(in)         :: carma       !! the carma object
  type(carmastate_type), intent(inout) :: cstate      !! the carma state object
  logical, intent(in)                  :: rescale     !! rescale the fall velocity for zmet change, this is instead of realculating
  real(kind=f), intent(in)             :: rplanet       !! Planetary radius [cm]
  integer, intent(inout)               :: rc          !! return code, negative indicates failure

  ! Local declarations
  !--
  ! For air viscosity calculations
  ! Air viscosity <rmu> is from Sutherland's equation (using Smithsonian
  !   Meteorological Tables, in which there is a misprint -- T is deg_K, not
  !   deg_C.
  !real(kind=f), parameter :: rmu_0 = 1.8325e-4_f  
  !real(kind=f), parameter :: rmu_t0 = 296.16_f     
  !real(kind=f), parameter :: rmu_c = 120._f  
  !real(kind=f), parameter :: rmu_const = rmu_0 * (rmu_t0 + rmu_c) 


  !PETER:  Values for CO2 viscosity parameters from Viscous Fluid Flow by
  !Frank M. White, 1974.  
 ! real(kind=f), parameter :: rmu_0 = 1.37e-4_f                     
 ! real(kind=f), parameter :: rmu_t0 = 273.15_f                       
 ! real(kind=f), parameter :: rmu_c = 222.2_f                          
 ! real(kind=f), parameter :: rmu_const = rmu_0 * (rmu_t0 + rmu_c)
  
  ! Values for H2 from Wikipedia! Probably have better source somewhere...
  real(kind=f), parameter :: rmu_0 = 0.876e-4_f                     
  real(kind=f), parameter :: rmu_t0 = 293.85_f                       
  real(kind=f), parameter :: rmu_c = 72._f                          
  real(kind=f), parameter :: rmu_const = rmu_0 * (rmu_t0 + rmu_c)
     

    
  integer :: ielem, ibin, i, j, ix, iy, iz, ie, ig, ip, igrp, jgrp, igroup


  ! Calculate the dry air density at each level, using the ideal gas
  ! law. This will be used to calculate zmet.
  rhoa(:) = p(:) / (RGAS/wtmol_air(:) * t(:))

  ! Calculate the dimensions and the dimensional metrics.
  dz(:) = abs(zl(2:NZP1) - zl(1:NZ))
  
  ! Horizontal Metrics
  select case(igridh)
    ! Cartesian
    case (I_CART)
      xmet(:) = 1._f
      ymet(:) = 1._f
    
    case (I_SIG)
      xmet(:) = 1._f
      ymet(:) = 1._f
    
    ! Latitude/Longitude
    case (I_LL)
      xmet(:) = rplanet * DEG2RAD * cos(DEG2RAD * yc(:))    !PETER: REARTH -> RPLANET
      ymet(:) = rplanet * DEG2RAD                           !PETER: REARTH -> RPLANET
      
    case default
      if (do_print) write(LUNOPRT,*) "setupatm:: ERROR - The specified horizontal grid type (", igridh, &
        ") is not supported."
      rc = -1
  end select

  
  ! Put the fall velocity back into cgs units, so that we can determine
  ! new metrics and then scale it back. This is optional and is done instead
  ! of recalculating everything from scratch to improve performance.
  if (rescale .and. (igridv /= I_CART) .and.(igridv /= I_LOGP)) then
    do ibin = 1, NBIN
      do igroup = 1, NGROUP
        vf(:, ibin, igroup)   = vf(:, ibin, igroup)  * zmetl(:)
        dkz(:, ibin, igroup)  = dkz(:, ibin, igroup) * (zmetl(:)**2)
        ekz(:)  = ekz(:) * (zmetl(:)**2) !DPOW CHECK
      end do
    end do
  end if

 
  ! Vertical Metrics
  select case(igridv)
    ! Cartesian
    case (I_CART)
      zmet(:) = 1._f 
      
    ! Log-Pressure
    case (I_LOGP)
      zmet(:) = t(:)/t0
        
      
    ! Sigma
    case (I_SIG)
      zmet(:) = abs(((pl(1:NZ) - pl(2:NZP1)) &
	/ (zl(1:NZ) - zl(2:NZP1))) / (grav(:) * rhoa(:)))
        
    ! Hybrid
    case (I_HYBRID)
      zmet(:) = abs(((pl(1:NZ) - pl(2:NZP1)) &
	/ (zl(1:NZ) - zl(2:NZP1))) / &
!        (grav(:) * (RPLANET/(RPLANET+zc(:)))**2._f * rhoa(:)))
        (grav(:) * rhoa(:)))
        
    case default
      if (do_print) write(LUNOPRT,*) "setupatm:: ERROR - The specified vertical grid type (", igridv, &
        ") is not supported."
      rc = -1
  end select
 
  ! Interpolate the z metric to the grid box edges.
  if (NZ == 1) then
    zmetl(:) = zmet(1)
  else
    
    ! Extrpolate the top and bottom.
    zmetl(1)    = zmet(1)  + (zmet(2) - zmet(1))    &
	 / (zc(2) - zc(1))     * (zl(1) - zc(1))
    zmetl(NZP1) = zmet(NZ) + (zmet(NZ) - zmet(NZ-1)) &
	/ (zc(NZ) - zc(NZ-1)) * (zl(NZP1) - zc(NZ))
    
    ! Interpolate the middles.
    if (NZ > 2) then
      do iz = 2, NZ
        zmetl(iz) = zmet(iz-1) + (zmet(iz) - zmet(iz-1)) &
	/ (zc(iz) - zc(iz-1)) * (zl(iz) - zc(iz-1))
      end do
    end if
  end if

 
  ! Determine the z metrics at the grid box edges and then use this to put the
  ! fall velocity back into /x/y/z units.
  if (rescale .and. (igridv /= I_CART) .and.(igridv /= I_LOGP) ) then !DPOW check this for logp
    do ibin = 1, NBIN
      do igroup = 1, NGROUP
        vf(:, ibin, igroup)   = vf(:, ibin, igroup)  / zmetl(:)
        dkz(:, ibin, igroup)  = dkz(:, ibin, igroup) / (zmetl(:)**2)
        ekz(:) = ekz(:) / (zmetl(:)**2) !DPOW check if I should do this here!
      end do
    end do
  end if
  
  
  ! Scale the density into the units carma wants (i.e. /x/y/z)     
  rhoa(:) = rhoa(:) * xmet(:) * ymet(:) * zmet(:)

  ! Use the pressure difference across the cell and the fact that the 
  ! atmosphere is hydrostatic to caclulate an average density in the
  ! grid box.
!  rhoa_wet(:) = abs((pl(2:NZP1) - pl(1:NZ))) / (grav(:) * (RPLANET/(RPLANET+zc(:)))**2._f)
!DPOW adjustment
  rhoa_wet(:) = rhoa(:)
   
  !rhoa_wet(:) = abs((pl(2:NZP1) - pl(1:NZ))) / grav(:)
  !rhoa_wet(:) = (rhoa_wet(:) * xmet(:) * ymet(:)) / dz(:)
  !write(*,*) 'rhoa_wet', rhoa_wet(:)
  !write(*,*) 'rhoa',rhoa(:)
  !write(*,*) 'diff', (rhoa_wet(:) - rhoa(:)) / rhoa(:)
  ! Calculate the thermal properties of the atmosphere.
  rmu(:)     = rmu_const / ( t(:) + rmu_c ) * (t(:) / rmu_t0 )**1.5_f
!  thcond(:)  = (5.69_f + .017_f*(t(:) - T0)) * 4.186e2_f
!  thcond(:)  = (1520._f + 7.65_f*(t(:) - T0))                  !PETER: thermal conductivity of C02 from handbook of chemistry and physics
  ! H2 atmosphere from Thermophysical Properties of Fluids section of the Handbook of Chemistry and Physics, by Eric W. Lemmon.
  thcond(:) = 100._f * (71.4857_f + 0.3912_f * t(:) + 3.1607e-5_f * t(:)**2._f)

end subroutine
