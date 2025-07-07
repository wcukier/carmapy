!! The CARMAGAS module contains configuration information about a gas used by CARMA.
!!
!!  @version May-2009 
!!  @author  Chuck Bardeen 
module carmagas_mod

  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_planet_mod
  use carma_condensate_mod
  use carma_types_mod

  ! CARMA explicitly declares all variables. 
  implicit none

  ! All CARMA variables and procedures are private except those explicitly declared to be public.
  private

  ! Declare the public methods.
  public CARMAGAS_Create
  public CARMAGAS_Destroy
  public CARMAGAS_Get
  public CARMAGAS_Print

contains

  !! Defines a gas used by CARMA for nucleation and growth of cloud and 
  !! aerosol particles.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMA_AddGas
  !! @see CARMAGAS_Destroy
  subroutine CARMAGAS_Create(carma, igas, name, wtmol, ivaprtn, icomposition,  &
    rho_cond, surften_0, coldia,  vp_offset, vp_tcoeff, rc,&
    is_type3,  surften_slope, vp_metcoeff, vp_logpcoeff, shortname, &
    dgc_threshold, ds_threshold, wtmol_dif)
    type(carma_type), intent(inout)       :: carma           !! the carma object
    integer, intent(in)                   :: igas            !! the gas index
    character(*), intent(in)              :: name            !! the gas name, maximum of 255 characters
    real(kind=f), intent(in)              :: wtmol           !! the gas molecular weight [g/mol]
    integer, intent(in)                   :: ivaprtn         !! vapor pressure routine for this gas
    integer, intent(in)                   :: icomposition    !! gas compound specification
    real(kind=f), intent(in)              :: rho_cond        !! density of condenstate [g/cm^3] WC
    real(kind=f), intent(in)              :: surften_0       !! surface tension at 0 K (assuming linear extrapolation) [dyne/cm] WC
    real(kind=f), intent(in)              :: coldia          !! collision diameter [cm] WC
    real(kind=f), intent(in)              :: vp_offset       !! vapor pressure = 10**(vp_offset - vp_tcoeff/Temperature - vp_metcoeff*metallicity - vp_logpcoeff * log10(Pressure)) WC
    real(kind=f), intent(in)              :: vp_tcoeff       !! vapor pressure = 10**(vp_offset - vp_tcoeff/Temperature - vp_metcoeff*metallicity - vp_logpcoeff * log10(Pressure)) WC
    integer, intent(out)                  :: rc              !! return code, negative indicates failure
    integer, optional, intent(in)         :: is_type3        !! 1 for type III reations, 0 otherwise (see Helling and Woitke 2006) WC
    real(kind=f), optional, intent(in)    :: surften_slope   !! slope of surface tension wrt temperature (surften = surften0 - surften_slope (t - 273.15)) [dyne/cm/K]
    real(kind=f), optional, intent(in)    :: vp_metcoeff     !! vapor pressure = 10**(vp_offset - vp_tcoeff/Temperature - vp_metcoeff*metallicity - vp_logpcoeff * log10(Pressure))
    real(kind=f), optional, intent(in)    :: vp_logpcoeff    !! vapor pressure = 10**(vp_offset - vp_tcoeff/Temperature - vp_metcoeff*metallicity - vp_logpcoeff * log10(Pressure))
    character(*), optional, intent(in)    :: shortname       !! the gas shortname, maximum of 6 characters
    real(kind=f), optional, intent(in)    :: dgc_threshold   !! convergence criteria for gas concentration [0 : off; > 0 : percentage change]
    real(kind=f), optional, intent(in)    :: ds_threshold    !! convergence criteria for gas saturation [0 : off; > 0 : percentage change; < 0 : amount past 0 crossing]
    real(kind=f), optional, intent(in)    :: wtmol_dif       !! gas molecular weight in vapor phase [g/mol]
    

    write(*, *) ivaprtn
    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough gases allocated.
    if (igas > carma%f_NGAS) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAGAS_GetCreate:: ERROR - The specifed gas (", &
        igas, ") is larger than the number of gases (", carma%f_NGAS, ")."
      rc = RC_ERROR
      return
    end if

    ! Save off the settings.
    carma%f_gas(igas)%f_name         = name
    carma%f_gas(igas)%f_wtmol        = wtmol
    carma%f_gas(igas)%f_ivaprtn      = ivaprtn
    carma%f_gas(igas)%f_icomposition = icomposition
    carma%f_gas(igas)%f_rho_cond     = rho_cond
    carma%f_gas(igas)%f_surften_0    = surften_0
    carma%f_gas(igas)%f_coldia       = coldia
    carma%f_gas(igas)%f_vp_offset    = vp_offset
    carma%f_gas(igas)%f_vp_tcoeff    = vp_tcoeff



    ! Defaults for optional parameters
    carma%f_gas(igas)%f_shortname       = ""
    carma%f_gas(igas)%f_dgc_threshold   = 0._f
    carma%f_gas(igas)%f_ds_threshold    = 0._f
    carma%f_gas(igas)%f_wtmol_dif       = wtmol
    carma%f_gas(igas)%f_is_type3        = 0
    carma%f_gas(igas)%f_surften_slope    = 0._f
    carma%f_gas(igas)%f_vp_metcoeff     = 0._f
    carma%f_gas(igas)%f_vp_logpcoeff    = 0._f

    ! Set optional parameters.
    if (present(shortname))     carma%f_gas(igas)%f_shortname      = shortname
    if (present(dgc_threshold)) carma%f_gas(igas)%f_dgc_threshold  = dgc_threshold
    if (present(ds_threshold))  carma%f_gas(igas)%f_ds_threshold   = ds_threshold
    if (present(wtmol_dif))  carma%f_gas(igas)%f_wtmol_dif         = wtmol_dif
    if (present(is_type3))    carma%f_gas(igas)%f_is_type3        = is_type3
    if (present(surften_slope))    carma%f_gas(igas)%f_surften_slope    = surften_slope
    if (present(vp_metcoeff))    carma%f_gas(igas)%f_vp_metcoeff     = vp_metcoeff
    if (present(vp_logpcoeff))    carma%f_gas(igas)%f_vp_logpcoeff    = vp_logpcoeff



    return
  end subroutine CARMAGAS_Create
    

  !! Deallocates the memory associated with a CARMAGAS object.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMAGAS_Create
  subroutine CARMAGAS_Destroy(carma, igas, rc)
    type(carma_type), intent(inout)    :: carma         !! the carma object
    integer, intent(in)                :: igas          !! the gas index
    integer, intent(out)               :: rc            !! return code, negative indicates failure

    ! Assume success.
    rc = RC_OK
    
    ! Make sure there are enough gases allocated.
    if (igas > carma%f_NGAS) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAGAS_Destroy:: ERROR - The specifed gas (", &
        igas, ") is larger than the number of gases (", carma%f_NGAS, ")."
      rc = RC_ERROR
      return
    end if

    return
  end subroutine CARMAGAS_Destroy


  !! Gets information about a gas.
  !!
  !! The group name and other properties are available after a call to
  !! CARMAGAS_Create().
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMAGAS_Create
  !! @see CARMA_GetGas
  subroutine CARMAGAS_Get(carma, igas, rc, name, shortname, wtmol, ivaprtn, icomposition, dgc_threshold, ds_threshold)
    type(carma_type), intent(in)                :: carma         !! the carma object
    integer, intent(in)                         :: igas          !! the gas index
    integer, intent(out)                        :: rc            !! return code, negative indicates failure
    character(len=*), optional, intent(out)     :: name          !! the gas name
    character(len=*), optional, intent(out)     :: shortname     !! the gas short name
    real(kind=f), optional, intent(out)         :: wtmol         !! the gas molecular weight [g/mol]
    integer, optional, intent(out)              :: ivaprtn       !! vapor pressure routine for this gas
    integer, optional, intent(out)              :: icomposition  !! gas compound specification
    real(kind=f), optional, intent(out)         :: dgc_threshold !! convergence criteria for gas concentration [fraction]
    real(kind=f), optional, intent(out)         :: ds_threshold  !! convergence criteria for gas saturation [fraction]

    ! Assume success.
    rc = RC_OK

    ! Make sure there are enough gases allocated.
    if (igas > carma%f_NGAS) then
      if (carma%f_do_print) write(carma%f_LUNOPRT, *) "CARMAGAS_Get:: ERROR - The specifed gas (", &
        igas, ") is larger than the number of gases (", carma%f_NGAS, ")."
      rc = RC_ERROR
      return
    end if

    ! Return any requested properties of the group.
    if (present(name))         name         = carma%f_gas(igas)%f_name
    if (present(shortname))    shortname    = carma%f_gas(igas)%f_shortname
    if (present(wtmol))        wtmol        = carma%f_gas(igas)%f_wtmol
    if (present(ivaprtn))      ivaprtn      = carma%f_gas(igas)%f_ivaprtn
    if (present(icomposition)) icomposition = carma%f_gas(igas)%f_icomposition
    if (present(dgc_threshold)) dgc_threshold = carma%f_gas(igas)%f_dgc_threshold
    if (present(ds_threshold)) ds_threshold = carma%f_gas(igas)%f_ds_threshold
        
    return
  end subroutine CARMAGAS_Get
  
  
  !! Prints information about a gas.
  !!
  !! @author  Chuck Bardeen
  !! @version May-2009
  !!
  !! @see CARMAGAS_Get
  subroutine CARMAGAS_Print(carma, igas, rc)
    type(carma_type), intent(in)              :: carma         !! the carma object
    integer, intent(in)                       :: igas          !! the gas index
    integer, intent(out)                      :: rc            !! return code, negative indicates failure
    
    ! Local variables
    character(len=CARMA_NAME_LEN)             :: name          !! name
    character(len=CARMA_SHORT_NAME_LEN)       :: shortname     !! shortname
    real(kind=f)                              :: wtmol         !! molecular weight (g/mol)
    integer                                   :: ivaprtn       !! vapor pressure routine for this gas
    integer                                   :: icomposition  !! gas compound specification
    real(kind=f)                              :: dgc_threshold !! convergence criteria for gas concentration [fraction]
    real(kind=f)                              :: ds_threshold  !! convergence criteria for gas saturation [fraction]

    ! Assume success.
    rc = RC_OK

    ! Test out the Get method.
    if (carma%f_do_print) then
      call CARMAGAS_Get(carma, igas, rc, name=name, shortname=shortname, wtmol=wtmol, &
                        ivaprtn=ivaprtn, icomposition=icomposition)
      if (rc < RC_OK) return

    
      write(carma%f_LUNOPRT,*) "    name          : ", trim(name)
      write(carma%f_LUNOPRT,*) "    shortname     : ", trim(shortname)
      write(carma%f_LUNOPRT,*) "    wtmol         : ", wtmol, " (g/mol)"
      write(carma%f_LUNOPRT,*) "    dgc_threshold : ", dgc_threshold
      write(carma%f_LUNOPRT,*) "    ds_threshold  : ", ds_threshold

      select case(ivaprtn)
        case (I_VAPRTN_H2O_BUCK1981)
          write(carma%f_LUNOPRT,*) "    ivaprtn       :    Buck [1981]"
        case (I_VAPRTN_H2O_MURPHY2005)
          write(carma%f_LUNOPRT,*) "    ivaprtn       :    Murphy & Koop [2005]"
        case default
          write(carma%f_LUNOPRT,*) "    ivaprtn       :    unknown, ", ivaprtn
      end select

      select case(icomposition)
        case (I_GCOMP_H2O)
          write(carma%f_LUNOPRT,*) "    icomposition  :    H2O"
        case default
          write(carma%f_LUNOPRT,*) "    icomposition  :    unknown, ", icomposition
      end select
    end if
    
    return
  end subroutine CARMAGAS_Print
end module
