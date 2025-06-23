!! @author Diana Powell
!! @version Oct-2019

program carma_day
  implicit none

  call test_day()

  write(*,*) "Done"
end program

subroutine test_day()
  use carma_precision_mod
  use carma_constants_mod
  use carma_enums_mod
  use carma_types_mod
  use carmaelement_mod
  use carmagroup_mod
  use carmastate_mod
  use carmagas_mod
  use carmasolute_mod
  use carma_mod

  implicit none

  integer, parameter    :: io = 90

  integer    :: NZ    
  integer    :: NZP1        
  integer    :: NELEM 
  integer    :: NBIN    
  integer    :: NGROUP     
  integer    :: NSOLUTE   
  integer   :: NGAS        
  integer    :: NWAVE  
  integer    :: NLONGITUDE 
  integer   :: NGROWTH, NNUC, NCOAG
  integer   :: IS_2D
  integer   :: igridv

  real(kind=f)   :: dtime 
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 200._f
  real(kind=f), parameter   :: zmin   = 0._f
  real(kind=f)   :: rplanet 
  real(kind=f)   :: velocity_avg 

 ! integer, parameter    :: irestart     = 1  ! =1 to restart
  integer    :: irestart    ! =1 to restart
!  integer, parameter    :: idiag     = 1  ! =1 to output diagnostic
  integer    :: idiag   ! =1 to output diagnostic
  integer    :: iappend
  integer    :: iskip        ! Output every iskip steps; no steps skipped if = 1
  integer    :: nstep        

  integer, parameter        :: I_KCL       = 1
  integer, parameter        :: I_ZNS       = 2
  integer, parameter        :: I_NA2S      = 3
  integer, parameter        :: I_MNS       = 4
  integer, parameter        :: I_CR        = 5
  integer, parameter        :: I_MG2SIO4   = 6
  integer, parameter        :: I_FE        = 7
  integer, parameter        :: I_TIO2      = 8
  integer, parameter        :: I_AL2O3     = 9

  type(carma_type), target            :: carma
  type(carma_type), pointer           :: carma_ptr
  type(carmastate_type)               :: cstate
  integer                             :: rc = 0

  real(kind=f), allocatable   :: xc(:)
  real(kind=f), allocatable   :: dx(:)
  real(kind=f), allocatable   :: yc(:)
  real(kind=f), allocatable   :: dy(:)
  real(kind=f), allocatable   :: zc(:)
  real(kind=f), allocatable   :: zl(:)
  real(kind=f), allocatable   :: dz_test(:)
  real(kind=f), allocatable   :: p(:)
  real(kind=f), allocatable   :: pl(:)
  real(kind=f), allocatable   :: t(:)
  real(kind=f), allocatable   :: rho_atm_cgs(:)

  real(kind=f), allocatable   :: mmr(:,:,:)
  real(kind=f), allocatable   :: numden(:,:,:)
  real(kind=f), allocatable   :: mmr_gas(:,:)
  real(kind=f), allocatable   :: mmr_gas_old(:,:)
  real(kind=f), allocatable   :: satliq(:,:)
  real(kind=f), allocatable   :: satice(:,:)
  real(kind=f), allocatable   :: satliq_old(:,:)
  real(kind=f), allocatable   :: satice_old(:,:)
  real(kind=f), allocatable   :: svpliq(:,:)
  real(kind=f), allocatable   :: wtpct(:,:)
  real(kind=f), allocatable   :: gflux(:,:)
  real(kind=f), allocatable   :: pflux(:,:,:)
  real(kind=f), allocatable   :: winds(:)
  ! real(kind=f), allocatable   :: ekz(:)
  real(kind=f), allocatable   :: prodrate(:,:,:)
  real(kind=f), allocatable   :: prodrate_mass(:,:,:)
  real(kind=f), allocatable   :: prodrate_gas(:,:)
  real(kind=f), allocatable   :: totmass(:)

  real(kind=f), allocatable   :: ftopp(:,:)
  real(kind=f), allocatable   :: fbotp(:,:)
  real(kind=f), allocatable   :: pctop(:,:)
  real(kind=f), allocatable   :: pcbot(:,:)
  real(kind=f), allocatable   :: gctop(:)
  real(kind=f), allocatable   :: gcbot(:)
  real(kind=f), allocatable   :: ftopg(:)
  real(kind=f), allocatable   :: fbotg(:)

  real(kind=f), allocatable   :: gasprod(:,:)
  real(kind=f), allocatable   :: rhompe(:,:,:)
  real(kind=f), allocatable   :: growpe(:,:,:)
  real(kind=f), allocatable   :: evappe(:,:,:)
  real(kind=f), allocatable   :: growlg(:,:,:)
  real(kind=f), allocatable   :: evaplg(:,:,:)

  real(kind=f), allocatable   :: zsubsteps(:)

  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: rlow(:)
  real(kind=f), allocatable   :: rup(:)
  real(kind=f), allocatable   :: dr(:)
  real(kind=f), allocatable   :: rmass(:,:)

  real(kind=f)          :: lat
  real(kind=f)          :: lon

  integer               :: i
  integer               :: j
  integer               :: istep, istep_old, itype
  integer               :: ifrom, ito, ievp2elem, is_het
  integer               :: igas
  integer               :: igroup
  integer               :: ielem
  integer               :: ibin
  integer               :: ibinm
  integer               :: iz
  integer               :: icomposition, iroutine
  integer, parameter    :: lunerr = 48
  integer, parameter    :: lun = 42
  integer, parameter    :: lunp = 43
  integer, parameter    :: lunf = 44
  integer, parameter    :: lunfp = 45
  integer, parameter    :: lunres = 46
  integer, parameter    :: lundiagn = 47
  integer, parameter    :: gas_in = 61
  integer               :: nsubsteps
  integer               :: binmultiple
  integer               :: lastsub = 0

  integer, parameter    :: lun_atm = 49
  integer, parameter    :: lun_kzz = 50

  integer, parameter    :: lunrates = 55
  integer, parameter    :: lunratesp = 56

  real(kind=f)          :: nretries
  real(kind=f)          :: lastret = 0._f
  real(kind=f)          :: t_orig

  real(kind=f)          :: t_in
  real(kind=f)          :: p_in
  real(kind=f)          :: pl_in
  real(kind=f)          :: z_in
  real(kind=f)          :: zl_in
  real(kind=f)          :: mu_in
  real(kind=f)          :: g_in
  real(kind=f)          :: kzz_in
  real(kind=f)          :: tio2_in
  real(kind=f)          :: t0_in

  real(kind=f)          :: met	! Metallicity = [Fe/H]

  real(kind=f)          :: startcd
  real(kind=f)          :: endcd
  real(kind=f)          :: inputrate
  real(kind=f)          :: vertpartflux
  real(kind=f)          :: vertgasflux

  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat, wtmol, wtmol_dif, mucos

  real(kind=f)          :: wtmol_air_set, grav_set
  real(kind=f)          :: rho 


  character(30)      	:: name
  character(30)  :: type_spec
  character(30)      	:: gname
  character(len=3)	:: fileprefix
  character(len=5)	:: temp = 'temp_'
  character(len=5)	:: flux = 'flux_'
  character(len=5)	:: diag = 'diag_'
  character(len=6)	:: rates = 'rates_'
  character(len=4)	:: filesuffix = '.txt'
  character(len=4)	:: filesuffix_restart = '.dat'
  character(len=100)	:: filename_restart
  character(len=100)	:: filename
  character(len=100)  :: nml_file = "inputs/input.nml"
  character(len=100)  :: gas_input_file
  character(len=100)  :: centers_file
  character(len=100)  :: levels_file
  character(len=100)  :: temps_file
  character(len=100)  :: groups_file
  character(len=100)  :: elements_file
  character(len=100)  :: gases_file
  character(len=100)  :: growth_file
  character(len=100)  :: nuc_file
  character(len=100)  :: coag_file
  character(len=20)   :: file_pos

  real(kind=f)          :: distance_btwn_elements, circumference, rotation_counter, slope, intercept
  real(kind=f)          :: current_distance, closeto_temp_profile, num_steps_btwn, current_step, RPLANET_DAT

  namelist / io_files / filename, filename_restart, fileprefix, gas_input_file, centers_file, levels_file,temps_file, groups_file, elements_file, gases_file, growth_file, nuc_file, coag_file
  namelist / physical_params / wtmol_air_set, grav_set, rplanet, velocity_avg
  namelist / input_params / NZ, NELEM, NGROUP, NGAS, NBIN, NSOLUTE, NWAVE, NLONGITUDE, irestart, idiag, iskip, nstep, dtime, NGROWTH, NNUC, NCOAG, IS_2D, igridv, iappend


  real(kind=f), allocatable ::tempr(:), pre(:), prel(:), alt(:), altl(:), wtmol_air(:), grav(:), ekz(:), ekzl(:), wtmol_gas(:)
  real(kind=f), allocatable ::temp_equator(:, :), p_equator_center(:), p_equator_level(:), velocity(:), longitudes(:)
  integer, allocatable :: elem2group(:)



  write(*,*) ""

  open(unit=10, file=nml_file, status='old')
    read(10, nml=input_params)
    read(10, nml=io_files)
    read(10, nml=physical_params)
  close(10)

  NZP1 = NZ + 1

  file_pos = "asis"
  if ((irestart .eq. 1).and.(iappend .eq. 1)) file_pos = "append"



  allocate(tempr(NZ), pre(NZ), prel(NZP1), alt(NZ), altl(NZP1), wtmol_air(NZ), grav(NZ), ekz(NZP1), ekzl(NZP1), wtmol_gas(NGAS))
  allocate(temp_equator(NZ, NLONGITUDE), p_equator_center(NZ), p_equator_level(NZP1), velocity(NLONGITUDE), longitudes(NLONGITUDE))
  allocate(elem2group(NELEM))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(12, file = centers_file)
  read(12, *)
  do i = 1, NZ
 	  read(12,*) alt(i), pre(i)
  end do
  close(12)

  open(12, file = levels_file) 
  read(12, *)

  do i=1, NZP1
    read(12,*) altl(i), prel(i), ekz(i)
  end do

  close(12)


  open(12, file=temps_file)
  do i=1, NZ
    if (IS_2D .eq. 1) then
      read (12, *) temp_equator(i, :)
      tempr(i) = temp_equator(i, 1)
    else
      read (12, *) tempr(i)
    end if
  end do

  close(12)





  wtmol_air(:) =wtmol_air_set 
  grav(:) = grav_set 
  met = 1._f

  t0_in = MAXVAL(tempr)

  circumference = 2 *PI * rplanet !cm
  distance_btwn_elements = circumference/NLONGITUDE
  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Open the output text file
  open(unit=lun,file = fileprefix // filename(1:len_trim(filename)) // filesuffix, status="unknown", position=file_pos)
  open(unit=lunf,file = fileprefix // flux // filename(1:len_trim(filename)) // filesuffix, status="unknown", position=file_pos)
  open(unit=lunrates,file = fileprefix // rates // filename(1:len_trim(filename)) // filesuffix, status="unknown", position=file_pos)

  ! Allocate the arrays that we need for the model
  allocate(xc(NZ), dx(NZ), yc(NZ), dy(NZ), &
           zc(NZ), zl(NZP1), dz_test(NZ), p(NZ), pl(NZP1), &
           t(NZ), rho_atm_cgs(NZ))
  allocate(mmr(NZ,NELEM,NBIN))
  allocate(numden(NZ,NELEM,NBIN))
  allocate(mmr_gas(NZ,NGAS))
  allocate(mmr_gas_old(NZ,NGAS))
  allocate(satliq(NZ,NGAS))
  allocate(satice(NZ,NGAS))
  allocate(satliq_old(NZ,NGAS))
  allocate(satice_old(NZ,NGAS))
  allocate(svpliq(NZ,NGAS))
  allocate(wtpct(NZ,NGAS))
  allocate(gflux(NZP1,NGAS))
  allocate(pflux(NZP1,NBIN,NELEM))
  allocate(winds(NZ))
!  allocate(ekz(NZP1))
  allocate(prodrate(NZ,NBIN,NELEM))
  allocate(prodrate_mass(NZ,NBIN,NELEM))
  allocate(prodrate_gas(NZ,NGAS))
  allocate(totmass(NZ))
  allocate(zsubsteps(NZ))
  allocate(r(NBIN))
  allocate(rlow(NBIN))
  allocate(rup(NBIN))
  allocate(dr(NBIN))
  allocate(rmass(NBIN,NGROUP))
  allocate(ftopp(NBIN,NELEM))
  allocate(fbotp(NBIN,NELEM))
  allocate(pctop(NBIN,NELEM))
  allocate(pcbot(NBIN,NELEM))
  allocate(gctop(NGAS))
  allocate(gcbot(NGAS))
  allocate(ftopg(NGAS))
  allocate(fbotg(NGAS))
  allocate(gasprod(NZ,NGAS))
  allocate(rhompe(NZ,NBIN,NELEM))
  allocate(growpe(NZ,NBIN,NELEM))
  allocate(evappe(NZ,NBIN,NELEM))
  allocate(growlg(NZ,NBIN,NGROUP))
  allocate(evaplg(NZ,NBIN,NGROUP))





  

  ! Define the particle-grid extent of the CARMA test
  write(*,*) "Create CARMA Object ..."

  if (idiag .eq. 1) then 
    call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=6, lundiag=lundiagn)
  else
    call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=6)
  end if


  if (rc < 0) stop "    *** FAILED in CARMA_Create ***"

	carma_ptr => carma

  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define the groups
  write(*,*) "  Add Group(s) ..."
  write(*,*) " "

  open(10, file = groups_file)
  read(10, *)
  do i = 1, NGROUP
    read(10, *) name, rmin

    write(*,*) "Add " //TRIM(name)//"..."

    call CARMAGROUP_Create(carma, i, name, rmin, 2._f, I_SPHERE, 1._f, &
    .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)

    if (rc < 0) stop "    *** FAILED ***"

  enddo
  close(10)

  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define the element
  write(*,*) "  Add Element(s) ..."
  write(*,*) " "

  open(10, file=elements_file)
  read(10, *)
  do i=1, NELEM
    read(10, *) igroup, name, rho, type_spec, icomposition


    write(*,*) "Add "// trim(name)// "..."

    if(trim(type_spec) == "Volatile") then
      itype = I_VOLATILE
    else if(trim(type_spec) == "Core Mass") then
      itype = I_COREMASS
    else
      stop "invalid element type"
    endif
      call CARMAELEMENT_Create(carma, i, igroup, name, rho, itype, icomposition, rc)
      elem2group(i) = igroup
      if (rc < 0) stop "    *** FAILED ***"

  enddo
  close(10)

  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define the solutes
  write(*,*) "  Add Solute(s) ..."
  write(*,*) " "

  write(*,*) "Add H2SO4 Solute ..."
  call CARMASOLUTE_Create(carma, 1, "Sulfuric Acid", 2, WTMOL_H2SO4,RHO_H2SO4, rc)
  if (rc < 0) stop "    *** FAILED ***"

  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define the gases
  write(*,*) "  Add Gas(es) ..."
  write(*,*) " "

  open(10, file=gases_file)
  read(10, *)
  do i=1, NGAS

    read(10, *) name, wtmol, iroutine, icomposition, wtmol_dif

    write(*,*) "Add "//trim(name) //" ..."
    if (wtmol == wtmol_dif) then
      call CARMAGAS_Create(carma, i, name, wtmol, INT(iroutine), icomposition, rc)
    else
      call CARMAGAS_Create(carma, i, name, wtmol, INT(iroutine), icomposition, rc, wtmol_dif=wtmol_dif)
    endif
    if (rc < 0) stop "    *** FAILED ***"

    wtmol_gas(i) = wtmol_dif
  enddo
  close(10)
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Setup the CARMA processes to exercise growth, nucleation, and coagulation.

  write(*,*) "  Add Nucleation ..."
  open(10, file=nuc_file)
  read(10, *)
  do i = 1, NNUC
    read(10, *) ifrom, ito, is_het, igas, ievp2elem, mucos 
    if (is_het == 1) then
        call CARMA_AddNucleation(carma, ifrom, ito, I_HETGEN, 0._f, rc, igas=igas, ievp2elem=ievp2elem, mucos=mucos )
    else
      call CARMA_AddNucleation(carma, ifrom, ito, I_HOMGEN, 0._f, rc, igas=igas)
    endif
    if (rc < 0) stop "    *** FAILED ***"
  enddo
  close(10)

  write(*,*) "  Add Growth ..."
  open(10, file=growth_file)
  read(10, *)
  do i = 1, NNUC
    read(10, *) ito, igas

    call CARMA_AddGrowth(carma, ito, igas, rc)
    if (rc < 0) stop "    *** FAILED ***"

  enddo
  close(10)

 write(*,*) "  Add Coagulation ..."
 open(10, file=coag_file)
 read(10, *)
 do i = 1, NCOAG
   read(10, *) igroup

   call CARMA_AddCoagulation(carma, igroup, igroup, igroup, I_COLLEC_FUCHS, rc) 
   if (rc < 0) stop "    *** FAILED ***"
   
 enddo
 close(10)



  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Setup the CARMA processes to exercise
  write(*,*) "Initialize CARMA ..."

  call CARMA_Initialize(carma, rc, do_cnst_rlh =.FALSE., do_coag=.FALSE., do_fixedinit=.TRUE., do_grow=.TRUE., &
                        do_explised=.FALSE., do_substep=.TRUE., do_print_init=.TRUE., &
                        do_vdiff=.TRUE., do_vtran=.TRUE., maxsubsteps=10, maxretries=20, &
                        itbnd_pc=I_FLUX_SPEC, ibbnd_pc=I_FIXED_CONC, itbnd_gc=I_FLUX_SPEC, ibbnd_gc=I_FIXED_CONC)
  if (rc < 0) stop "    *** FAILED CARMA_Initialize ***"

  write(*,*) " "

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! For simplicity of setup, do a case with Cartesian coordinates,
  ! which are specified in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom
  ! of the model (e.g. z = 0), while for sigma and hybrid coordinates
  ! the first level is the top of the model.

  ! Setup model atmosphere

  lat = 45.0_f
  lon = -105.0_f

  ! Horizonal centers
  dx(:) = deltax
  xc(:) = dx(:) / 2._f
  dy(:) = deltay
  yc(:) = dy(:) / 2._f

  ! Vertical center
  do i = 1, NZ
    zc(i) = alt(i)
    t(i) = tempr(i)
    p(i) = pre(i)
  end do

  write(*,'(a6, 3a12)') "level", "zc", "p", "t"
  write(*,'(a6, 3a12)') "", "(m)", "(Pa)", "(K)"
  do i = 1, NZ
    write(*,'(i6,3f12.3)') i, zc(i), p(i), t(i)
  end do

  ! Vertical edge
  do i = 1, NZP1
    zl(i) = altl(i)
    pl(i) = prel(i)
  end do

  write(*,*) ""
  write(*,'(a6, 2a12)') "level", "zl", "pl"
  write(*,'(a6, 2a12)') "", "(m)", "(Pa)"
  do i = 1, NZP1
    write(*,'(i6,2f12.3)') i, zl(i), pl(i)
  end do

  ! ekz(:) = ekzl(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Setup up a mass mixing ratio of H20 and H2SO4 vapor

  mmr_gas(:,:) = 1e-50_f

  mmr_gas_old = mmr_gas

  ! Setup an intial mmr for condensation nuclei
  mmr(:,:,:) = 0._f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Write output for the test
  if (idiag .eq. 1) then
    open(unit=lundiagn,file = fileprefix // diag // filename(1:len_trim(filename)) // filesuffix, status="unknown")
    write(lundiagn,*) 'PART 0: HEADER'
    write(lundiagn,*) 'SIMULATION TYPE:'
    write(lundiagn,*) 'VENUS'
    write(lundiagn,*) 'DIMENSIONS:'
    write(lundiagn,'(6A9)') 'NZ', 'NGROUP', 'NELEM', 'NBIN', 'NGAS', 'NSTEP'
    write(lundiagn,'(6I9)') NZ, NGROUP, NELEM, NBIN, NGAS, nstep + 1
    write(lundiagn,*) 'GROUP INFORMATION:'
    write(lundiagn,'(A13,2A10,A17,A14,A19,A15)') 'IGROUP', 'IBIN', 'R', 'MASS', 'dR', 'R_LOWBOUND', 'R_UPBOUND' 
    write(lundiagn,'(A10,A12,A15,A11,A19,2A15)') '','','(microns)', '(g)','(microns)','(microns)','(microns)'
  end if 





  write(lun,'(7i10)') NZ, NGROUP, NELEM, NBIN, NGAS, nstep + 1, iskip

  do igroup = 1, NGROUP
    call CARMAGROUP_Get(carma, igroup, rc, r=r, rlow=rlow, rup=rup, dr=dr, rmass=rmass(:,igroup))
    if (rc < 0) stop "    *** FAILED CARMAGROUP_Get ***"

    do ibin = 1, NBIN
      write(lun,'(2i4,5e15.5)') igroup, ibin, r(ibin) * 1e4_f, rmass(ibin,igroup), dr(ibin) * 1e4_f, rlow(ibin) * 1e4_f, rup(ibin) * 1e4_f
    
      if (idiag .eq. 1) then
        write(lundiagn,'(i10,i12,5e15.5)') igroup, ibin, r(ibin) * 1e4_f, &
	      rmass(ibin,igroup), dr(ibin) * 1e4_f, rlow(ibin) * 1e4_f, rup(ibin) * 1e4_f
      end if

    end do

  end do

  
  if (idiag .eq. 1) then  
    write(lundiagn,*) 'ELEMENT INFORMATION:'
    write(lundiagn,'(2A10,A18)') 'IELEM', 'IGROUP', 'NAME' 

    do ielem = 1, NELEM
      call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup, name=name)
      if (rc < 0) stop "    *** FAILED ***"
      write(lundiagn,'(i8,i10,A35)') ielem, igroup, name
    end do

    write(lundiagn,*) 'GAS INFORMATION:'
    write(lundiagn,'(A6,A15,A35)') 'IGAS', 'NAME', 'WTMOL (g/mol)'

    do igas = 1, NGAS
      call CARMAGAS_Get(carma, igas, rc, name=gname, wtmol=wtmol)
      if (rc < 0) stop "    *** FAILED ***"
      write(lundiagn,'(i4,A40,e10.3)') igas, gname, wtmol
    end do

    write(lundiagn,*) 'ATMOSPHERE INFORMATION:'
    write(lundiagn,'(A5,4A15)') 'Z', 'ALTITUDE', 'dALT', 'PRESSURE', 'TEMPERATURE'
    write(lundiagn,'(A5,4A15)') '', '(km)', '(km)', '(mbars)', '(K)'  
  end if



  do i = 1, NZ
    write(lun,'(i3,5e15.5)') i, zc(i), zl(i+1)-zl(i), p(i) * 10._f, t(i), ekz(i)
  
    if (idiag .eq. 1) then
      write(lundiagn,'(i5,4e15.5)') i, zc(i) / 1000._f, (zl(i+1)-zl(i)) / 1000._f, p(i) / 100._f, t(i)
    end if

  end do

  if (idiag .eq. 1) then
    write(lundiagn,*) ' '
  end if

  write(*,*) ""

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define production rates of particles and gases

  prodrate(:,:,:) = 0._f
  prodrate_gas(:,:) = 0._f


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define constant wind speeds (temporally varying winds defined within step loop below)

  rho_atm_cgs = p(:) / (RGAS/wtmol_air(:) * t(:))


  ! Initialize longitudinal steps
  
  closeto_temp_profile = int(0)
  current_distance = 0
  rotation_counter = 0
  current_step = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Boundary Conditions

  ftopp(:,:) = 0._f
  fbotp(:,:) = 0._f
  pctop(:,:) = 0._f
  pcbot(:,:) = 0._f
  gctop(:) = 1e-50_f
  gcbot(:) = 1e-50_f
  ftopg(:) = 0._f
  fbotg(:) = 0._f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do ielem = 1, NELEM-1
    call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup, name=name)
    if (rc < 0) stop "    *** FAILED CARMA_ELEMENT_Get ***"
    write(lun,'(A35)', advance="no") name
  end do

  call CARMAELEMENT_Get(carma, NELEM, rc, igroup=igroup, name=name)
  if (rc < 0) stop "    *** FAILED CARMA_ELEMENT_Get ***"
  write(lun,'(A35)') name


  write(lun,*) 0
  write(lunf,*) 0
  do j = 1, NBIN
   do i = 1, NZ
    write(lun, '(i3,i4)', advance="no") j, i

    do ielem = 1, NELEM
      write(lun, '(e11.3)', advance="no") real(mmr(i,ielem,j) * rho_atm_cgs(i) / rmass(j,elem2group(ielem)))
    end do

    do igas = 1, NGAS
      write(lun, '(2e11.3)', advance="no") &
        real(mmr_gas(i,igas) * 1.0e6_f / (wtmol_gas(igas) / wtmol_air(i))), &
        0._f
    end do

    write(lun, '(f8.0)') 0.0
      
    end do
  enddo

  binmultiple = int(NBIN / 10._f)

  if (irestart .eq. 1) then
    open(unit=lunres,file=fileprefix // filename_restart(1:len_trim(filename_restart)) // filesuffix_restart,form='unformatted',status="unknown")
    read(lunres) istep_old, xc, dx, yc, dy, &
         zc, zl, p, pl, t, rho_atm_cgs, &
         mmr, mmr_gas, mmr_gas_old, satliq, &
         satice, satliq_old, satice_old, svpliq, wtpct, &
         zsubsteps, r, rlow, rup, dr, rmass, pflux, gflux, &
	 prodrate, prodrate_mass, prodrate_gas, totmass, &
	 winds, ekz, ftopp, fbotp, pctop, pcbot, gctop, &
	 gcbot, ftopg, fbotg, gasprod, rhompe, growpe, &
         evappe, growlg, evaplg
     write(*,*)'read restart file'
     rewind(lunres)
     !     istep = istep + 1
    close(unit=lunres)
  endif

  call CARMASTATE_CreateFromReference(cstate, carma_ptr, time, dtime, NZ, &
                         	      ! igridv, I_CART, lat, lon, &
                                 I_CART, I_CART, lat, lon, &

                         	      xc(:), dx(:), &
                         	      yc(:), dy(:), &
                         	      zc(:), zl(:), p(:), &
                         	      pl(:), t(:), wtmol_air(:), grav(:), rplanet, rc, winds=winds(:), ekz=ekz(:), met=met, t0 = t0_in)
  if (rc < 0) stop "    *** FAILED CARMASTATE_CreateFromReference ***"

  ! Iterate the model over a few time steps.
  do istep = 1, nstep

    !open(unit=lunres,file="venus_atmosphere_300z_45bins_noinit_lf_h2ofixed_10nmCNs_t0030et_2.dat",form='unformatted',status="unknown")
    open(unit=lunres,file=fileprefix // filename_restart(1:len_trim(filename_restart)) // filesuffix_restart,form='unformatted',status="unknown")
    open(unit=lunp,file = fileprefix // temp // filename(1:len_trim(filename)) // filesuffix, status="unknown", position=file_pos)
    open(unit=lunfp,file = fileprefix // temp // flux // filename(1:len_trim(filename)) // filesuffix, status="unknown", position=file_pos)
    open(unit=lunratesp,file = fileprefix // temp // rates // filename(1:len_trim(filename)) // filesuffix, status="unknown", position=file_pos)

    ! Calculate the model time.
    time = (istep - 1) * dtime

    write(*,*) 'istep, time', istep, time, filename(1:len_trim(filename))

    if (idiag .eq. 1) then
      do iz = 1, NZ
        totmass(iz) = sum((mmr(iz,1,:)+mmr(iz,2,:)) * abs((pl(iz+1) - pl(iz))) / 88.7_f) / deltaz / 100._f
      end do
      write(lundiagn,'(A6,I10,A12,f15.2,A8)') 'STEP:', istep, 'TIME:', istep*dtime, 'SECONDS'
      write(lundiagn,*) ' '		
      write(lundiagn,*) 'PART 1: TOTAL MASS AT START OF TIME STEP'
      write(lundiagn,*) '****************************************'						
      write(lundiagn,*) 'Z: ALTITUDE LEVEL INDEX'								
      write(lundiagn,*) 'PARTMASS: TOTAL MASS DENSITY OF GAS AT Z [g/cm3]'
      write(lundiagn,*) 'GASMASS: TOTAL MASS DENSITY OF PARTICLES AT Z [g/cm3]'
      write(lundiagn,*) '****************************************'						
      write(lundiagn,'(A6,2A15)') 'Z', 'PARTMASS', 'GASMASS'
      do iz = 1, NZ
	      write(lundiagn,'(i6,2e15.3)') iz, totmass(iz), (mmr_gas(iz,1)+mmr_gas(iz,2)) * abs((pl(iz+1) - pl(iz))) / 88.7_f / deltaz / 100._f
      end do
      write(lundiagn,*) ' '
      write(lundiagn,'(A6,2e15.3)') 'TOT:',sum(totmass), &
                                    sum((mmr_gas(:,1)+mmr_gas(:,2)) * abs(pl(2:NZP1) - pl(1:NZ)) / 88.7_f / deltaz / 100._f ) 
      write(lundiagn,*) ' '
      write(lundiagn,*) 'COLUMN TOTAL MASS = TOT * NZ * deltaz (in meters) * (100 cm / m)'
      startcd = sum(totmass) * NZ * deltaz * 100._f + sum((mmr_gas(:,1)+mmr_gas(:,2)) * abs(pl(2:NZP1) - pl(1:NZ)) / 88.7_f ) * NZ
      write(lundiagn,'(A28,e28.15)') 'COLUMN TOTAL MASS [g/cm2]: ', startcd
      write(lundiagn,*) ' '
    end if


    if (IS_2D .eq. 1) then
      current_distance = (velocity_avg*time)/distance_btwn_elements !in grid space
      rotation_counter = int(current_distance/NLONGITUDE)
      current_step = current_distance - (NLONGITUDE*rotation_counter)
      closeto_temp_profile = int(current_step)!int(current_step)+1

      write(*,*) rotation_counter, closeto_temp_profile, current_step

      if (closeto_temp_profile .eq. NLONGITUDE) then
        t(:) = temp_equator(:,1)
      else
        t(:) = temp_equator(:,closeto_temp_profile+1)
      end if

     end if


    ! To do: change gas input rate; add gaussian distribution to size of CNs being added; change nucleation rate with
    ! substepping; look closer at eddy diffusion machinery, share with Bardeen; read Bardeen's emails more!

    ! Create a CARMASTATE for this column.
    call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                          ! igridv, I_CART, lat, lon, &

                           I_CART, I_CART, lat, lon, &
                           xc(:), dx(:), &
                           yc(:), dy(:), &
                           zc(:), zl(:), p(:), &
                           pl(:), t(:), wtmol_air(:), grav(:), rplanet, rc, told=t(:), winds=winds(:), ekz=ekz(:), &
			   ftopp=ftopp,fbotp=fbotp,pctop=pctop,pcbot=pcbot, &
			   gctop=gctop,gcbot=gcbot,ftopg=ftopg,fbotg=fbotg,met=met)
    if (rc < 0) stop "    *** FAILED CARMASTATE_Create ***"


    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc)
        if (rc < 0) stop "    *** FAILED CARMASTATE_SetBin ***"
      end do
    end do

    do igas = 1, NGAS
      call CARMASTATE_SetGas(cstate, igas, mmr_gas(:,igas), rc, mmr_old=mmr_gas_old(:,igas), &
                             satice_old=satice(:,igas), satliq_old=satliq(:,igas))
      if (rc < 0) stop "    *** FAILED CARMASTATE_SetGas ***"
    end do

    ! Execute the step
    call CARMASTATE_Step(cstate, rc)
    if (rc < 0) stop "    *** FAILED CARMASTATE_Step ***"

    ! Get the retry stats and the updated temperature.
    call CARMASTATE_Get(cstate, rc, nsubstep=nsubsteps, nretry=nretries, zsubsteps=zsubsteps)
    if (rc < 0) stop "    ***  FAILED CARMASTATE_Get ***"

    ! Get the updated bin mmr.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,ielem,ibin), rc, numberDensity=numden(:,ielem,ibin), pflux=pflux(:,ibin,ielem))
        if (rc < 0) stop "    *** FAILED CARMASTATE_GetBin ***"
      end do
    end do

    do igas = 1, NGAS
      call CARMASTATE_GetGas(cstate, igas, mmr_gas(:,igas), rc, satliq=satliq(:,igas), &
                             satice=satice(:,igas), eqliq = svpliq(:,igas), wtpct=wtpct(:,igas), gflux=gflux(:,igas))
      if (rc < 0) stop "    *** FAILED CARMASTATE_GetGas ***"
    end do

    call CARMASTATE_GetDiag(cstate,rc,rhompe_tot=rhompe,growpe_tot=growpe,evappe_tot=evappe, &
                            growlg_tot=growlg,evaplg_tot=evaplg,gasprod_tot=gasprod)
    if (rc < 0) stop "    *** FAILED CARMASTATE_GetDiag ***"

    mmr_gas_old = mmr_gas

    ! Get the updated temperature.
!    call CARMASTATE_GetState(cstate, rc, t=t(:))
!    if (rc < 0) stop "    *** FAILED ***"

    ! Write output

 ! write(*,*) rmass


     if ((istep .eq. 1).and.(irestart .eq. 0)) then !PARAM
      open(unit=gas_in, file=gas_input_file, status='old')
      read(gas_in, *)

      do i = 1, NZ
        READ(gas_in, *) mmr_gas(i, :)
        do j = 1, NGAS
          mmr_gas(i,j) = min(mmr_gas(i,j) * (wtmol_gas(j)/ wtmol_air(1)) * 10._f**met,svpliq(i,j) * (wtmol_gas(j) / wtmol_air(1)))
        end do
      end do

      close(unit=gas_in)

      do i=1, NGAS
        gcbot(i) = mmr_gas(1, i) + 1e-50_f
      end do

     endif


    if (MOD (istep, iskip) .eq. 0) then

      write(*,*) 'Recorded'
      if (IS_2D .eq. 1) then
        write(lun,*) (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
        write(lunp,*) (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
        write(lunf,*) (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
        write(lunfp,*) (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
        write(lunrates,*) (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
        write(lunratesp,*) (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
      else
        write(lun,*) (istep)*dtime
        write(lunp,*) (istep)*dtime
        write(lunf,*) (istep)*dtime
        write(lunfp,*) (istep)*dtime
        write(lunrates,*) (istep)*dtime
        write(lunratesp,*) (istep)*dtime
      end if

      do j = 1, NBIN
        do i = 1, NZ

          write(lun, '(i3,i4)', advance="no") j, i
          write(lunp, '(i3,i4)', advance="no") j, i

          do ielem = 1, NELEM
            write(lun, '(e11.3)', advance="no") real(numden(i, ielem, j))
            write(lunp, '(e11.3)', advance="no") real(numden(i, ielem, j))
          end do

          do igas = 1, NGAS
            write(lun, '(2e11.3)', advance="no") &
            real(mmr_gas(i,igas) * 1.0e6_f / (wtmol_gas(igas) / wtmol_air(i))), &
            real(svpliq(i,igas) * 1.0e6_f)
            write(lunp, '(2e11.3)', advance="no") &
            real(mmr_gas(i,igas) * 1.0e6_f / (wtmol_gas(igas) / wtmol_air(i))), &
            real(svpliq(i,igas) * 1.0e6_f)
          end do

          write(lun, '(f8.0)') zsubsteps(i)
          write(lunp, '(f8.0)') zsubsteps(i)

        end do
      end do

    !endif

      if (idiag .eq. 1) then
        do iz = 1, NZ
          totmass(iz) = sum((mmr(iz,1,:)+mmr(iz,2,:)) * abs((pl(iz+1) - pl(iz))) / 88.7_f) / deltaz / 100._f
        end do		

        write(lundiagn,*) 'PART 6: TOTAL MASS AT END OF TIME STEP'
        write(lundiagn,*) '****************************************'						
        write(lundiagn,*) 'Z: ALTITUDE LEVEL INDEX'								
        write(lundiagn,*) 'PARTMASS: TOTAL MASS DENSITY OF GAS AT Z [g/cm3]'
        write(lundiagn,*) 'GASMASS: TOTAL MASS DENSITY OF PARTICLES AT Z [g/cm3]'
        write(lundiagn,*) '****************************************'						
        write(lundiagn,'(A6,2A15)') 'Z', 'PARTMASS', 'GASMASS'

        do iz = 1, NZ
          write(lundiagn,'(i6,2e15.3)') iz, totmass(iz), (mmr_gas(iz,1)+mmr_gas(iz,2)) * abs((pl(iz+1) - pl(iz))) / 88.7_f / deltaz / 100._f
        end do

        write(lundiagn,*) ' '
        write(lundiagn,'(A6,2e15.3)') 'TOT:',sum(totmass), &
                                      sum((mmr_gas(:,1)+mmr_gas(:,2)) * abs(pl(2:NZP1) - pl(1:NZ)) / 88.7_f / deltaz / 100._f ) 
        write(lundiagn,*) ' '
        write(lundiagn,*) 'COLUMN TOTAL MASS = TOT * NZ * deltaz (in meters) * (100 cm / m)'
          endcd = sum(totmass) * NZ * deltaz * 100._f + sum((mmr_gas(:,1)+mmr_gas(:,2)) * abs(pl(2:NZP1) - pl(1:NZ)) / 88.7_f ) * NZ
        write(lundiagn,'(A28,e28.15)') 'COLUMN TOTAL MASS [g/cm2]: ', endcd
        write(lundiagn,*) ' '
        write(lundiagn,*) 'PART 7: MASS CONSERVATION SUMMARY'
        write(lundiagn,*) '******************************************************************************'
        write(lundiagn,*) ''    
        write(lundiagn,'(A62)') 'MASS CONSERVATION: Ab - Aa = (B + C + D) * E'
        write(lundiagn,*) ''
        write(lundiagn,'(A45,e28.15)') 'Aa. TOTAL COLUMN DENSTY AT START (g/cm2):', startcd
        write(lundiagn,'(A45,e28.15)') 'Ab. TOTAL COLUMN DENSTY AT END (g/cm2):', endcd
        write(lundiagn,'(A45,e28.15)') 'B. TOTAL INPUT RATE (g/cm2/s):', inputrate
        write(lundiagn,'(A45,e28.15)') 'C. TOTAL VERTICAL PARTICLE FLUX (g/cm2/s):', vertpartflux
        write(lundiagn,'(A45,e28.15)') 'D. TOTAL VERTICAL GAS FLUX (g/cm2/s):', vertgasflux
        write(lundiagn,'(A45,e28.15)') 'E. TIME STEP (s):', dtime
        write(lundiagn,*) ''
        write(lundiagn,'(A45,e28.15)') 'Ab - Aa = ', endcd - startcd
        write(lundiagn,'(A45,e28.15)') '(B + C + D) * E = ', (inputrate + vertgasflux + vertpartflux) * dtime
        write(lundiagn,'(A45,e28.15)') 'Ab - Aa -[(B + C + D) * E] = ', endcd - startcd - &
          [(inputrate + vertpartflux + vertgasflux) * dtime]
        write(lundiagn,*) ''    
        write(lundiagn,*) '************************************************************10000000******************'
        write(lundiagn,*) ''    
        write(lundiagn,*) ''  
      end if





    !if (istep .ge. 89000) then
       write(lunres) istep+1, xc, dx, yc, dy, &
         zc, zl, p, pl, t, rho_atm_cgs, &
         mmr, mmr_gas, mmr_gas_old, satliq, &
         satice, satliq_old, satice_old, svpliq, wtpct, &
         zsubsteps, r, rlow, rup, dr, rmass, pflux, gflux, &
	 prodrate, prodrate_mass, prodrate_gas, totmass, &
	 winds, ekz, ftopp, fbotp, pctop, pcbot, gctop, &
	 gcbot, ftopg, fbotg, gasprod, rhompe, growpe, &
         evappe, growlg, evaplg
       write(*,*)'write restart file'
       rewind(lunres)
    !endif
    endif

    close(unit=lunp)
    close(unit=lunfp)
    close(unit=lunratesp)
    close(unit=lunres)

  end do   ! time loop

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc < 0) stop "    *** FAILED CARMASTATE_Destroy ***"

  ! Close the output file
  close(unit=lun)
  close(unit=lunf)
  close(unit=lunrates)
  if (idiag .eq. 1) then
    close(unit=lundiagn)
  end if



  if (rc < 0) stop "    *** FAILED ***"

  write(*,*)  ""
  write(*,*) "CARMA_Destroy() ..."
  call CARMA_Destroy(carma, rc)
  if (rc < 0) stop "    *** FAILED CARMA_Destroy ***"

end subroutine
