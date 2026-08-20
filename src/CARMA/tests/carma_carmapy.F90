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
  use carma_rce

  implicit none

  integer, parameter    :: io = 90

  ! Griz sizes
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
  integer   :: t_evolves
  integer   :: igridv
  integer   :: idocoag
  logical   :: do_coag

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

  ! Boundary Conditions
  real(kind=f), allocatable   :: ftopp(:,:)
  real(kind=f), allocatable   :: fbotp(:,:)
  real(kind=f), allocatable   :: pctop(:,:)
  real(kind=f), allocatable   :: pcbot(:,:)
  real(kind=f), allocatable   :: gctop(:)
  real(kind=f), allocatable   :: gcbot(:)
  real(kind=f), allocatable   :: ftopg(:)
  real(kind=f), allocatable   :: fbotg(:)

  ! Microphysical rates
  real(kind=f), allocatable   :: gasprod(:,:)
  real(kind=f), allocatable   :: rhompe(:,:,:)
  real(kind=f), allocatable   :: rnucpe(:,:,:)
  real(kind=f), allocatable   :: rnuclg(:,:,:)
  real(kind=f), allocatable   :: growpe(:,:,:)
  real(kind=f), allocatable   :: growlg(:,:,:)
  real(kind=f), allocatable   :: evappe(:,:,:)
  real(kind=f), allocatable   :: evaplg(:,:,:)


  real(kind=f), allocatable   :: corefrac(:,:,:)
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
  integer, parameter    :: lunrad = 57

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
  character(30)      	:: hill_formula
  character(30)  :: type_spec
  character(30)      	:: gname
  character(len=3)	:: fileprefix
  character(len=5)	:: temp = 'temp_'
  character(len=5)	:: flux = 'flux_'
  character(len=5)	:: diag = 'diag_'
  character(len=6)	:: rates = 'rates_'
  character(len=12)	:: radpre = 'temperature_'
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
  character(len=100)  :: winds_file
  character(len=100)  :: optics_file

  character(len=100)  ::  g_boundary_file, p_boundary_file


  character(len=20)   :: file_pos

  real(kind=f)  :: rmu_1, rmu_2, rmu_3, rmu_4, thcond_0, thcond_1, thcond_2, CP
  real(kind=f)  :: distance_btwn_elements, circumference, rotation_counter, slope, intercept
  real(kind=f)  :: current_distance, num_steps_btwn, current_step, RPLANET_DAT, restart_distance
  integer       :: closeto_temp_profile
  integer       :: ibbnd_pc, itbnd_pc, ibbnd_gc, itbnd_gc

  namelist / io_files / filename, filename_restart, fileprefix, gas_input_file,&
         centers_file, levels_file,temps_file, groups_file, elements_file, &
         gases_file, growth_file, nuc_file, coag_file, winds_file, &
         g_boundary_file, p_boundary_file, optics_file
  
  namelist / physical_params / wtmol_air_set, grav_set, rplanet, velocity_avg,&
         met, rmu_1, rmu_2, rmu_3, rmu_4, thcond_0, thcond_1, thcond_2, CP

  ! Radiative-convective coupling. Absent from input.nml unless the run was
  ! configured with Carma.set_radiation(), and inert when do_radiation is 0.
  type(rce_type)      :: rce
  integer             :: do_radiation, rad_mode, rad_gap_max, ios_rad
  ! Set in Python; the clouds' opacity is zeroed in optics_file rather than
  ! branched on here, so this only records which case the run was.
  integer             :: cloud_rad
  character(len=256)  :: ck_table_file
  real(kind=f)        :: t_int, rad_accel, rad_dT_max
  real(kind=f)        :: t_irr, t_star, rad_mu0, rad_w_surf
  real(kind=f)        :: rad_dT_tol, rad_dtau_tol
  ! Which adiabat the convective interior follows: 0 = Parmentier fit,
  ! 1 = the tabulated H/He gradient read from adiabat_file.
  integer             :: adiabat
  character(len=256)  :: adiabat_file
  ! Whether the eddy diffusion profile follows the evolving temperature
  ! profile: 0 = hold the levels_file values, 1 = mixing length theory.
  integer             :: kzz_mode
  ! Multiplier on the mixing length. 1 is the parameterisation as PICASO
  ! has it; Kzz goes as L^(4/3), so this moves it by its 4/3 power.
  real(kind=f)        :: kzz_mixl_scale
  integer             :: nsolve_last   !! rce%nsolve as of the previous output

  namelist / radiation / do_radiation, ck_table_file, t_int, &
         t_irr, t_star, rad_mu0, rad_w_surf, &
         rad_mode, rad_accel, rad_dT_max, rad_dT_tol, rad_dtau_tol, &
         rad_gap_max, adiabat, adiabat_file, kzz_mode, kzz_mixl_scale, cloud_rad

  namelist / input_params / NZ, NELEM, NGROUP, NGAS, NBIN, NSOLUTE, NWAVE, &
         NLONGITUDE, irestart, idiag, iskip, nstep, dtime, NGROWTH, NNUC, &
         NCOAG, IS_2D, t_evolves, igridv, iappend, idocoag, itbnd_pc, ibbnd_pc, itbnd_gc, &
         ibbnd_gc

  real(kind=f) ::rho_cond, surften_0, coldia, vp_offset, vp_tcoeff, surften_slope, vp_metcoeff, vp_logpcoeff, lat_heat_e, desorption
  integer :: is_type3, stofact
  real(kind=f) :: wtmol_core
  integer :: igrp

  real(kind=f), allocatable ::tempr(:), pre(:), prel(:), alt(:), altl(:), wtmol_air(:), grav(:), ekz(:), ekzl(:), wtmol_gas(:)
  real(kind=f), allocatable ::temp_equator(:, :), p_equator_center(:), p_equator_level(:), velocity(:), longitudes(:)
  integer, allocatable :: elem2group(:)

  ! Condensate properties are written per-group by CARMApy (new API), but the
  ! engine still owns them per-gas. These arrays buffer the group properties so
  ! they can be translated onto the corresponding gas (the gas grown by each
  ! group's "Volatile" element) when the gases are created below.
  real(kind=f), allocatable :: grp_wtmol(:), grp_rho_cond(:), grp_surften_0(:), grp_coldia(:)
  real(kind=f), allocatable :: grp_vp_offset(:), grp_vp_tcoeff(:), grp_surften_slope(:)
  real(kind=f), allocatable :: grp_vp_metcoeff(:), grp_vp_logpcoeff(:), grp_lat_heat_e(:)
  real(kind=f), allocatable :: grp_desorption(:), grp_wtmol_core(:)
  integer, allocatable :: grp_iroutine(:), grp_is_type3(:), grp_stofact(:)
  integer, allocatable :: gas2group(:)

  ! Per-group cloud optics, read from optics_file when NWAVE > 0. Generated in
  ! Python on exactly this bin grid and wavelength set, so they are handed to
  ! CARMAGROUP_Create as-is; see carmapy.radiation.write_optics_file.
  real(kind=f), allocatable :: grp_qext(:,:,:), grp_ssa(:,:,:), grp_asym(:,:,:)
  integer :: opt_ngroup, opt_nwave, opt_nbin, opt_ig, opt_iw, opt_ib

  logical, allocatable :: elem_is_number(:)
  real(kind=f), allocatable :: grp_r(:,:)   !! bin radii per group [cm]
  integer :: opt_skip1, opt_skip2, opt_skip3


  fileprefix = trim(fileprefix)

  if (idocoag .eq. 0) then
    DO_COAG = .False.
  else
    DO_COAG = .True.
  endif


  write(*,*) ""

  ! Defaults for namelist params that may be missing from older input.nml files.
  t_evolves = 0   ! 0 = T,P fixed across the run -> fixed initialization, setupgkern caches its outputs
  optics_file = ""  ! only read when NWAVE > 0, i.e. radiatively coupled runs

  do_radiation  = 0        ! 0 = fixed T,P -- the whole radiation group is inert
  ck_table_file = ""
  t_int         = 0._f
  t_irr         = 0._f      ! 0 = isolated object, no incident beam
  t_star        = 0._f      ! 0 = shape the beam like t_irr itself
  rad_mu0       = 0.5_f     ! dayside average
  rad_w_surf    = 0._f      ! gas giant: no surface to reflect
  rad_mode      = I_RCE_EQUILIBRIUM
  rad_accel     = 1._f
  rad_dT_max    = 0.5_f
  rad_dT_tol    = 1._f
  rad_dtau_tol  = 0.02_f
  rad_gap_max   = 100
  adiabat       = 0
  adiabat_file  = ""
  kzz_mode      = I_KZZ_STATIC
  kzz_mixl_scale = 1._f
  cloud_rad     = 1

  open(unit=10, file=nml_file, status='old')
    read(10, nml=input_params)
    read(10, nml=io_files)
    read(10, nml=physical_params)
    ! Optional: older input.nml files have no radiation group at all, which
    ! leaves the defaults above and the run behaves exactly as before.
    read(10, nml=radiation, iostat=ios_rad)
  close(10)

  NZP1 = NZ + 1

  file_pos = "asis"
  if ((irestart .eq. 1).and.(iappend .eq. 1)) file_pos = "append"



  allocate(tempr(NZ), pre(NZ), prel(NZP1), alt(NZ), altl(NZP1), wtmol_air(NZ), grav(NZ), ekz(NZP1), ekzl(NZP1), wtmol_gas(NGAS))
  allocate(temp_equator(NZ, NLONGITUDE), p_equator_center(NZ), p_equator_level(NZP1), velocity(NLONGITUDE), longitudes(NLONGITUDE))
  allocate(elem2group(NELEM), elem_is_number(NELEM), grp_r(NBIN, NGROUP))
  allocate(grp_wtmol(NGROUP), grp_rho_cond(NGROUP), grp_surften_0(NGROUP), grp_coldia(NGROUP))
  allocate(grp_vp_offset(NGROUP), grp_vp_tcoeff(NGROUP), grp_surften_slope(NGROUP))
  allocate(grp_vp_metcoeff(NGROUP), grp_vp_logpcoeff(NGROUP), grp_lat_heat_e(NGROUP))
  allocate(grp_desorption(NGROUP), grp_wtmol_core(NGROUP))
  allocate(grp_iroutine(NGROUP), grp_is_type3(NGROUP), grp_stofact(NGROUP))
  allocate(gas2group(NGAS))
  gas2group = 0
  allocate(winds(NZ))

  ! Cloud optics for the radiatively coupled path. NWAVE is 0 for an ordinary
  ! run, in which case nothing is read and the group optics stay unallocated.
  if (NWAVE > 0) then
    allocate(grp_qext(NWAVE, NBIN, NGROUP), grp_ssa(NWAVE, NBIN, NGROUP), &
             grp_asym(NWAVE, NBIN, NGROUP))

    if (len_trim(optics_file) == 0) then
      write(*,*) "*** NWAVE > 0 but no optics_file was given in input.nml ***"
      stop
    end if

    open(10, file = optics_file, status='old')
    read(10, *) opt_ngroup, opt_nwave, opt_nbin

    ! A silent shape mismatch here would put the wrong wavelength's optics on
    ! every bin, so it is checked rather than trusted.
    if (opt_ngroup /= NGROUP .or. opt_nwave /= NWAVE .or. opt_nbin /= NBIN) then
      write(*,*) "*** optics_file is (ngroup,nwave,nbin) = ", opt_ngroup, opt_nwave, opt_nbin, &
                 " but input.nml says ", NGROUP, NWAVE, NBIN, " ***"
      stop
    end if

    read(10, *)   ! column header
    do opt_ig = 1, NGROUP
      do opt_iw = 1, NWAVE
        do opt_ib = 1, NBIN
          read(10, *) opt_skip1, opt_skip2, opt_skip3, &
                      grp_qext(opt_iw, opt_ib, opt_ig), &
                      grp_ssa(opt_iw, opt_ib, opt_ig), &
                      grp_asym(opt_iw, opt_ib, opt_ig)
        end do
      end do
    end do
    close(10)

    write(*,*) "  Read cloud optics: ", NWAVE, " bands x ", NBIN, " bins x ", NGROUP, " groups"
  end if

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


  open(12, file = winds_file) 

  do i=1, NZ
    read(12,*) winds(i)
  end do

  close(12)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Boundary Conditions

  allocate(ftopp(NBIN,NELEM))
  allocate(fbotp(NBIN,NELEM))
  allocate(pctop(NBIN,NELEM))
  allocate(pcbot(NBIN,NELEM))
  allocate(gctop(NGAS))
  allocate(gcbot(NGAS))
  allocate(ftopg(NGAS))
  allocate(fbotg(NGAS))

  ftopp(:,:) = 0._f
  fbotp(:,:) = 0._f
  pctop(:,:) = 0._f
  pcbot(:,:) = 0._f
  gctop(:) = 1e-50_f
  gcbot(:) = 1e-50_f
  ftopg(:) = 0._f
  fbotg(:) = 0._f

  open(12, file = g_boundary_file) 
  read(12, *)
  do i=1, NGAS
    read(12, *) name, gctop(i), gcbot(i), ftopg(i), fbotg(i)
  enddo
  close(12)

  open(12, file = p_boundary_file) 
  read(12, *)
  do i=1, NBIN 
    do j=1, NELEM
      read(12, *) name, ibin, pctop(i,j), pcbot(i,j), ftopp(i,j), fbotp(i,j)
    enddo
  enddo
  close(12)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  wtmol_air(:) =wtmol_air_set 
  grav(:) = grav_set 

  t0_in = MAXVAL(tempr)

  circumference = 2 *PI * rplanet !cm
  distance_btwn_elements = circumference/NLONGITUDE
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Open the output text file
  open(unit=lun,file =  filename(1:len_trim(filename)) // filesuffix, status="unknown", position=file_pos)
  open(unit=lunf,file =  flux // filename(1:len_trim(filename)) // filesuffix, status="unknown", position=file_pos)
  open(unit=lunrates,file =  rates // filename(1:len_trim(filename)) // filesuffix, status="unknown", position=file_pos)

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


  ! Microphysical Rates
  allocate(gasprod(NZ,NGAS))
  allocate(rhompe(NZ,NBIN,NELEM))
  allocate(rnucpe(NZ,NBIN,NELEM))
  allocate(rnuclg(NZ,NBIN,NGROUP))
  allocate(growpe(NZ,NBIN,NELEM))
  allocate(growlg(NZ,NBIN,NGROUP))
  allocate(evappe(NZ,NBIN,NELEM))
  allocate(evaplg(NZ,NBIN,NGROUP))

  allocate(corefrac(NZ,NBIN,NGROUP))


  ! Define the particle-grid extent of the CARMA test
  write(*,*) "Create CARMA Object ..."

  if (idiag .eq. 1) then 
    open(lundiagn)
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
    ! New CARMApy group file carries the condensate properties. They are
    ! buffered here and applied to the corresponding gas (see gas creation).
    read(10, *) name, rmin, grp_wtmol(i), grp_iroutine(i), grp_rho_cond(i), &
     grp_surften_0(i), grp_coldia(i), grp_vp_offset(i), grp_vp_tcoeff(i), &
     grp_is_type3(i), grp_surften_slope(i), grp_vp_metcoeff(i), grp_vp_logpcoeff(i), &
     grp_lat_heat_e(i), grp_desorption(i), grp_stofact(i), grp_wtmol_core(i)

    write(*,*) "Add " //TRIM(name)//"..."

    if (NWAVE > 0) then
      call CARMAGROUP_Create(carma, i, name, rmin, 2._f, I_SPHERE, 1._f, &
      .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE., &
      qext=grp_qext(:,:,i), ssa=grp_ssa(:,:,i), asym=grp_asym(:,:,i))
    else
      call CARMAGROUP_Create(carma, i, name, rmin, 2._f, I_SPHERE, 1._f, &
      .FALSE., rc, do_vtran=.TRUE., is_sulfate=.FALSE.)
    end if

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

    ! WC TODO: Handle this in python
    if(trim(type_spec) == "Volatile") then
      itype = I_VOLATILE
    else if(trim(type_spec) == "Core Mass") then
      itype = I_COREMASS
    else
      stop "invalid element type"
    endif
      call CARMAELEMENT_Create(carma, i, igroup, name, rho, itype, icomposition, rc)
      elem2group(i) = igroup
      ! Core-mass elements hold a mass concentration in the same array, so only
      ! the others carry the group's particle number density.
      elem_is_number(i) = (itype /= I_COREMASS)
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

  ! Build the gas -> group map from the growth pathways: each growth line pairs a
  ! volatile element (ito) with the gas it exchanges with (igas); elem2group then
  ! gives the group that owns that gas's condensate properties. Core (Core Mass)
  ! elements have no growth pathway, so they never appear here. (This file is read
  ! again below for CARMA_AddGrowth.)
  open(10, file=growth_file)
  read(10, *)
  do i = 1, NGROWTH
    read(10, *) ito, igas
    gas2group(igas) = elem2group(ito)
  enddo
  close(10)

  open(10, file=gases_file)
  read(10, *)
  do i=1, NGAS

    ! New CARMApy gas file holds gas-phase properties only. The condensate
    ! properties are pulled from the group that condenses this gas. hill_formula
    ! is read but unused by the engine (it is a Python-side fastchem label).
    read(10, *) name, wtmol_dif, icomposition, hill_formula

    igrp = gas2group(i)

    write(*,*) "Add "//trim(name) //" ..."
    if (igrp > 0) then
      call CARMAGAS_Create(carma, i, name, grp_wtmol(igrp), grp_iroutine(igrp), icomposition, &
       grp_rho_cond(igrp), grp_surften_0(igrp), grp_coldia(igrp), grp_vp_offset(igrp), &
       grp_vp_tcoeff(igrp), rc, is_type3=grp_is_type3(igrp), surften_slope=grp_surften_slope(igrp), &
       vp_metcoeff=grp_vp_metcoeff(igrp), vp_logpcoeff=grp_vp_logpcoeff(igrp), wtmol_dif=wtmol_dif, &
       lat_heat_e=grp_lat_heat_e(igrp), stofact=grp_stofact(igrp), desorption=grp_desorption(igrp))
    else
      ! Background gas with no condensate (e.g. H2O when there is no water cloud).
      ! No group grows it, so its condensate properties are never used; supply
      ! harmless placeholders and pick the vapor-pressure routine from composition.
      if (icomposition == I_GCOMP_H2O) then
        iroutine = I_VAPRTN_H2O_MURPHY2005
      else
        iroutine = I_VAPRTN_USER
      end if
      call CARMAGAS_Create(carma, i, name, wtmol_dif, iroutine, icomposition, &
       1._f, 0._f, 1.0e-8_f, 0._f, 0._f, rc, wtmol_dif=wtmol_dif)
    end if

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

  ! do_fixedinit initializes the fall velocities, diffusion coefficients and
  ! growth kernels once from the reference atmosphere. That is only valid while
  ! T and P are fixed, so a T-evolving run pays for per-step initialization.
  call CARMA_Initialize(carma, rc, do_cnst_rlh =.FALSE., do_coag=DO_COAG, &
                        do_fixedinit=(t_evolves .eq. 0), do_grow=.TRUE., &
                        do_explised=.FALSE., do_substep=.TRUE., do_print_init=.TRUE., &
                        do_vdiff=.TRUE., do_vtran=.TRUE., maxsubsteps=10, maxretries=10, &
                        itbnd_pc=itbnd_pc, ibbnd_pc=ibbnd_pc, itbnd_gc=itbnd_gc, ibbnd_gc=ibbnd_gc, &
                        do_t_evolves=(t_evolves .ne. 0))
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

  write(lun,'(7i10)') NZ, NGROUP, NELEM, NBIN, NGAS, nstep + 1, iskip

  do igroup = 1, NGROUP
    call CARMAGROUP_Get(carma, igroup, rc, r=r, rlow=rlow, rup=rup, dr=dr, rmass=rmass(:,igroup))
    if (rc < 0) stop "    *** FAILED CARMAGROUP_Get ***"

    do ibin = 1, NBIN
      write(lun,'(2i4,5e15.5)') igroup, ibin, r(ibin) * 1e4_f, rmass(ibin,igroup), dr(ibin) * 1e4_f, rlow(ibin) * 1e4_f, rup(ibin) * 1e4_f
      grp_r(ibin, igroup) = r(ibin)
    end do
  end do

  



  do i = 1, NZ
    write(lun,'(i3,5e15.5)') i, zc(i), zl(i+1)-zl(i), p(i) * 10._f, t(i), ekz(i)
  end do


  write(*,*) ""

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define production rates of particles and gases

  prodrate(:,:,:) = 0._f
  prodrate_gas(:,:) = 0._f


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define constant wind speeds (temporally varying winds defined within step loop below)

  rho_atm_cgs = p(:) / (RGAS/wtmol_air(:) * t(:))

  ! Set up the radiative-convective coupling. The optics read above are the
  ! cloud half of it; the ck table is the gas half.
  if (do_radiation .eq. 1) then
    if (NWAVE .le. 0) then
      write(*,*) "*** do_radiation = 1 requires NWAVE > 0 and an optics_file ***"
      stop
    end if

    call rce_init(rce, NZ, NBIN, NGROUP, NWAVE, igridv, trim(ck_table_file), &
                  t_int, t_irr, t_star, rad_mu0, rad_w_surf, &
                  CP, grav_set, wtmol_air_set, &
                  adiabat, trim(adiabat_file), kzz_mode, kzz_mixl_scale, &
                  rad_mode, rad_accel, rad_dT_max, rad_dT_tol, rad_dtau_tol, &
                  rad_gap_max, rc)
    if (rc < 0) stop "    *** FAILED rce_init ***"

    ! The evolving column gets its own output file, so the existing outputs
    ! keep their format and a fixed-T run produces no extra file at all.
    open(unit=lunrad, file = radpre // filename(1:len_trim(filename)) // filesuffix, &
         status="unknown", position=file_pos)
    write(lunrad,'(3i10,3e25.15)') NZ, nstep / iskip, iskip, t_int, &
                                   t_irr, rad_mu0
    nsolve_last = 0

    write(*,*) "Radiative coupling on: T_int =", t_int, &
               " K; T_irr =", t_irr, " K; cloud opacity =", cloud_rad
  end if


  ! Initialize longitudinal steps
  
  closeto_temp_profile = int(0)
  current_distance = 0
  rotation_counter = 0
  current_step = 0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do ielem = 1, NELEM-1
    call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup, name=name)
    if (rc < 0) stop "    *** FAILED CARMA_ELEMENT_Get ***"
    write(lun,'(A35)', advance="no") name
  end do

  call CARMAELEMENT_Get(carma, NELEM, rc, igroup=igroup, name=name)
  if (rc < 0) stop "    *** FAILED CARMA_ELEMENT_Get ***"
  write(lun,'(A35)') name


  write(lun,'(i1)') 0
  write(lunf,'(i1)') 0
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
    open(unit=lunres,file= filename_restart(1:len_trim(filename_restart)) // filesuffix_restart,form='unformatted',status="unknown")
    read(lunres) istep_old, xc, dx, yc, dy, &
         zc, zl, p, pl, t, rho_atm_cgs, &
         mmr, mmr_gas, mmr_gas_old, satliq, &
         satice, satliq_old, satice_old, svpliq, wtpct, &
         zsubsteps, r, rlow, rup, dr, rmass, pflux, gflux, &
	 prodrate, prodrate_mass, prodrate_gas, totmass, &
	 winds, ekz, ftopp, fbotp, pctop, pcbot, gctop, &
	 gcbot, ftopg, fbotg, gasprod, rhompe, growpe, &
         evappe, growlg, evaplg, restart_distance
     write(*,*)'read restart file'
     rewind(lunres)
     !     istep = istep + 1
    close(unit=lunres)
  endif

  call CARMASTATE_CreateFromReference(cstate, carma_ptr, time, dtime, NZ, &
                         	      igridv, I_CART, lat, lon, &
                         	      xc(:), dx(:), &
                         	      yc(:), dy(:), &
                         	      zc(:), zl(:), p(:), &
                         	      pl(:), t(:), wtmol_air(:), grav(:), rplanet, &
                                rmu_1, rmu_2, rmu_3, rmu_4, thcond_0, thcond_1, thcond_2, CP, &
                         	      rc, winds=winds(:), ekz=ekz(:), met=met, t0 = t0_in)
  if (rc < 0) stop "    *** FAILED CARMASTATE_CreateFromReference ***"
  
  ! Iterate the model over a few time steps.
  do istep = 1, nstep

    open(unit=lunres,file= filename_restart(1:len_trim(filename_restart)) // filesuffix_restart,form='unformatted',status="unknown")
    open(unit=lunp,file =  temp // filename(1:len_trim(filename)) // filesuffix, status="unknown", position=file_pos)
    open(unit=lunfp,file =  temp // flux // filename(1:len_trim(filename)) // filesuffix, status="unknown", position=file_pos)
    open(unit=lunratesp,file =  temp // rates // filename(1:len_trim(filename)) // filesuffix, status="unknown", position=file_pos)

    ! Calculate the model time.
    time = (istep - 1) * dtime

    write(*,*) 'istep, time', istep, time, filename(1:len_trim(filename))


    if (IS_2D .eq. 1) then
      current_distance = (velocity_avg*time)/distance_btwn_elements + restart_distance!in grid space
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


    ! Evolve the temperature profile and the altitude grid before the state is
    ! built, so this step's microphysics sees the updated column. The cloud
    ! field is the one left by the previous step, which is what makes the
    ! coupling explicit.
    if (do_radiation .eq. 1) then
      call rce_update(rce, dtime, p(:), pl(:), t(:), zc(:), zl(:), &
                      NELEM, numden, elem2group, elem_is_number, &
                      grp_r, grp_qext, grp_ssa, grp_asym)

      ! Mix with the eddy diffusion the new profile implies rather than the
      ! one the run started from. ekz is handed to CARMASTATE_Create below, so
      ! this is the whole of the coupling back into the microphysics.
      if (kzz_mode .eq. I_KZZ_MIXING_LENGTH) then
        ekz(:) = rce%kzz(:)
      end if
    end if

    ! To do: change gas input rate; add gaussian distribution to size of CNs being added; change nucleation rate with
    ! substepping; look closer at eddy diffusion machinery, share with Bardeen; read Bardeen's emails more!

    ! Create a CARMASTATE for this column.
    call CARMASTATE_Create(&
              cstate, carma_ptr, time, dtime, NZ, &
              igridv, I_CART, lat, lon, &
              xc(:), dx(:), &
              yc(:), dy(:), &
              zc(:), zl(:), p(:), &
              pl(:), t(:), wtmol_air(:), grav(:), rplanet, &
              rmu_1, rmu_2, rmu_3, rmu_4, thcond_0, thcond_1, thcond_2, CP, &
              rc, told=t(:), winds=winds(:), ekz=ekz(:), &
              ftopp=ftopp,fbotp=fbotp,pctop=pctop,pcbot=pcbot, &
              gctop=gctop,gcbot=gcbot,ftopg=ftopg,fbotg=fbotg,met=met, t0=t0_in)
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
          mmr_gas(i,j) = min(mmr_gas(i,j) * (wtmol_gas(j)/ wtmol_air(1)), svpliq(i,j) * (wtmol_gas(j) / wtmol_air(1)))
        end do
      end do


      close(unit=gas_in)

      ! do i=1, NGAS
      !   gcbot(i) = mmr_gas(1, i) + 1e-50_f
      ! end do

     endif


    if (MOD (istep, iskip) .eq. 0) then

      call CARMASTATE_GetDiag(cstate, &
      rc, &
      rhompe_tot=rhompe, &
      rnucpeup_tot=rnucpe, &
      rnuclg_tot=rnuclg, &
      growpe_tot=growpe, &
      growlg_tot=growlg, &
      evappe_tot=evappe, &
      evaplg_tot=evaplg, &
      corefrac=corefrac, &
      gasprod_tot=gasprod)

      if (rc < 0) stop "    *** FAILED CARMASTATE_GetDiag ***"



      write(*,*) 'Recorded'
      if (IS_2D .eq. 1) then
        write(lun,'(5e25.5)') (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
        write(lunp,'(5e25.5)') (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
        write(lunf,'(5e25.5)') (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
        write(lunfp,'(5e25.5)') (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
        write(lunrates,'(5e25.5)') (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
        write(lunratesp,'(5e25.5)') (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
      else
        write(lun,'(f25.5)') (istep)*dtime !TODO -- change to just istep?
        write(lunp,'(f25.5)') (istep)*dtime
        write(lunf,'(f25.5)') (istep)*dtime
        write(lunfp,'(f25.5)') (istep)*dtime
        write(lunrates,'(f25.5)') (istep)*dtime
        write(lunratesp,'(f25.5)') (istep)*dtime
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
            write(lun, '(2e25.15)', advance="no") &
            real(mmr_gas(i,igas) * 1.0e6_f / (wtmol_gas(igas) / wtmol_air(i))), &
            real(svpliq(i,igas) * 1.0e6_f)
            write(lunp, '(2e11.3)', advance="no") &
            real(mmr_gas(i,igas) * 1.0e6_f / (wtmol_gas(igas) / wtmol_air(i))), &
            real(svpliq(i,igas) * 1.0e6_f)
          end do

          write(lun, '(f8.0)') zsubsteps(i)
          write(lunp, '(f8.0)') zsubsteps(i)

          do ielem = 1, NELEM
            write(lunrates, '(3i4,7e13.3e3)') j, &
                              i, &
                              ielem, &
                              rhompe(i, j, ielem), &
                              rnucpe(i, j, ielem), &
                              growpe(i, j, ielem), &
                              evappe(i, j, ielem)
          enddo
          do igroup = 1, NGROUP
            write(lunrates, '(3i4,7e13.3e3)') j, &
                              i, &
                              igroup, &
                              rnuclg(i, j, igroup), &
                              growlg(i, j, igroup), &
                              evaplg(i, j, igroup), &
                              corefrac(i, j, igroup)
          enddo




        end do
      end do

      ! The radiative column: layer temperatures, heating rates and which
      ! layers convection left neutrally stable, then the level net fluxes the
      ! heating was differenced from. F_net at the top is the emergent flux,
      ! which is what sigma*t_int^4 is checked against.
      if (do_radiation .eq. 1) then
        write(lunrad,'(7e25.15)') (istep)*dtime, rce%fnet(NZP1), &
                                  real(rce%nsolve - nsolve_last, f), &
                                  real(rce%nzone, f), rce%conv_resid, &
                                  rce%absorbed_sw, rce%reflected_sw
        nsolve_last = rce%nsolve

        do i = 1, NZ
          write(lunrad,'(i5,3e25.15,i3)') i, p(i) * 10._f, t(i), rce%dtdt(i), &
                                          merge(1, 0, rce%is_conv(i))
        end do

        do i = 1, NZP1
          write(lunrad,'(i5,3e25.15)') i, pl(i) * 10._f, rce%fnet(i), ekz(i)
        end do
      end if

    !endif

        write(lunres) istep+1, xc, dx, yc, dy, &
          zc, zl, p, pl, t, rho_atm_cgs, &
          mmr, mmr_gas, mmr_gas_old, satliq, &
          satice, satliq_old, satice_old, svpliq, wtpct, &
          zsubsteps, r, rlow, rup, dr, rmass, pflux, gflux, &
          prodrate, prodrate_mass, prodrate_gas, totmass, &
          winds, ekz, ftopp, fbotp, pctop, pcbot, gctop, &
          gcbot, ftopg, fbotg, gasprod, rhompe, growpe, &
          evappe, growlg, evaplg, current_distance
       write(*,*)'write restart file'
       rewind(lunres)
    !endif
    endif

    close(unit=lunp)
    close(unit=lunfp)
    close(unit=lunratesp)
    close(unit=lunres)

  end do   ! time loop

  if (do_radiation .eq. 1) then
    write(*,*) "Radiation: solves =", rce%nsolve, " of", nstep, &
               " steps; dT clamp fired", rce%nclamp, "times"
    close(unit=lunrad)
    call rce_destroy(rce)
  end if

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
