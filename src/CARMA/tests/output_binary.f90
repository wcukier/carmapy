! Binary output module for CARMApy.
!
! Writes all output arrays for one timestep to a single binary file using
! Fortran stream I/O (no record markers), which numpy can read directly.
! The file is replaced on every call; Python reads it after each "Recorded"
! sentinel and appends the data to an HDF5 file.
!
! Binary layout (all floats are float64, ints are default int = int32):
!   8 x int32  : NZ, NGROUP, NELEM, NBIN, NGAS, istep, IS_2D, NZP1
!   1 x float64: time  (= istep * dtime)
!   NBIN*NGROUP float64: r(NBIN,NGROUP)    [Fortran column-major]
!   NBIN*NGROUP float64: rmass(NBIN,NGROUP)
!   NZ*NELEM*NBIN float64: numden(NZ,NELEM,NBIN)
!   NZ*NGAS float64: mmr_gas(NZ,NGAS)
!   NZ*NGAS float64: svpliq(NZ,NGAS)
!   NZ float64: zsubsteps(NZ)
!   NZ*NBIN*NELEM float64 x4: rhompe, rnucpe, growpe, evappe
!   NZ*NBIN*NGROUP float64 x4: rnuclg, growlg, evaplg, corefrac
!   4 x float64: current_distance, rotation_counter, current_step, lon_frac
!                (all zero for 1-D runs)
!   NZP1*NBIN*NELEM float64: pflux(NZP1,NBIN,NELEM)
!   NZP1*NGAS float64: gflux(NZP1,NGAS)

module output_binary_mod
  use carma_precision_mod
  implicit none

  integer, parameter, private :: lun_bin = 97

contains

  subroutine write_binary_output(bin_filename, &
      NZ, NGROUP, NELEM, NBIN, NGAS, NLONGITUDE, istep, IS_2D, dtime, &
      current_distance, rotation_counter, current_step, &
      r, rmass, numden, mmr_gas, svpliq, zsubsteps, &
      rhompe, rnucpe, growpe, evappe, rnuclg, growlg, evaplg, corefrac, &
      pflux, gflux)

    character(len=*), intent(in) :: bin_filename
    integer,      intent(in) :: NZ, NGROUP, NELEM, NBIN, NGAS, NLONGITUDE
    integer,      intent(in) :: istep, IS_2D
    real(kind=f), intent(in) :: dtime
    real(kind=f), intent(in) :: current_distance, rotation_counter, current_step
    real(kind=f), intent(in) :: r(NBIN, NGROUP)
    real(kind=f), intent(in) :: rmass(NBIN, NGROUP)
    real(kind=f), intent(in) :: numden(NZ, NELEM, NBIN)
    real(kind=f), intent(in) :: mmr_gas(NZ, NGAS)
    real(kind=f), intent(in) :: svpliq(NZ, NGAS)
    real(kind=f), intent(in) :: zsubsteps(NZ)
    real(kind=f), intent(in) :: rhompe(NZ, NBIN, NELEM)
    real(kind=f), intent(in) :: rnucpe(NZ, NBIN, NELEM)
    real(kind=f), intent(in) :: growpe(NZ, NBIN, NELEM)
    real(kind=f), intent(in) :: evappe(NZ, NBIN, NELEM)
    real(kind=f), intent(in) :: rnuclg(NZ, NBIN, NGROUP)
    real(kind=f), intent(in) :: growlg(NZ, NBIN, NGROUP)
    real(kind=f), intent(in) :: evaplg(NZ, NBIN, NGROUP)
    real(kind=f), intent(in) :: corefrac(NZ, NBIN, NGROUP)
    real(kind=f), intent(in) :: pflux(NZ+1, NBIN, NELEM)
    real(kind=f), intent(in) :: gflux(NZ+1, NGAS)

    integer      :: NZP1
    real(kind=f) :: time, lon_frac

    NZP1 = NZ + 1
    time = real(istep, kind=f) * dtime
    if (NLONGITUDE > 0) then
      lon_frac = current_step / real(NLONGITUDE, kind=f) * 360._f
    else
      lon_frac = 0._f
    end if

    open(unit=lun_bin, file=bin_filename, access='stream', form='unformatted', status='replace')

    write(lun_bin) NZ, NGROUP, NELEM, NBIN, NGAS, istep, IS_2D, NZP1
    write(lun_bin) time
    write(lun_bin) r
    write(lun_bin) rmass
    write(lun_bin) numden
    write(lun_bin) mmr_gas
    write(lun_bin) svpliq
    write(lun_bin) zsubsteps
    write(lun_bin) rhompe
    write(lun_bin) rnucpe
    write(lun_bin) growpe
    write(lun_bin) evappe
    write(lun_bin) rnuclg
    write(lun_bin) growlg
    write(lun_bin) evaplg
    write(lun_bin) corefrac
    write(lun_bin) current_distance, rotation_counter, current_step, lon_frac
    write(lun_bin) pflux
    write(lun_bin) gflux

    close(unit=lun_bin)

  end subroutine write_binary_output

end module output_binary_mod
