module output_ascii_mod
  use carma_precision_mod
  implicit none

contains

  subroutine write_ascii_output(lun, lunp, lunf, lunfp, lunrates, lunratesp, &
      NZ, NELEM, NBIN, NGAS, NGROUP, istep, IS_2D, NLONGITUDE, dtime, &
      current_distance, rotation_counter, current_step, &
      numden, mmr_gas, svpliq, zsubsteps, wtmol_gas, wtmol_air, &
      rhompe, rnucpe, growpe, evappe, rnuclg, growlg, evaplg, corefrac)

    integer,      intent(in) :: lun, lunp, lunf, lunfp, lunrates, lunratesp
    integer,      intent(in) :: NZ, NELEM, NBIN, NGAS, NGROUP
    integer,      intent(in) :: istep, IS_2D, NLONGITUDE
    real(kind=f), intent(in) :: dtime
    real(kind=f), intent(in) :: current_distance, rotation_counter, current_step
    real(kind=f), intent(in) :: numden(NZ, NELEM, NBIN)
    real(kind=f), intent(in) :: mmr_gas(NZ, NGAS)
    real(kind=f), intent(in) :: svpliq(NZ, NGAS)
    real(kind=f), intent(in) :: zsubsteps(NZ)
    real(kind=f), intent(in) :: wtmol_gas(NGAS)
    real(kind=f), intent(in) :: wtmol_air(NZ)
    real(kind=f), intent(in) :: rhompe(NZ, NBIN, NELEM)
    real(kind=f), intent(in) :: rnucpe(NZ, NBIN, NELEM)
    real(kind=f), intent(in) :: growpe(NZ, NBIN, NELEM)
    real(kind=f), intent(in) :: evappe(NZ, NBIN, NELEM)
    real(kind=f), intent(in) :: rnuclg(NZ, NBIN, NGROUP)
    real(kind=f), intent(in) :: growlg(NZ, NBIN, NGROUP)
    real(kind=f), intent(in) :: evaplg(NZ, NBIN, NGROUP)
    real(kind=f), intent(in) :: corefrac(NZ, NBIN, NGROUP)

    integer :: i, j, ielem, igas, igroup

    if (IS_2D .eq. 1) then
      write(lun,     '(5e25.5)') (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
      write(lunp,    '(5e25.5)') (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
      write(lunf,    '(5e25.5)') (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
      write(lunfp,   '(5e25.5)') (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
      write(lunrates,'(5e25.5)') (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
      write(lunratesp,'(5e25.5)') (istep)*dtime, current_distance, rotation_counter, current_step, current_step/NLONGITUDE * 360
    else
      write(lun,      '(f25.5)') (istep)*dtime
      write(lunp,     '(f25.5)') (istep)*dtime
      write(lunf,     '(f25.5)') (istep)*dtime
      write(lunfp,    '(f25.5)') (istep)*dtime
      write(lunrates, '(f25.5)') (istep)*dtime
      write(lunratesp,'(f25.5)') (istep)*dtime
    end if

    do j = 1, NBIN
      do i = 1, NZ

        write(lun,  '(i3,i4)', advance="no") j, i
        write(lunp, '(i3,i4)', advance="no") j, i

        do ielem = 1, NELEM
          write(lun,  '(e11.3)', advance="no") real(numden(i, ielem, j))
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

        write(lun,  '(f8.0)') zsubsteps(i)
        write(lunp, '(f8.0)') zsubsteps(i)

        do ielem = 1, NELEM
          write(lunrates, '(3i4,7e13.3e3)') j, i, ielem, &
            rhompe(i, j, ielem), rnucpe(i, j, ielem), &
            growpe(i, j, ielem), evappe(i, j, ielem)
        end do

        do igroup = 1, NGROUP
          write(lunrates, '(3i4,7e13.3e3)') j, i, igroup, &
            rnuclg(i, j, igroup), growlg(i, j, igroup), &
            evaplg(i, j, igroup), corefrac(i, j, igroup)
        end do

      end do
    end do

  end subroutine write_ascii_output

end module output_ascii_mod
