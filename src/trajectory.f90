!
! Copyright (C) 2001-2019 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!-----------------------------------------------------------------------
SUBROUTINE trajectoryXYZ
!-----------------------------------------------------------------------
  USE constants,      ONLY : bohr_radius_angs
  USE io_global,      ONLY : ionode
  USE io_files,       ONLY : prefix     
  USE ions_base,      ONLY : nat, ityp, tau, atm
  USE control_flags,  ONLY : istep
  USE ener,           ONLY : etot
  USE cell_base,      ONLY : alat
  USE dynamics_module,ONLY : vel
  integer :: ia
   
  if(ionode) then
     open(117,file="trajectory-"//trim(prefix)//".xyz",status="unknown",position='APPEND')
     write(117,'(I5)') nat
     write(117,'("# Step: ",I5,5x,"Total energy: ",F17.8,5x,"Ry")') istep-1, etot
       do ia = 1, nat
          write( 117, '(A3,3X,6F14.9)') atm(ityp(ia)),tau(:,ia)*alat*bohr_radius_angs, vel(:,ia)
       enddo
     close(117)
  endif

  return
!-----------------------------------------------------------------------
END SUBROUTINE trajectoryXYZ
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
SUBROUTINE save_rho(istep)
!-----------------------------------------------------------------------
  USE kinds,            ONLY : dp
  USE cell_base,        ONLY : at, bg, omega, alat, celldm, ibrav
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE run_info,         ONLY : title 
  USE fft_base,         ONLY : dfftp
  USE scatter_mod,      ONLY : gather_grid
  USE io_global,        ONLY : stdout, ionode
  USE io_files,         ONLY : prefix     
  USE scf,              ONLY : rho, vltot, v
  implicit none
  integer, intent(in) :: istep
  character(256) :: filename
  real(dp), allocatable :: raux(:)
#if defined(__MPI)
  real(dp), allocatable :: raux1(:)
#endif

  allocate(raux(dfftp%nnr))
#if defined(__MPI)
  allocate(raux1(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
#endif

  if (ionode) then
     write(filename,'(''rho-'',A,''-'',I9.9,''.xsf'')') trim(prefix), istep
     write(stdout,'(5X,''writing density to file: '',A)') trim(filename)
     open(unit=118,file=trim(filename),status='unknown')
     call xsf_struct (alat, at, nat, tau, atm, ityp, 118)
  endif

  raux(:) = rho%of_r(:,1)
#if defined(__MPI)
  call gather_grid (dfftp, raux, raux1)
#endif

  if (ionode) then
#if defined(__MPI)
    call xsf_fast_datagrid_3d(raux1, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                              dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, at, alat, 118)
#else
    call xsf_fast_datagrid_3d(raux, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                              dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, at, alat, 118)
#endif
  endif

  deallocate(raux)
#if defined(__MPI)
  deallocate(raux1)
#endif

!-----------------------------------------------------------------------
END SUBROUTINE save_rho
!-----------------------------------------------------------------------

