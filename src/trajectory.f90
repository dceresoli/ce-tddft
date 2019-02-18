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





