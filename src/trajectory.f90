!
! Copyright (C) 2001-2019 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!-----------------------------------------------------------------------
SUBROUTINE sum_energies
!-----------------------------------------------------------------------
  USE kinds,              ONLY : dp
  USE paw_variables,      ONLY : okpaw
  USE ldaU,               ONLY : lda_plus_u, eth
  USE control_flags,      ONLY : llondon, ldftd3, lxdm, ts_vdw
  USE xdm_module,         ONLY : energy_xdm
  USE extfield,           ONLY : tefield, etotefield
  USE tsvdw_module,       ONLY : EtsvdW
  USE plugin_variables,   ONLY : plugin_etot
  USE pwcom
  implicit none
  real(dp) :: eext = 0.d0

  if (okpaw) etot = etot + epaw
  if (lda_plus_u) etot = etot + eth

  if (llondon) then
     etot = etot + elondon
     hwf_energy = hwf_energy + elondon
  endif

  if (lxdm) then
     exdm = energy_xdm()
     etot = etot + exdm
     hwf_energy = hwf_energy + exdm
  endif

  if (ts_vdw) then
     ! factor 2 converts from Ha to Ry units
     etot = etot + 2.0d0*EtsvdW
     hwf_energy = hwf_energy + 2.0d0*EtsvdW
  endif

  if (tefield) then
     etot = etot + etotefield
     hwf_energy = hwf_energy + etotefield
  endif

  ! adds possible external contribution from plugins to the energy
  etot = etot + plugin_etot

  return

!-----------------------------------------------------------------------
END SUBROUTINE sum_energies
!-----------------------------------------------------------------------


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





