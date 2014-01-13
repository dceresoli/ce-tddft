!
! Copyright (C) 2001-2014 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE apply_electric_field(tddft_psi)
  !-----------------------------------------------------------------------
  !
  ! ... Apply an electric field impuse at t = 0, a homogeneous phase-shift
  ! ... to each band
  !
  USE kinds,        ONLY : dp
  USE mp_global,    ONLY : me_pool, intra_pool_comm
  USE mp,           ONLY : mp_sum
  USE fft_base,     ONLY : dffts
  USE ions_base,    ONLY : nat, tau, ityp, zv
  USE cell_base,    ONLY : at, bg, alat
  USE gvecs,        ONLY : nls
  USE wvfct,        ONLY : current_k, igk, npw, npwx, nbnd
  USE io_files,     ONLY : nwordwfc, iunwfc
  USE buffers,      ONLY : save_buffer
  USE wavefunctions_module, ONLY : evc
  USE fft_base,       ONLY : dffts
  USE fft_interfaces, ONLY : invfft, fwfft
  USE tddft_module
  implicit none

  complex(dp), intent(out) :: tddft_psi(npwx,nbnd)
  integer :: ik, ibnd, ir
  complex(dp) :: psic(dffts%nnr)
  real(dp) :: phase 
  
  ik = current_k

  do ibnd = 1, nbnd_occ(ik)
    ! transform wavefunction from reciprocal space into real space
    psic = (0.d0, 0.d0)
    psic(nls(igk(1:npw))) = evc(1:npw, ibnd)
    call invfft ('Wave', psic, dffts)  

    do ir = 1, dffts%nnr
      phase = e_strength * r_pos_s(e_direction,ir) * alat  !! CAREFUL WITH ULTRASOFT
      psic(ir) = psic(ir) * cmplx(cos(phase),sin(phase),dp)
    enddo

    call fwfft ('Wave', psic, dffts)  
       
    evc(1:npw,ibnd) = psic(nls(igk(1:npw)))
  enddo
    
  call save_buffer (evc, nwordwfc, iunevcn, ik)
    
  tddft_psi(1:npwx,1:nbnd_occ(ik)) = evc(1:npwx,1:nbnd_occ(ik))
    
END SUBROUTINE apply_electric_field
  
