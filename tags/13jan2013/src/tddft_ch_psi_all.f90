!
! Copyright (C) 2001-2010 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine tddft_ch_psi_all (n, h, ah, ee, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( S + ee * H), where, ee = i * dt/2
  ! to a vector h. The result is given in Ah.
  !
  USE kinds,        ONLY : dp
  USE wvfct,        ONLY : npwx, nbnd
  USE uspp,         ONLY : vkb
  USE becmod,       ONLY : becp, calbec
  USE mp_global,    ONLY : intra_pool_comm
  USE mp,           ONLY : mp_sum
  USE tddft_module, ONLY : nbnd_occ, alpha_pv

  implicit none

  integer :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point

  complex(DP) :: ee
  ! input: the eigenvalue

  complex(DP) :: h (npwx, m), ah (npwx, m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !
  integer :: ibnd, ig

  complex(dp), allocatable :: hpsi (:,:)
  
  call start_clock ('ch_psi')
  allocate (hpsi( npwx , m))
  hpsi (:,:) = cmplx(0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian H and ultrasoft projection operator S with the h vector
  !
  call h_psi (npwx, n, m, h, hpsi)
  call s_psi (npwx, n, m, h, ah)
  !
  call start_clock ('last')
  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = ah (ig, ibnd) +  ee * hpsi (ig, ibnd)
     enddo
  enddo

  deallocate (hpsi)
  call stop_clock ('last')
  call stop_clock ('ch_psi')
  return

end subroutine tddft_ch_psi_all
