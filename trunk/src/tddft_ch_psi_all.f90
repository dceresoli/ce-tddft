!
! Copyright (C) 2001-2015 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
#ifdef __BANDS
subroutine tddft_ch_psi_all (n, h, ah, ee, ik, m, ibnd_start, ibnd_end)
#else
subroutine tddft_ch_psi_all (n, h, ah, ee, ik, m)
#endif
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( S + ee * H), where, ee = i * dt/2
  ! to a vector h. The result is given in Ah.
  !
  USE kinds,        ONLY : dp
  USE wvfct,        ONLY : npwx, nbnd
  USE uspp,         ONLY : vkb, nkb
  USE becmod,       ONLY : becp, calbec
  USE tddft_module, ONLY : nbnd_occ, alpha_pv
  USE mp_pools,     ONLY : intra_pool_comm
  USE mp,           ONLY : mp_sum
#ifdef __BANDS
  USE mp_bands,     ONLY : intra_bgrp_comm
#endif
  implicit none

  integer :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point
#ifdef __BANDS
  integer, intent(in) :: ibnd_start, ibnd_end
#endif
  complex(DP) :: ee
  ! input: i*dt/2

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

  !   compute the product of the hamiltonian H and ultrasoft projection operator S with the h vector
  !   TODO: band version
  call h_psi (npwx, n, m, h, hpsi)
  call s_psi (npwx, n, m, h, ah)

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
