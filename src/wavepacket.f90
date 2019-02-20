!
! Copyright (C) 2001-2019 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE setup_wavepacket
  !-----------------------------------------------------------------------
  !
  ! ... Setup a wavepacket in real space. The wavepacket is moving
  ! ... in the -z direction.
  !
  USE kinds,            ONLY : dp
  USE mp,               ONLY : mp_sum
  USE mp_pools,         ONLY : inter_pool_comm
  USE mp_bands,         ONLY : me_bgrp, intra_bgrp_comm
  USE fft_base,         ONLY : dffts
  USE cell_base,        ONLY : alat, at, bg
  USE io_global,        ONLY : stdout
  USE io_files,         ONLY : nwordwfc, iunwfc
  USE wavefunctions,    ONLY : evc, psic
  USE wvfct,            ONLY : npw, wg, nbnd
  USE klist,            ONLY : ngk, igk_k
  USE fft_interfaces,   ONLY : invfft, fwfft
  USE buffers
  USE tddft_module
  IMPLICIT NONE
  real(dp), parameter :: bohrradius = 0.52917721d0
  real(dp), parameter :: rydberg = 13.605693d0
  real(dp) :: wp_k
  real(dp) :: inv_nr1s, inv_nr2s, inv_nr3s
  integer :: ia, i, j, k, j0, k0, idx, ir, ir_end, ipol, ik
  real(dp) :: wp, phase, r(3), norm
  complex(dp), external :: zdotc

  ! convert to atomic units
  wp_pos = wp_pos / bohrradius
  wp_d = wp_d / bohrradius
  wp_ekin = wp_ekin/rydberg
  wp_k = sqrt(wp_ekin)

  ! put WP into bands (TODO: spin-polarized case) check nbnd is enough
  ik = 1
  nbnd_occ(ik) = nbnd_occ(ik) + 1
  wp_ibnd = nbnd_occ(ik)

  ! force fixed occupations
  tfixed_occ = .true.
  wg(wp_ibnd,ik) = 1.d0 ! non-polarized
  allocate(f_inp(nbnd,ik))
  f_inp = wg
  
  ! print info
  write(stdout,*)
  write(stdout,'(5X,''===== Electron wavepacket ======'')')
  write(stdout,'(5X,''initial position (au) :'',3(F12.6,2X))') wp_pos
  write(stdout,'(5X,''spread           (au) :'',3(F12.6,2X))') wp_d
  write(stdout,'(5X,''kinetic energy (Ry)   :'',F12.6)') wp_ekin
  write(stdout,'(5X,''momentum//z (au^-1)   :'',F12.6)') -wp_k
  write(stdout,'(5X,''wavepacket band       :'',I5)') wp_ibnd
  write(stdout,*)

  ! wavefunction (smooth) grid
  inv_nr1s = 1.d0 / real(dffts%nr1,dp)
  inv_nr2s = 1.d0 / real(dffts%nr2,dp)
  inv_nr3s = 1.d0 / real(dffts%nr3,dp)

  ! loop in the charge array
  psic = (0.d0,0.d0)
  j0 = dffts%my_i0r2p; k0 = dffts%my_i0r3p
  do ir = 1, dffts%nr1x*dffts%my_nr2p*dffts%my_nr3p
     ! ... three dimensional indexes
     idx = ir - 1
     k   = idx / (dffts%nr1x*dffts%my_nr2p)
     idx = idx - (dffts%nr1x*dffts%my_nr2p)*k
     k   = k + k0
     j   = idx / dffts%nr1x
     idx = idx - dffts%nr1x * j
     j   = j + j0
     i   = idx

     do ipol = 1, 3
       r(ipol) = real(i,dp)*inv_nr1s*at(ipol,1) + &
                 real(j,dp)*inv_nr2s*at(ipol,2) + &
                 real(k,dp)*inv_nr3s*at(ipol,3)
     enddo
     r = r * alat

     ! gaussian envelope
     wp = 1.d0
     wp = wp*exp(-0.5d0 * (r(1)-wp_pos(1))**2 / wp_d(1)**2)
     wp = wp*exp(-0.5d0 * (r(2)-wp_pos(2))**2 / wp_d(2)**2)
     wp = wp*exp(-0.5d0 * (r(3)-wp_pos(3))**2 / wp_d(3)**2)

     ! wp momentum
     phase = -wp_k * r(3)
     psic(ir) = wp * cmplx(cos(phase),sin(phase),dp)
  enddo

  ! transform to real space
  call fwfft ('Wave', psic, dffts)

  ! load evc from disk
  call get_buffer (evc, nwordwfc, iunwfc, ik)
  npw = ngk(ik)
  evc(1:npw,wp_ibnd) = psic(dffts%nl(igk_k(1:npw,ik)))
  
  ! normalize
  norm = dble( zdotc(npw, evc(1,wp_ibnd), 1, evc(1,wp_ibnd),1) )
  call mp_sum(norm, intra_bgrp_comm)
  evc(1:npw,wp_ibnd) = evc(1:npw,wp_ibnd) / sqrt(real(norm))

  ! add wavepacket to evc
  call save_buffer (evc, nwordwfc, iunwfc, ik)

END SUBROUTINE setup_wavepacket

