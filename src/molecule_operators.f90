!
! Copyright (C) 2001-2010 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!==============================================================================
! Setup the position operator in real space. The origin is set to center
! of ionic charge. (r is in units of alat)
!==============================================================================   
SUBROUTINE setup_position_operator
  USE kinds,        ONLY : dp
  USE mp_global,    ONLY : me_pool, intra_pool_comm
  USE mp,           ONLY : mp_sum
  USE fft_base,     ONLY : dfftp, dffts
  USE ions_base,    ONLY : nat, tau, ityp, zv
  USE cell_base,    ONLY : at, bg, alat
  USE tddft_module, ONLY : r_pos, r_pos_s
  IMPLICIT NONE

  real(dp) :: zvtot, x0(3), r(3)
  real(dp) :: inv_nr1, inv_nr2, inv_nr3
  real(dp) :: inv_nr1s, inv_nr2s, inv_nr3s
  integer :: ia, i, j, k, index, index0, ir, ipol

  ! calculate the center of charge
  zvtot = 0.d0
  x0 = 0.d0
  do ia = 1, nat
     zvtot = zvtot + zv(ityp(ia))
     x0(:) = x0(:) + tau(:,ia)*zv(ityp(ia))
  enddo
  x0 = x0 / zvtot

  ! density (hard) grid
  inv_nr1 = 1.d0 / real(dfftp%nr1,dp)
  inv_nr2 = 1.d0 / real(dfftp%nr2,dp)
  inv_nr3 = 1.d0 / real(dfftp%nr3,dp)

  index0 = 0
#ifdef __PARA
  do i = 1, me_pool
    index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  enddo
#endif

  ! loop over real space grid
  do ir = 1, dfftp%nnr
    index = index0 + ir - 1
    k     = index / (dfftp%nr1x*dfftp%nr2x)
    index = index - (dfftp%nr1x*dfftp%nr2x)*k
    j     = index / dfftp%nr1x
    index = index - dfftp%nr1x*j
    i     = index

    do ipol = 1, 3
      r(ipol) = real(i,dp)*inv_nr1*at(ipol,1) + &
                real(j,dp)*inv_nr2*at(ipol,2) + &
                real(k,dp)*inv_nr3*at(ipol,3)
    enddo

    ! minimum image convenction
    r = r - x0
    call cryst_to_cart( 1, r, bg, -1 )
    r = r - anint(r)
    call cryst_to_cart( 1, r, at, 1 )
    
    r_pos(1:3,ir) = r(1:3)
  enddo


  ! wavefunction (smooth) grid
  inv_nr1s = 1.d0 / real(dffts%nr1,dp)
  inv_nr2s = 1.d0 / real(dffts%nr2,dp)
  inv_nr3s = 1.d0 / real(dffts%nr3,dp)

  index0 = 0
#ifdef __PARA
  do i = 1, me_pool
    index0 = index0 + dffts%nr1x * dffts%nr2x * dffts%npp(i)
  enddo
#endif

  ! loop over real space grid
  do ir = 1, dffts%nnr
    index = index0 + ir - 1
    k     = index / (dffts%nr1x*dffts%nr2x)
    index = index - (dffts%nr1x*dffts%nr2x)*k
    j     = index / dffts%nr1x
    index = index - dffts%nr1x*j
    i     = index

    do ipol = 1, 3
      r(ipol) = real(i,dp)*inv_nr1s*at(ipol,1) + &
                real(j,dp)*inv_nr2s*at(ipol,2) + &
                real(k,dp)*inv_nr3s*at(ipol,3)
    enddo

    ! minimum image convenction
    r = r - x0
    call cryst_to_cart( 1, r, bg, -1 )
    r = r - anint(r)
    call cryst_to_cart( 1, r, at, 1 )
    
    r_pos_s(1:3,ir) = r(1:3)
  enddo


END SUBROUTINE setup_position_operator


!==============================================================================
! Compute electron dipole moment using total charge density (Xiaofeng Qian)
!==============================================================================   
SUBROUTINE compute_electron_dipole(charge, dip)
  USE kinds,        ONLY : dp
  USE mp_global,    ONLY : me_pool, intra_pool_comm
  USE mp,           ONLY : mp_sum
  USE fft_base,     ONLY : dfftp
  USE cell_base,    ONLY : omega, alat
  USE scf,          ONLY : rho
  USE lsda_mod,     ONLY : nspin
  USE tddft_module, ONLY : r_pos
  IMPLICIT NONE
  real(dp), intent(out) :: charge(nspin), dip(3,nspin)
  integer :: ispin, ipol, nrp

  call start_clock('dipole')
  do ispin = 1, nspin
    charge(ispin) = sum(rho%of_r(:,ispin))
    do ipol = 1, 3
      dip(ipol,ispin) = sum(r_pos(ipol,:)*rho%of_r(:,ispin))
    enddo
  enddo
 
  nrp = dfftp%nr1 * dfftp%nr2 * dfftp%nr3 
  charge = charge * omega / real(nrp, dp)  
  dip = dip * omega / real(nrp, dp) * alat

#ifdef __PARA
  call mp_sum(charge, intra_pool_comm)
  call mp_sum(dip, intra_pool_comm)
#endif
  call stop_clock('dipole')
    
END SUBROUTINE compute_electron_dipole


!==============================================================================
! Compute electron quadrupoledipole moment using total charge density
!==============================================================================   
SUBROUTINE compute_electron_quadrupole(quad)
  USE kinds,        ONLY : dp
  USE mp_global,    ONLY : me_pool, intra_pool_comm
  USE mp,           ONLY : mp_sum
  USE fft_base,     ONLY : dfftp
  USE cell_base,    ONLY : omega, alat
  USE scf,          ONLY : rho
  USE lsda_mod,     ONLY : nspin
  USE tddft_module, ONLY : r_pos
  IMPLICIT NONE
  real(dp), intent(out) :: quad(3,3,nspin)
  integer :: ispin, ipol, jpol, nrp

  call start_clock('quadrupole')
  do ispin = 1, nspin
    do ipol = 1, 3
      do jpol = 1, 3
        quad(ipol,jpol,ispin) = sum(r_pos(ipol,:)*r_pos(jpol,:)*rho%of_r(:,ispin))
      enddo
    enddo
  enddo
    
  nrp = dfftp%nr1 * dfftp%nr2 * dfftp%nr3 
  quad = quad * omega / real(nrp, dp) * (alat*alat)

#ifdef __PARA
  call mp_sum(quad, intra_pool_comm)
#endif
  call stop_clock('quadrupole')
    
END SUBROUTINE compute_electron_quadrupole


