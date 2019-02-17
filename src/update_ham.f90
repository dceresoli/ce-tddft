!
! Copyright (C) 2001-2014 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE update_hamiltonian(istep)
  !-----------------------------------------------------------------------
  !
  ! ... Update the hamiltonian
  !
  USE kinds,         ONLY : dp
  USE ldaU,          ONLY : lda_plus_U
  USE scf,           ONLY : rho, rho_core, rhog_core, vltot, v, kedtau, vrs
  USE fft_base,      ONLY : dfftp
  USE gvecs,         ONLY : doublegrid
  USE io_global,     ONLY : stdout
  USE lsda_mod,      ONLY : nspin
  USE uspp,          ONLY : okvan, nkb
  USE dfunct,        ONLY : newd
  USE tddft_module,  ONLY : nupdate_Dnm, iverbosity
  USE becmod,        ONLY : becp, allocate_bec_type, deallocate_bec_type
  USE wvfct,         ONLY : nbnd
  USE scf,           ONLY : rho
  USE cell_base,     ONLY : tpiba2, alat, omega, at, bg
  USE ions_base,     ONLY : nsp, atm, zv, nat, tau, ityp
  USE gvect,         ONLY : ngm, gstart, g, gg, gcutm, eigts1, eigts2, eigts3
  USE control_flags, ONLY : gamma_only
  USE pwcom
  implicit none
  integer, intent(in) :: istep
  real(dp) :: charge, eth, etotefield
  real(dp), external :: ewald

  call start_clock('updateH')
  
  ! calculate total charge density
  call deallocate_bec_type(becp)
  rho%of_g(:,:) = (0.d0,0.d0)
  rho%of_r(:,:) = 0.d0
  call sum_band()
  call allocate_bec_type(nkb, nbnd, becp)

  if (lda_plus_U) then
    call new_ns
    if (iverbosity > 10) call write_ns()
  end if
    
  ! calculate HXC-potential
  call v_of_rho( rho, rho_core, rhog_core, ehart, etxc, vtxc, eth, etotefield, charge, v )
    
  ! calculate total local potential (external + scf)
  call set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)    
  
  ! calculate new D_nm matrix for ultrasoft pseudopotential
  if (okvan) then
    if (istep == -1 .or. ( (nupdate_Dnm /= 0 .and. mod(istep,nupdate_Dnm) == 0) ) ) then
      call newd()
      if (iverbosity > 10) write(stdout,'(5X,''call newd'')')
    endif
  endif

  ! calculate band energy and Ewald energy
  call calculate_eband( (istep == -1) )
  ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )

  ! calculate new energy
  etot = eband + ( etxc - etxcc ) + ewld + ehart
  call sum_energies
    
  call stop_clock('updateH')
    
END SUBROUTINE update_hamiltonian


!-----------------------------------------------------------------------
 SUBROUTINE calculate_eband(first)
  !-----------------------------------------------------------------------
  !
  ! ... Calculate band energy
  !
  USE constants,                   ONLY : rytoev
  USE kinds,                       ONLY : dp
  USE io_global,                   ONLY : stdout
  USE io_files,                    ONLY : nwordwfc, iunwfc
  USE klist,                       ONLY : nks, wk
  USE wvfct,                       ONLY : nbnd, npwx, npw, wg, current_k, et, g2kin
  USE lsda_mod,                    ONLY : current_spin, isk, nspin
  USE becmod,                      ONLY : becp, calbec, allocate_bec_type, &
                                          is_allocated_bec_type, deallocate_bec_type
  USE mp,                          ONLY : mp_sum, mp_barrier
  USE mp_pools,                    ONLY : intra_pool_comm, inter_pool_comm
  USE buffers,                     ONLY : get_buffer, save_buffer
  USE uspp,                        ONLY : nkb, vkb
  USE cell_base,                   ONLY : tpiba2, omega
  USE klist,                       ONLY : igk_k
  USE pwcom
  USE tddft_module
  implicit none
  logical, intent(in) :: first
  complex(dp), allocatable :: hpsi(:,:)
  complex(dp), allocatable :: evc_temp(:,:)
  integer :: ik, ibnd
  complex(dp), external :: zdotc
 
  if (.not. is_allocated_bec_type(becp)) call allocate_bec_type(nkb, nbnd, becp)
  allocate(hpsi(npwx,nbnd))
  allocate(evc_temp(npwx,nbnd))

  ! loop over k-points     
  eband = 0.d0

  do ik = 1, nks
    current_k = ik
    current_spin = isk(ik)
    npw = ngk(ik)
  
    call g2_kin(ik) 
    call init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)

    evc_temp = (0.d0, 0.d0)
    if (first) then 
        call get_buffer(evc_temp, nwordwfc, iunwfc, ik)
    else
        call get_buffer(evc_temp, nwordwfc, iunevcn, ik)
    endif

    hpsi = (0.d0,0.d0)
    call calbec(npw, vkb, evc_temp, becp)
    call h_psi(npwx, npw, nbnd, evc_temp, hpsi)
    
    do ibnd = 1, nbnd
        et(ibnd,ik) = real(zdotc(npw,hpsi(1,ibnd),1,evc_temp(1,ibnd),1),dp)
        eband = eband + wg(ibnd,ik)*et(ibnd,ik)
    enddo
  enddo

  call mp_sum(et, intra_pool_comm)
  call mp_sum(eband, intra_pool_comm)
  call mp_sum(eband, inter_pool_comm)
  
  eband = eband + delta_e()
   
  deallocate(hpsi)
  deallocate(evc_temp)
  call deallocate_bec_type(becp)
  
  CONTAINS
  
   !-----------------------------------------------------------------------
   FUNCTION delta_e()
     !-----------------------------------------------------------------------
     ! ... delta_e = - \int rho%of_r(r)  v%of_r(r)
     !               - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
     !               - \sum rho%ns       v%ns       [for LDA+U]
     !               - \sum becsum       D1_Hxc     [for PAW]
     USE scf,              ONLY : scf_type, rho, v
     USE funct,            ONLY : dft_is_meta
     USE fft_base,         ONLY : dfftp
     USE noncollin_module, ONLY : noncolin
     USE mp_bands,         ONLY : intra_bgrp_comm
     USE paw_variables,    ONLY : okpaw, ddd_paw
     USE ldaU,             ONLY : lda_plus_U
     IMPLICIT NONE
     REAL(DP) :: delta_e, delta_e_hub
     !
     delta_e = - SUM( rho%of_r(:,:)*v%of_r(:,:) )
     !
     IF ( dft_is_meta() ) &
        delta_e = delta_e - SUM( rho%kin_r(:,:)*v%kin_r(:,:) )
     !
     delta_e = omega * delta_e / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
     !
     CALL mp_sum( delta_e, intra_bgrp_comm )
     !
     if (lda_plus_u) then
       if (noncolin) then
         delta_e_hub = - real(SUM (rho%ns_nc(:,:,:,:)*v%ns_nc(:,:,:,:)),kind=dp)
         delta_e = delta_e + delta_e_hub
       else
         delta_e_hub = - SUM (rho%ns(:,:,:,:)*v%ns(:,:,:,:))
         if (nspin==1) delta_e_hub = 2.d0 * delta_e_hub
         delta_e = delta_e + delta_e_hub
       endif
     end if
     !
     IF (okpaw) delta_e = delta_e - SUM(ddd_paw(:,:,:)*rho%bec(:,:,:))
     !
     RETURN
     !
   END FUNCTION delta_e

END SUBROUTINE calculate_eband

