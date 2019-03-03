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
  USE uspp,          ONLY : okvan
  USE dfunct,        ONLY : newd
  USE tddft_module,  ONLY : nupdate_Dnm, iverbosity
  USE wvfct,         ONLY : nbnd
  USE scf,           ONLY : rho
  USE cell_base,     ONLY : alat, omega, at, bg
  USE ions_base,     ONLY : nsp, zv, nat, tau, ityp
  USE gvect,         ONLY : ngm, gstart, g, gg, gcutm
  USE control_flags, ONLY : gamma_only
  USE becmod,        ONLY : becp, calbec, allocate_bec_type, &
                            is_allocated_bec_type, deallocate_bec_type
  USE uspp,          ONLY : nkb
  USE pwcom
  implicit none
  integer, intent(in) :: istep
  real(dp) :: charge, eth, etotefield
  real(dp), external :: ewald, delta_eband

  call start_clock('updateH')
  
  ! calculate total charge density
  rho%of_g(:,:) = (0.d0,0.d0)
  rho%of_r(:,:) = 0.d0
  if (okvan .and. is_allocated_bec_type(becp)) call deallocate_bec_type(becp)
  call sum_band()

  if (lda_plus_U) then
    call new_ns
    if (iverbosity > 10) call write_ns()
  end if
    
  ! calculate HXC-potential
  call v_of_rho( rho, rho_core, rhog_core, ehart, etxc, vtxc, eth, etotefield, charge, v )
    
  ! calculate total local potential (external + scf)
  call setlocal
  call set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)    
  
  ! calculate new D_nm matrix for ultrasoft pseudopotential
  if (okvan) then
    if (istep == -1 .or. ( (nupdate_Dnm /= 0 .and. mod(istep,nupdate_Dnm) == 0) ) ) then
      call newd()
      if (iverbosity > 10) write(stdout,'(5X,''call newd'')')
    endif
  endif

  ! calculate band energy and Ewald energy
  deband = delta_eband()
  ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )

  ! calculate new energy
  etot = eband + deband + ( etxc - etxcc ) + ewld + ehart
  call sum_energies

  ! intialization step
  if (istep == -1) &
      write(stdout,'(''ENERGY '',2X,I6,5F16.8)') istep, etot, eband + deband, ehart, etxc+etxcc, ewld
    
  call stop_clock('updateH')

END SUBROUTINE update_hamiltonian


!-----------------------------------------------------------------------
FUNCTION delta_eband() RESULT(delta_e)
  !-----------------------------------------------------------------------
  ! ... delta_e = - \int rho%of_r(r)  v%of_r(r)
  !               - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
  !               - \sum rho%ns       v%ns       [for LDA+U]
  !               - \sum becsum       D1_Hxc     [for PAW]
  USE kinds,            ONLY : dp
  USE scf,              ONLY : scf_type, rho, v
  USE funct,            ONLY : dft_is_meta
  USE fft_base,         ONLY : dfftp
  USE noncollin_module, ONLY : noncolin
  USE mp,               ONLY : mp_sum, mp_barrier
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE paw_variables,    ONLY : okpaw, ddd_paw
  USE ldaU,             ONLY : lda_plus_U
  USE cell_base,        ONLY : omega
  USE lsda_mod,         ONLY : nspin
  IMPLICIT NONE
  REAL(DP) :: delta_e, delta_e_hub
  integer :: ir

  delta_e = 0.d0
  IF ( nspin==2 ) THEN
     !
     DO ir = 1,dfftp%nnr
       delta_e = delta_e - ( rho%of_r(ir,1) + rho%of_r(ir,2) ) * v%of_r(ir,1) &  ! up
                         - ( rho%of_r(ir,1) - rho%of_r(ir,2) ) * v%of_r(ir,2)    ! dw
     ENDDO
     delta_e = 0.5_dp*delta_e
     !
  ELSE
     delta_e = - SUM( rho%of_r(:,:)*v%of_r(:,:) )
  ENDIF
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
      delta_e_hub = - SUM (rho%ns_nc(:,:,:,:)*v%ns_nc(:,:,:,:))
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
END FUNCTION delta_eband


!-----------------------------------------------------------------------
SUBROUTINE sum_energies
!-----------------------------------------------------------------------
  USE kinds,              ONLY : dp
  USE paw_variables,      ONLY : okpaw
  USE ldaU,               ONLY : lda_plus_u, eth
  USE control_flags,      ONLY : llondon, lxdm, ts_vdw
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
  etot = etot + plugin_etot + eext

  return

!-----------------------------------------------------------------------
END SUBROUTINE sum_energies
!-----------------------------------------------------------------------
