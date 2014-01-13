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
  implicit none
  integer, intent(in) :: istep
  real(dp) :: charge, ehart, etxc, vtxc, eth, etotefield

  call start_clock('updateH')
  
  ! calculate total charge density
  call deallocate_bec_type(becp)
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
    
  call stop_clock('updateH')
    
END SUBROUTINE update_hamiltonian

