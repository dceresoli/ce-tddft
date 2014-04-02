!
! Copyright (C) 2001-2014 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine molecule_optical_absorption
  !----------------------------------------------------------------------
  !  ... Compute optical absorption spectrum by real-time TDDFT 
  !  ... References:
  !      (1) Phys. Rev. B 73, 035408 (2006)
  !      (2) http://www.netlib.org/linalg/html_templates/Templates.html
  !                                             Xiaofeng Qian, MIT (2008)
  !----------------------------------------------------------------------
  USE kinds,                       ONLY : dp
  USE io_global,                   ONLY : stdout, ionode
  USE io_files,                    ONLY : nwordwfc, iunwfc, iunigk
  USE ions_base,                   ONLY : nat, ntyp => nsp, ityp
  USE cell_base,                   ONLY : at, bg, omega, tpiba, tpiba2, alat
  USE wavefunctions_module,        ONLY : evc
  USE klist,                       ONLY : nks, nkstot, wk, xk, nelec, ngk
  USE wvfct,                       ONLY : nbnd, npwx, npw, igk, wg, g2kin, current_k, ecutwfc
  USE lsda_mod,                    ONLY : current_spin, lsda, isk, nspin
  USE becmod,                      ONLY : becp  
  USE mp_pools,                    ONLY : my_pool_id, inter_pool_comm, intra_pool_comm
  USE mp,                          ONLY : mp_sum, mp_barrier
  USE gvect,                       ONLY : ngm, g
  USE gvecs,                       ONLY : nls
  USE fft_base,                    ONLY : dfftp, dfftp
  USE buffers,                     ONLY : get_buffer, save_buffer
  USE fixed_occ,                   ONLY : tfixed_occ 
  USE uspp,                        ONLY : nkb, vkb, deeq
  USE ldaU,                        ONLY : lda_plus_U
  USE uspp_param,                  ONLY : nh
  USE scf,                         ONLY : rho, rho_core, rhog_core, vltot, v, vrs
  USE control_flags,               ONLY : tqr
  USE tddft_module

  IMPLICIT NONE

  !-- tddft variables ----------------------------------------------------
  complex(dp), allocatable :: tddft_psi(:,:,:), b(:,:)
  complex(dp), allocatable :: tddft_hpsi(:,:), tddft_spsi(:,:)
  real(dp), allocatable :: charge(:), dipole(:,:), quadrupole(:,:,:)
  complex(dp), allocatable :: circular(:,:), circular_local(:)

  integer :: istep, lter, flag_global
  integer :: ik, is, ibnd
  complex(dp) :: ee                     ! i*dt/2
  real(dp) :: anorm
  integer, external :: find_free_unit
  external tddft_ch_psi_all

  ! TODO: gk_sort
  ! TODO: restart

  ! allocate memory
  call allocate_optical()

  ee = i_complex * dt / 2.d0  ! i*dt/2: do not change
  
  evc = cmplx(0.d0,0.d0)
  call tddft_cgsolver_initialize(npwx, nbnd_occ_max)
  if (iverbosity > 0) then
    write(stdout,'(5X,''Done with tddft_cgsolver_initialize'')')
    call flush_unit(stdout)
  endif
 
  ! print the legend
  if (ionode) call print_legend
  
  ! check if we are restarting
  if (l_tddft_restart) then

     if (nks > 1) rewind (iunigk)
     do ik = 1, nks
        current_k = ik
        current_spin = isk(ik)
        
        ! initialize at k-point k 
        call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
        g2kin = g2kin * tpiba2
        call init_us_2(npw, igk, xk(1,ik), vkb)
        
        ! read wfcs from file and compute becp
        evc = (0.d0, 0.d0)
        call get_buffer (evc, nwordwfc, iunevcn, ik)
     end do
     
     call update_hamiltonian(-1)
  end if
 
  if (iverbosity > 0) then
    write(stdout,'(5X,''Done with restart'')')
    call flush_unit(stdout)
  endif


  ! enter the main TDDFT loop 
  do istep = 1, nstep
     
    ! calculate dipole moment along x, y, and z direction
    call molecule_compute_dipole( charge, dipole )
    !call molecule_compute_quadrupole( quadrupole )

    ! loop over k-points     
    if (nks > 1) rewind (iunigk)
    do ik = 1, nks
      current_k = ik
      current_spin = isk(ik)
        
      ! initialize at k-point k 
      call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
      g2kin = g2kin * tpiba2
      call init_us_2(npw, igk, xk(1,ik), vkb)
        
      ! read wfcs from file and compute becp
      evc = (0.d0, 0.d0)
      if (istep == 1) then
        call get_buffer (evc, nwordwfc, iunwfc, ik)
      else
        call get_buffer (evc, nwordwfc, iunevcn, ik)
      endif
      if ( (istep > 1) .or. (l_tddft_restart .and. (istep == 1)) ) then
        call get_buffer (tddft_psi, nwordtdwfc, iuntdwfc, ik)
      endif

      ! apply electric field
      if (istep == 1 .and. (.not. l_tddft_restart)) &
        call apply_electric_field(tddft_psi)
        
      ! calculate circular dichroism along x, y, and z direction
      if (l_circular_dichroism)  then
        circular_local = (0.d0, 0.d0)
        circular = (0.d0, 0.d0)
        call compute_circular_dichroism(circular_local)
        circular(1:3, current_spin) = circular_local(1:3)
      end if
        
      ! guess the wavefunction at the next timestep
      tddft_psi(:,:,1) = (0.d0,0.d0)
      do ibnd = 1, nbnd_occ(ik)
        tddft_psi(:,ibnd,1) = 2.d0*evc(:,ibnd) - tddft_psi(:,ibnd,2)
      enddo

      ! calculate H |psi_current>, S |psi_current>
      call h_psi(npwx, npw, nbnd_occ(ik), evc, tddft_hpsi)
      call s_psi(npwx, npw, nbnd_occ(ik), evc, tddft_spsi)
        
      ! calculate (S - H*dt*i/2) |\psi_current>
      b = (0.d0, 0.d0)
      b(1:npw, 1:nbnd_occ(ik)) = tddft_spsi(1:npw,1:nbnd_occ(ik)) - ee * tddft_hpsi(1:npw,1:nbnd_occ(ik))
        
      ! solve A * x = b
      call tddft_cgsolver(tddft_ch_psi_all, b, tddft_psi(:,:,1), npwx, npw, &
                          conv_threshold, ik, lter, flag_global, anorm, &
                          nbnd_occ(ik), ee)
        
      ! update the wavefunctions
      tddft_psi(:,:,2) = tddft_psi(:,:,1)
      evc(:,1:nbnd_occ(ik)) = tddft_psi(:,1:nbnd_occ(ik),1)

      ! save wavefunctions to disk
      call save_buffer (evc, nwordwfc, iunevcn, ik)
      call save_buffer (tddft_psi, nwordtdwfc, iuntdwfc, ik)
        
    enddo ! ik

#ifdef __PARA
    ! reduce over k-points
    if (l_circular_dichroism) call mp_sum(circular, inter_pool_comm)
    call mp_sum(charge, inter_pool_comm)
    call mp_sum(dipole, inter_pool_comm)
#endif

    ! print observables
    if (ionode) then
      do is = 1, nspin
        write(stdout,'(''CHARGE '',I1,1X,I6,3E16.6)') is, istep, charge(is)
        write(stdout,'(''DIP    '',I1,1X,I6,3E16.6)') is, istep, dipole(:,is)
        !write(stdout,'(''QUAD   '',I1,1X,I6,9E18.9)') is, istep, quadrupole(:,:,is)
        !write(stdout,'(''ANG    '',I1,1X,I6,3E16.6)') is, istep, circular(:,is)
      enddo
    endif
     
    ! update the hamiltonian (recompute charge and potential)
    call update_hamiltonian(istep)
     
    call flush_unit(stdout)
     
  enddo      ! end of TDDFT loop

  ! finish  
  call tddft_cgsolver_finalize()
  call deallocate_optical()
  
    
CONTAINS

  !====================================================================
  ! Print the legend key
  !====================================================================    
  SUBROUTINE print_legend
    write(stdout,'(5X,''Output quantities:'')')
    write(stdout,'(5X,''  CHARGE spin  istep  charge'')')
    write(stdout,'(5X,''  DIP    spin  istep  dipole(1:3)'')')
    write(stdout,'(5X,''  QUAD   spin  istep  quadrupole(1:3,1:3)'')')
    write(stdout,'(5X,''  ANG    spin  istep  Re[L(1:3)]  Im[L(1:3)]'')')
    write(stdout,*)
    call flush_unit(stdout)
  END SUBROUTINE print_legend

  
  !====================================================================
  ! Initialize and allocate memory
  !====================================================================    
  SUBROUTINE allocate_optical()
    USE becmod, ONLY : becp, allocate_bec_type
    IMPLICIT NONE
    integer :: ik
    
    nbnd_occ_max = 0
    do ik = 1, nks
      if (nbnd_occ(ik) > nbnd_occ_max) nbnd_occ_max = nbnd_occ(ik)
    enddo

    call allocate_bec_type(nkb, nbnd, becp)
   
    allocate (tddft_psi (npwx,nbnd,2))
    allocate (tddft_hpsi(npwx,nbnd_occ_max))
    allocate (tddft_spsi(npwx,nbnd_occ_max))
    allocate (b(npwx,nbnd_occ_max))
    tddft_psi = (0.d0,0.d0)
    tddft_hpsi = (0.d0,0.d0)
    tddft_spsi = (0.d0,0.d0)
    b = (0.d0,0.d0)

    allocate (charge(nspin), dipole(3,nspin), quadrupole(3,3,nspin))
    allocate (circular(3,nspin), circular_local(3))
    charge = 0.d0
    dipole = 0.d0
    quadrupole = 0.d0
    circular = (0.d0, 0.d0)
    circular_local = (0.d0, 0.d0)

    allocate (r_pos(3,dfftp%nnr), r_pos_s(3,dfftp%nnr))
    call molecule_setup_r
    
  END SUBROUTINE allocate_optical
  
  
  !====================================================================
  ! Deallocate memory
  !====================================================================    
  SUBROUTINE deallocate_optical()
    USE becmod, ONLY : becp, deallocate_bec_type
    IMPLICIT NONE

    call deallocate_bec_type(becp)
    deallocate (tddft_psi, tddft_hpsi, tddft_spsi, b)
    deallocate (charge, dipole, quadrupole, circular, circular_local)
    deallocate (r_pos, r_pos_s)

  END SUBROUTINE deallocate_optical
   
   
  !====================================================================
  ! compute circular dichroism (EXPERIMENTAL, NORM-CONSERVING ONLY)
  !====================================================================      
  subroutine compute_circular_dichroism(cd)
    USE fft_base,               ONLY : dfftp
    USE fft_interfaces,         ONLY : invfft
    USE mp_global,              ONLY : me_pool
    IMPLICIT NONE
    REAL(DP) :: xx(dfftp%nnr), yy(dfftp%nnr), zz(dfftp%nnr), gk
    INTEGER  :: ik, ibnd, i, ii, jj, kk, index0, index, ir, ipol, ind, i_current_spin, ig
    complex(dp) :: p_psi(npwx), p_psi_r(dfftp%nnr, 3), cd(3, nspin), psic1(dfftp%nnr)
    
    xx(:) = 0.d0
    yy(:) = 0.d0
    zz(:) = 0.d0
    
    index0 = 0

#ifdef __PARA
  do i = 1, me_pool
    index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  enddo
#endif

  ! loop over real space grid
  do ir = 1, dfftp%nnr
    index = index0 + ir - 1
    kk     = index / (dfftp%nr1x*dfftp%nr2x)
    index = index - (dfftp%nr1x*dfftp%nr2x)*kk
    jj     = index / dfftp%nr1x
    index = index - dfftp%nr1x*jj
    ii     = index

             xx(ir) = &
                  dble( ii-1 )/dble(dfftp%nr1) * at(1,1) * alat + &
                  dble( jj-1 )/dble(dfftp%nr2) * at(1,2) * alat + &
                  dble( kk-1 )/dble(dfftp%nr3) * at(1,3) * alat
             
             yy(ir) = &
                  dble( ii-1 )/dble(dfftp%nr1) * at(2,1) * alat + &
                  dble( jj-1 )/dble(dfftp%nr2) * at(2,2) * alat + &
                  dble( kk-1 )/dble(dfftp%nr3) * at(2,3) * alat
             
             zz(ir) = &
                  dble( ii-1 )/dble(dfftp%nr1) * at(3,1) * alat + &
                  dble( jj-1 )/dble(dfftp%nr2) * at(3,2) * alat + &
                  dble( kk-1 )/dble(dfftp%nr3) * at(3,3) * alat
             
    end do
    
    cd(:,:) = (0.d0, 0.d0)
    
    do ik = 1, nks
       
       if (nbnd_occ(ik) > 0) then
          
       i_current_spin = isk(ik)
       call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
       g2kin(:) = g2kin(:) * tpiba2
       
       do ibnd = 1, nbnd_occ(ik)
          
          p_psi_r(:, :) = (0.d0, 0.d0)
          do ipol = 1, 3
             p_psi(:) = (0.d0, 0.d0)
             do ig = 1, npw
                gk = xk(ipol,ik) + g(ipol,igk(ig))
                p_psi(ig) = gk * tpiba * tddft_psi(ig, ibnd, ik)
             end do
             psic1(:) = (0.d0, 0.d0)
             psic1(nls(igk(1:npw))) = p_psi(:)
             call invfft('Wave', psic1, dfftp)
             p_psi_r(:,ipol) = psic1(:)
          end do
          
          ! transform wavefunction from reciprocal space into real space
          psic1(:) = (0.d0, 0.d0)
          psic1(nls(igk(1:npw))) = tddft_psi(1:npw, ibnd, ik)
          call invfft('Wave', psic1, dfftp)
          
          do ind = 1, dfftp%nnr
             cd(1, i_current_spin) = cd(1, i_current_spin) + &
                  conjg(psic1(ind)) * ( yy(ind) * p_psi_r(ind,3) - zz(ind) * p_psi_r(ind,2) )
             cd(2, i_current_spin) = cd(2, i_current_spin) + &
                  conjg(psic1(ind)) * ( zz(ind) * p_psi_r(ind,1) - xx(ind) * p_psi_r(ind,3) )
             cd(3, i_current_spin) = cd(3, i_current_spin) + &
                  conjg(psic1(ind)) * ( xx(ind) * p_psi_r(ind,2) - yy(ind) * p_psi_r(ind,1) )
          end do
          
          
       end do

       end if
       
    end do
    cd = cd  / dble(dfftp%nnr)
    
    
    RETURN
  end subroutine compute_circular_dichroism

END SUBROUTINE molecule_optical_absorption
 

