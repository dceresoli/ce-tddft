!
! Copyright (C) 2001-2019 Quantum-ESPRESSO group
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
  USE constants,                   ONLY : rytoev
  USE io_global,                   ONLY : stdout, ionode
  USE io_files,                    ONLY : nwordwfc, iunwfc
  USE ions_base,                   ONLY : nat, if_pos
  USE cell_base,                   ONLY : at, tpiba, alat
  USE wavefunctions,               ONLY : evc
  USE klist,                       ONLY : nks, nkstot, wk, xk, nelec, ngk, igk_k
  USE wvfct,                       ONLY : nbnd, npwx, wg, g2kin, current_k
  USE lsda_mod,                    ONLY : current_spin, lsda, isk, nspin
  USE becmod,                      ONLY : becp, calbec, allocate_bec_type, &
                                          is_allocated_bec_type, deallocate_bec_type
  USE mp_pools,                    ONLY : inter_pool_comm
  USE mp,                          ONLY : mp_sum, mp_barrier
  USE gvect,                       ONLY : g
  USE fft_base,                    ONLY : dfftp, dffts
  USE buffers,                     ONLY : get_buffer, save_buffer
  USE fixed_occ,                   ONLY : tfixed_occ 
  USE uspp,                        ONLY : nkb, vkb
  USE uspp_init,                   ONLY : init_us_2
  USE ener,                        ONLY : ef
  USE dynamics_module,             ONLY : vel, verlet, allocate_dyn_vars, deallocate_dyn_vars
  USE pwcom
  USE tddft_module
  IMPLICIT NONE

  !-- tddft variables ----------------------------------------------------
  complex(dp), allocatable :: tddft_psi(:,:,:), b(:,:)
  complex(dp), allocatable :: tddft_hpsi(:,:), tddft_spsi(:,:)
  complex(dp), allocatable :: tddft_Ppsi(:,:)        ! PAW correction to forces (Ehrenfest)
  real(dp), allocatable :: charge(:), dipole(:,:), quadrupole(:,:,:)
  complex(dp), allocatable :: circular(:,:), circular_local(:)

  integer :: istep, lter, flag_global
  integer :: ik, is, ibnd
  complex(dp) :: ee                     ! i*dt/2
  real(dp) :: anorm, wclock
  integer, external :: find_free_unit
  real(dp), external :: get_clock
  external tddft_ch_psi_all

  ! TODO: restart

  ! allocate memory
  call allocate_optical()

  ee = i_complex * dt / 2.d0  ! i*dt/2: do not change
  
  evc = cmplx(0.d0,0.d0)
  call tddft_cgsolver_initialize(npwx, nbnd_occ_max)
  if (iverbosity > 0) then
    write(stdout,'(5X,''Done with tddft_cgsolver_initialize'')')
    flush(stdout)
  endif
 
  ! print the legend
  if (ionode) call print_legend
  
  ! check if we are restarting
  if (l_tddft_restart) then

     do ik = 1, nks
        current_k = ik
        current_spin = isk(ik)
        npw = ngk(ik)
 
        ! initialize at k-point k
        call g2_kin(ik)
        call init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
        
        ! read wfcs from file and compute becp
        evc = (0.d0, 0.d0)
        call get_buffer (evc, nwordwfc, iunevcn, ik)
     end do
     
     call update_hamiltonian(-1)
 
     if (iverbosity > 0) write(stdout,'(5X,''Done with restart'')')
  endif


  if (ehrenfest) then
     call allocate_dyn_vars()
     vel(:,:) = 0.d0
     allocate(if_pos(3,nat)) ! Ehrenfest work around
     if_pos(:,:) = 1
  endif

  if (isave_rho /= 0) call save_rho(0)

  ! enter the main TDDFT loop
  wclock = get_clock('TDDFT')
  do istep = 1, nstep
     
    ! calculate dipole moment along x, y, and z direction
    call molecule_compute_dipole( charge, dipole )
    !call molecule_compute_quadrupole( quadrupole )

    ! loop over k-points     
    do ik = 1, nks
      current_k = ik
      current_spin = isk(ik)
      npw = ngk(ik)
 
      ! initialize at k-point k
      call g2_kin(ik)
      call init_us_2(npw, igk_k(1,ik), xk(1,ik), vkb)
      
      ! read wfcs from file and compute becp
      evc = (0.d0, 0.d0)
      if (istep == 1) then
        call get_buffer (evc, nwordwfc, iunwfc, ik)
      else
        call get_buffer (evc, nwordwfc, iunevcn, ik)
      endif
      if (.not. is_allocated_bec_type(becp)) call allocate_bec_type(nkb, nbnd, becp)
      call calbec( npw, vkb, evc, becp )
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
      call save_buffer (tddft_psi(:,:,1:2), nwordtdwfc, iuntdwfc, ik)
        
    enddo ! ik

#ifdef __MPI
    ! reduce over k-points
    if (l_circular_dichroism) call mp_sum(circular, inter_pool_comm)
    call mp_sum(charge, inter_pool_comm)
    call mp_sum(dipole, inter_pool_comm)
#endif

    ! update the hamiltonian (recompute charge and potential)
    call update_hamiltonian(istep)

    ! print observables
    if (ionode) then
      do is = 1, nspin
        write(stdout,'(''ENERGY '',2X,I6,5F16.8)') istep, etot, eband + deband, ehart, etxc+etxcc, ewld
        if (degauss > 0.d0) write(stdout,'(''EFERMI '',F16.8)') ef*rytoev
        write(stdout,'(''CHARGE '',I1,1X,I6,3E16.6)') is, istep, charge(is)
        write(stdout,'(''DIP    '',I1,1X,I6,3E16.6)') is, istep, dipole(:,is)
        if (iverbosity > 11) write(stdout,'(''CPUTIME'',F16.6)') get_clock('TDDFT') - wclock
        wclock = get_clock('TDDFT')
        !write(stdout,'(''QUAD   '',I1,1X,I6,9E18.9)') is, istep, quadrupole(:,:,is)
        !write(stdout,'(''ANG    '',I1,1X,I6,3E16.6)') is, istep, circular(:,is)
      enddo
    endif

    ! output fields
    if (isave_rho /= 0) then
        if (mod(istep,isave_rho) == 0) call save_rho(istep)
    endif

    ! Ehrenfest dynamics
    if (ehrenfest) then
       if (is_allocated_bec_type(becp)) call deallocate_bec_type(becp)
       call forces()
       call verlet()
       call trajectoryXYZ()
    endif
     
    flush(stdout)
     
  enddo      ! end of TDDFT loop
  write(stdout,*)

  ! finish  
  call tddft_cgsolver_finalize()
  call deallocate_optical()
  if (ehrenfest) call deallocate_dyn_vars() 
    
CONTAINS

  !====================================================================
  ! Print the legend key
  !====================================================================    
  SUBROUTINE print_legend
    write(stdout,'(5X,''Output quantities:'')')
    write(stdout,'(5X,''  ENERGY istep etot eband ehart etxc ewld'')')
    write(stdout,'(5X,''  CHARGE spin  istep  charge'')')
    write(stdout,'(5X,''  DIP    spin  istep  dipole(1:3)'')')
    !write(stdout,'(5X,''  QUAD   spin  istep  quadrupole(1:3,1:3)'')')
    !write(stdout,'(5X,''  ANG    spin  istep  Re[L(1:3)]  Im[L(1:3)]'')')
    write(stdout,*)
    flush(stdout)
  END SUBROUTINE print_legend

  
  !====================================================================
  ! Initialize and allocate memory
  !====================================================================    
  SUBROUTINE allocate_optical()
    USE becmod, ONLY : becp, allocate_bec_type, is_allocated_bec_type
    IMPLICIT NONE
    integer :: ik
    
    nbnd_occ_max = 0
    do ik = 1, nks
      if (nbnd_occ(ik) > nbnd_occ_max) nbnd_occ_max = nbnd_occ(ik)
    enddo

    if (.not. is_allocated_bec_type(becp)) call allocate_bec_type(nkb, nbnd, becp)
   
    allocate (tddft_psi (npwx,nbnd,2))
    allocate (tddft_hpsi(npwx,nbnd_occ_max))
    allocate (tddft_spsi(npwx,nbnd_occ_max))
    allocate (b(npwx,nbnd_occ_max))
    tddft_psi = (0.d0,0.d0)
    tddft_hpsi = (0.d0,0.d0)
    tddft_spsi = (0.d0,0.d0)
    b = (0.d0,0.d0)

    if (ehrenfest) allocate(tddft_Ppsi(npwx,nbnd))
 
    allocate (charge(nspin), dipole(3,nspin), quadrupole(3,3,nspin))
    allocate (circular(3,nspin), circular_local(3))
    charge = 0.d0
    dipole = 0.d0
    quadrupole = 0.d0
    circular = (0.d0, 0.d0)
    circular_local = (0.d0, 0.d0)

    allocate (r_pos(3,dfftp%nnr), r_pos_s(3,dffts%nnr))
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
    if (ehrenfest) deallocate(tddft_Ppsi)

  END SUBROUTINE deallocate_optical
   
   
  !====================================================================
  ! compute circular dichroism (EXPERIMENTAL, NORM-CONSERVING ONLY)
  !====================================================================      
  subroutine compute_circular_dichroism(cd)
    USE fft_base,               ONLY : dfftp
    USE fft_interfaces,         ONLY : invfft
    IMPLICIT NONE
    REAL(DP) :: xx(dfftp%nnr), yy(dfftp%nnr), zz(dfftp%nnr), gk
    INTEGER  :: ik, ibnd, i, j, k, j0, k0, idx, ir, ipol, ind, i_current_spin, ig
    complex(dp) :: p_psi(npwx), p_psi_r(dfftp%nnr, 3), cd(3, nspin), psic1(dfftp%nnr)
    
    xx(:) = 0.d0
    yy(:) = 0.d0
    zz(:) = 0.d0
    
    ! Loop in the charge array
    j0 = dfftp%my_i0r2p ; k0 = dfftp%my_i0r3p
    DO ir = 1, dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p
       !
       ! ... three dimensional indexes
       !
       idx = ir -1
       k   = idx / (dfftp%nr1x*dfftp%my_nr2p)
       idx = idx - (dfftp%nr1x*dfftp%my_nr2p)*k
       k   = k + k0
       j   = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x * j
       j   = j + j0
       i   = idx

             xx(ir) = &
                  dble( i-1 )/dble(dfftp%nr1) * at(1,1) * alat + &
                  dble( j-1 )/dble(dfftp%nr2) * at(1,2) * alat + &
                  dble( k-1 )/dble(dfftp%nr3) * at(1,3) * alat
             
             yy(ir) = &
                  dble( i-1 )/dble(dfftp%nr1) * at(2,1) * alat + &
                  dble( j-1 )/dble(dfftp%nr2) * at(2,2) * alat + &
                  dble( k-1 )/dble(dfftp%nr3) * at(2,3) * alat
             
             zz(ir) = &
                  dble( i-1 )/dble(dfftp%nr1) * at(3,1) * alat + &
                  dble( j-1 )/dble(dfftp%nr2) * at(3,2) * alat + &
                  dble( k-1 )/dble(dfftp%nr3) * at(3,3) * alat
             
    end do
    
    cd(:,:) = (0.d0, 0.d0)
    
    do ik = 1, nks
       
       if (nbnd_occ(ik) > 0) then
          
       i_current_spin = isk(ik)
       
       do ibnd = 1, nbnd_occ(ik)
          
          p_psi_r(:, :) = (0.d0, 0.d0)
          do ipol = 1, 3
             p_psi(:) = (0.d0, 0.d0)
             do ig = 1, npw
                gk = xk(ipol,ik) + g(ipol,igk_k(ig,ik))
                p_psi(ig) = gk * tpiba * tddft_psi(ig, ibnd, ik)
             end do
             psic1(:) = (0.d0, 0.d0)
             psic1(dffts%nl(igk_k(1:npw,ik))) = p_psi(:)
             call invfft('Wave', psic1, dfftp)
             p_psi_r(:,ipol) = psic1(:)
          end do
          
          ! transform wavefunction from reciprocal space into real space
          psic1(:) = (0.d0, 0.d0)
          psic1(dffts%nl(igk_k(1:npw,ik))) = tddft_psi(1:npw, ibnd, ik)
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
 

