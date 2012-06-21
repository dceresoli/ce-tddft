!

! Copyright (C) 2001-2010 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
! Read input file
!-----------------------------------------------------------------------
SUBROUTINE tddft_readin()
  USE io_files,      ONLY : nd_nmbr, prefix, tmp_dir  
  USE io_global,     ONLY : ionode, stdout
  USE constants,     ONLY : bohr_radius_angs, au_sec
  USE tddft_module
  IMPLICIT NONE
  integer :: ios
  namelist /inputtddft/ job, prefix, tmp_dir, conv_threshold, iverbosity, &
                        dt, e_strength, e_direction, nstep, nupdate_Dnm, &
                        l_circular_dichroism, l_tddft_restart

  ! set default values
  job          = ''
  prefix       = 'pwscf'
  tmp_dir      = './scratch/'    
  iverbosity   = 0
  dt           = 2.d0                      ! time step (default: 2 attosecond)
  e_strength   = 0.01d0                    ! impulse electric field strength (default: 0.01/Ang)
  e_direction  = 1                         ! impulse electric field direction: 1-x 2-y 3-z
  conv_threshold = 1.0d-12                 ! convergence threshold    
  nstep        = 1000                      ! total time steps
  nupdate_Dnm  = 1                         ! update USPP Dnm every step
  l_circular_dichroism = .false.
  l_tddft_restart      = .false.

  if (ionode) then
    call input_from_file()
    read(5, inputtddft, err = 200, iostat = ios)
200 call errore('tddft_readin', 'reading inputtddft namelist', abs(ios))

    ! convert to atomic units
    e_strength = e_strength * bohr_radius_angs ! change from 1/Ang to 1/Bohr
    dt = dt * 1.d-18 / (2.d0*au_sec)           ! change from femto-second to a.u. (Rydberg unit)
  endif

  
  call tddft_bcast_input

END SUBROUTINE tddft_readin


  
!-----------------------------------------------------------------------
! Broadcast input data to all processors 
!-----------------------------------------------------------------------
SUBROUTINE tddft_bcast_input
  USE mp,            ONLY : mp_bcast
  USE io_files,      ONLY : prefix, tmp_dir
  USE tddft_module
  IMPLICIT NONE
  integer, parameter :: root = 0    

#ifdef __PARA
  call mp_bcast(job, root)
  call mp_bcast(prefix, root)
  call mp_bcast(tmp_dir, root)
  call mp_bcast(dt, root)
  call mp_bcast(e_strength, root)
  call mp_bcast(e_direction, root)
  call mp_bcast(conv_threshold, root)
  call mp_bcast(nstep, root)
  call mp_bcast(nupdate_Dnm , root)
  call mp_bcast(l_circular_dichroism, root)
  call mp_bcast(l_tddft_restart, root)
  call mp_bcast(iverbosity, root)
#endif

END SUBROUTINE tddft_bcast_input
  
  

!-----------------------------------------------------------------------
! Allocate memory for TDDFT
!-----------------------------------------------------------------------
SUBROUTINE tddft_allocate
  USE lsda_mod,      ONLY : nspin, lsda
  USE ions_base,     ONLY : ntyp => nsp
  USE klist,         ONLY : nkstot
  USE wvfct,         ONLY : btype, nbndx
  USE tddft_module
  IMPLICIT NONE

  ! needed by sum_band
  allocate(btype(nbndx,nkstot))
  btype = 1
    
END SUBROUTINE tddft_allocate
  
  

!-----------------------------------------------------------------------
! Print a short summary of the calculation
!-----------------------------------------------------------------------
SUBROUTINE tddft_summary
  USE io_global,     ONLY : stdout
  USE tddft_module
  IMPLICIT NONE
  
  write(stdout,*)
  write(stdout,'(5X,''Calculation type      : '',A12)') job
  write(stdout,'(5X,''Number or steps       : '',I12)') nstep
  write(stdout,'(5X,''Time step             : '',F12.4,'' rydberg_atomic_time'')') dt
  write(stdout,'(5X,''Electric field dir.   : '',I12,'' (1=x,2=y,3=z)'')') e_direction
  write(stdout,'(5X,''Electric field impulse: '',F12.4,'' bohrradius^-1'')') e_strength
  write(stdout,*)

  call flush_unit( stdout )

END SUBROUTINE tddft_summary
  
  

!-----------------------------------------------------------------------
! Open files needed for TDDFT
!-----------------------------------------------------------------------
SUBROUTINE tddft_openfil
  USE io_global,        ONLY : stdout
  USE io_files,         ONLY : diropn, seqopn
  USE basis,            ONLY : natomwfc, starting_wfc
  USE wvfct,            ONLY : nbnd, npwx
  USE ldaU,             ONLY : lda_plus_u
  USE klist,            ONLY : nks
  USE io_files,         ONLY : prefix, iunat, iunsat, iunwfc, iunigk, &
                               nwordwfc, nwordatwfc, tmp_dir, wfc_dir
  USE noncollin_module, ONLY : npol
  USE mp_global,        ONLY : kunit
  USE buffers,          ONLY : open_buffer
  USE control_flags,    ONLY : io_level 
  USE tddft_module   
  IMPLICIT NONE  
  LOGICAL            :: exst

  ! ... iunwfc=10: read/write wfc from/to file
  ! ... iunwfc=-1: copy wfc to/from RAM 
  if ( io_level > 0 ) then
    iunwfc = 10
  else
    iunwfc = -1
  end if
  call open_buffer( iunwfc, 'wfc', nwordwfc, nks, exst )
  call open_buffer( iunevcn, 'evcn', nwordwfc, nks, exst )
    
  nwordtdwfc = nbnd*npwx*npol
  call open_buffer( iuntdwfc, 'tdwfc', nwordtdwfc, nks, exst )
    
  ! ... Needed for LDA+U
  ! ... iunat  contains the (orthogonalized) atomic wfcs 
  ! ... iunsat contains the (orthogonalized) atomic wfcs * S
  ! ... iunocc contains the atomic occupations computed in new_ns
  ! ... it is opened and closed for each reading-writing operation  
  nwordatwfc = 2*npwx*natomwfc*npol
  if ( lda_plus_u ) then
     call diropn( iunat,  'atwfc',  nwordatwfc, exst )
     call diropn( iunsat, 'satwfc', nwordatwfc, exst )
  end if

  ! ... iunigk contains the number of PW and the indices igk
  call seqopn( iunigk, 'igk', 'UNFORMATTED', exst )

END SUBROUTINE tddft_openfil



!-----------------------------------------------------------------------
! Close files created by TDDFT
!-----------------------------------------------------------------------
SUBROUTINE tddft_closefil
  USE buffers,          ONLY : close_buffer
  USE tddft_module
  IMPLICIT NONE
    
  call close_buffer(iunevcn, 'KEEP')
  call close_buffer(iuntdwfc, 'KEEP')
 
END SUBROUTINE tddft_closefil
  


!-----------------------------------------------------------------------
! Print timings
!-----------------------------------------------------------------------
SUBROUTINE print_clock_tddft
  USE io_global,  ONLY : stdout
  IMPLICIT NONE

  write(stdout,*)
  call print_clock ('TDDFT') 
  write(stdout,*) '    INITIALIZATION: '
  call print_clock ('tddft_setup')
  write(stdout,*)
  call print_clock ('updateH')
  call print_clock ('greenf')
  call print_clock ('cgsolve')
  call print_clock ('ch_psi')
  call print_clock ('h_psi')
  call print_clock ('s_psi')
  call print_clock ('dipole')
  call print_clock ('quadrupole')
  write(stdout,*)
  write(stdout,*) '     General routines'   
  call print_clock( 'calbec' )
  call print_clock( 'fft' )
  call print_clock( 'ffts' )
  call print_clock( 'fftw' )
  call print_clock( 'interpolate' )
  call print_clock( 'davcio' )
  write( stdout, * )
#ifdef __PARA
  write( stdout,  * ) '     Parallel routines'
  call print_clock ('reduce')
  call print_clock( 'fft_scatter' )  
#endif

END SUBROUTINE print_clock_tddft



!-----------------------------------------------------------------------
! TDDFT setup
!-----------------------------------------------------------------------
SUBROUTINE tddft_setup
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, ionode
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, atm
  USE atom,          ONLY : rgrid
  USE wvfct,         ONLY : nbnd, et, wg, npwx, ecutwfc
  USE lsda_mod,      ONLY : nspin, lsda
  USE gvect,         ONLY : ngm, g
  USE gvecs,         ONLY : doublegrid
  USE fft_base,      ONLY : dfftp
  USE klist,         ONLY : xk, degauss, ngauss, nks, nelec
  USE constants,     ONLY : degspin, pi
  USE symm_base,     ONLY : nsym, s
  USE uspp_param,    ONLY : upf
  USE mp_global,     ONLY : inter_pool_comm 
  USE mp,            ONLY : mp_max, mp_min 
  USE ldaU,          ONLY : lda_plus_u
  USE scf,           ONLY : rho, rho_core, rhog_core, vltot, v, vrs, kedtau
  USE io_files,      ONLY : iunigk
  USE cell_base,     ONLY : tpiba2
  USE wvfct,         ONLY : npw, g2kin, igk
  USE dfunct,        ONLY : newd
  USE tddft_module   

  IMPLICIT none
  integer :: ik, ibnd
  real(dp) :: emin, emax
    
  call start_clock ('tddft_setup')
    
  ! initialize pseudopotentials
  call init_us_1
  call init_at_1

  ! computes the total local potential (external+scf) on the smooth grid
  call setlocal
  call set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
    
  ! compute the D for the pseudopotentials
  call newd
    
  ! computes the number of occupied bands for each k point
  allocate(nbnd_occ(nks))
  nbnd_occ(:) = 0
  do ik = 1, nks
     do ibnd = 1, nbnd
        if ( wg(ibnd,ik) > 1e-6 ) then
           nbnd_occ(ik) = ibnd
        end if
     end do
  end do
    
  ! computes alpha_pv
  emin = et (1, 1)
  do ik = 1, nks
    do ibnd = 1, nbnd
      emin = min (emin, et (ibnd, ik) )
    enddo
  enddo
#ifdef __PARA
  ! find the minimum across pools
  call mp_min( emin, inter_pool_comm )
#endif
  if (degauss.ne.0.0_dp) then
    call errore('tddft_setup', 'implemented only for insulators', -1)
  else
    emax = et (1, 1)
    do ik = 1, nks
      do ibnd = 1, nbnd
        emax = max (emax, et (ibnd, ik) )
      enddo
    enddo
#ifdef __PARA
    ! find the maximum across pools
    call mp_max( emax, inter_pool_comm )
#endif
    alpha_pv = 2.0_dp * (emax - emin)
  endif
  ! avoid zero value for alpha_pv
  alpha_pv = max (alpha_pv, 1.0d-2)

  rewind(iunigk)
  do ik = 1, nks
     call gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
     if (nks > 1) write(iunigk) igk
  end do
    
  call stop_clock('tddft_setup')
    
END SUBROUTINE tddft_setup

