!
! Copyright (C) 2001-2018 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
PROGRAM tddft_main
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the real time TDDFT propagation.
  ! ... Authors: Xiaofeng Qian and Davide Ceresoli
  ! ...
  ! ... References:
  ! ...   Xiaofeng Qian, Ju Li, Xi Lin, and Sidney Yip, PRB 73, 035408 (2006)
  ! ...
  USE kinds,           ONLY : DP
  USE io_files,        ONLY : tmp_dir !, create_directory
  USE io_global,       ONLY : stdout, meta_ionode, meta_ionode_id
  USE mp,              ONLY : mp_bcast
  USE cell_base,       ONLY : tpiba
  USE tddft_module,    ONLY : job, molecule, max_seconds, tddft_exit_code
  USE check_stop,      ONLY : check_stop_init
  USE control_flags,   ONLY : io_level, gamma_only, use_para_diag
  USE mp_global,       ONLY : mp_startup, nproc_pool_file
  USE mp_bands,        ONLY : nbgrp
  USE mp_pools,        ONLY : nproc_pool
  USE mp_world,        ONLY : world_comm
  USE environment,     ONLY : environment_start, environment_end
  USE lsda_mod,        ONLY : nspin
  USE wvfct,           ONLY : nbnd
  USE uspp,            ONLY : okvan
  USE wvfct,           ONLY : nbnd
  USE io_global,       ONLY : stdout
  USE noncollin_module,ONLY : noncolin
  ! for pluginization
  USE input_parameters, ONLY : nat_ => nat, ntyp_ => ntyp
  USE input_parameters, ONLY : assume_isolated_ => assume_isolated, &
                               ibrav_ => ibrav
  USE ions_base,        ONLY : nat, ntyp => nsp
  USE cell_base,        ONLY : ibrav
  USE tddft_version
  USE iotk_module  
  !------------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER (LEN=9)   :: code = 'TDDFT'
  CHARACTER (LEN=10)  :: dirname = 'dummy'
  LOGICAL, EXTERNAL  :: check_para_diag
  !------------------------------------------------------------------------

  ! begin with the initialization part
#ifdef __MPI
  call mp_startup(start_images=.true.)
#else
  call mp_startup(start_images=.false.)
#endif
  call environment_start (code)

  ! read plugin command line arguments, if any
  if (meta_ionode) call plugin_arguments()
  call plugin_arguments_bcast( meta_ionode_id, world_comm )

#ifndef __BANDS
  if (nbgrp > 1) &
    call errore('tddft_main', 'configure and recompile TDDFT with --enable-band-parallel', 1)
#endif

  write(stdout,*)
  write(stdout,'(5X,''***** This is TDDFT git revision '',A,'' *****'')') tddft_git_revision
  write(stdout,'(5X,''***** you can cite: X. Qian et al. Phys. Rev. B 73, 035408 (2006)         *****'')')
  write(stdout,'(5X,''***** in publications or presentations arising from this work.            *****'')')
  write(stdout,*)

  call tddft_readin()
  call check_stop_init( max_seconds )

  io_level = 1
 
  ! read ground state wavefunctions
  call read_file
#ifdef __MPI
  use_para_diag = check_para_diag(nbnd)
#else
  use_para_diag = .false.
#endif

  call tddft_openfil

  if (gamma_only) call errore ('tdddft_main', 'Cannot run TDFFT with gamma_only == .true. ', 1)
#ifdef __BANDS
  if (nbgrp > 1 .and. (twfcollect .eqv. .false.)) &
    call errore('tddft_main', 'Cannot use band-parallelization without wf_collect in SCF', 1)
#endif
  if (noncolin) call errore('tdddft_main', 'non-collinear not supported yet', 1)

  nat_ = nat
  ntyp_ = ntyp
  ibrav_ = ibrav
  assume_isolated_ = 'none'
  call plugin_read_input()
  call tddft_allocate()
  call tddft_setup()
  call tddft_summary()

#ifdef __BANDS
  call init_parallel_over_band(inter_bgrp_comm, nbnd)
#endif

  ! calculation
  select case (trim(job))
  case ('optical')
     if (molecule) then
        call molecule_optical_absorption
     else
        call errore('tddft_main', 'solids are not yet implemented', 1)
     endif

  case default
     call errore('tddft_main', 'wrong or undefined job in input', 1)

  end select
  
  ! print timings and stop the code
  call tddft_closefil
  call print_clock_tddft
  call environment_end(code)
  call stop_code( .true. )
  
  STOP
  
END PROGRAM tddft_main

