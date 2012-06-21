!
! Copyright (C) 2001-2010 Quantum-ESPRESSO group
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
  USE control_flags,   ONLY : io_level, gamma_only, use_para_diag
  USE mp_global,       ONLY : mp_startup
  USE environment,     ONLY : environment_start
  USE klist,           ONLY : nelec
  USE tddft_module
  !------------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER (LEN=9)   :: code = 'TDDFT'
  !------------------------------------------------------------------------

#ifdef __PARA
  call mp_startup()
#endif

  call environment_start(code)

  ! ... Intel compilers v .ge.8 allocate a lot of stack space
  ! ... Stack limit is often small, thus causing SIGSEGV and crash
  CALL remove_stack_limit()

  CALL tddft_readin()
  io_level = 1
#ifdef __PARA
  IF ( use_para_diag )  CALL check_para_diag( nelec )
#else
  use_para_diag = .FALSE.
#endif
  ! read ground state wavefunctions
  call read_file
  call tddft_openfil
  
  if ( gamma_only ) call errore('tddft_main', 'gamma_only == .true.', 1)
  
  call tddft_allocate()
  call tddft_setup()
  call tddft_summary()
  call print_clock('TDDFT')
  
  ! calculation
  select case (trim(job))
  case ('optical')
     call tddft_optical_absorption
  case default
     call errore('tddft_main', 'wrong or undefined job in input', 1)
  end select
  
  ! print timings and stop the code
  call print_clock_tddft()
  call stop_code(.true.)
  
  STOP
  
END PROGRAM tddft_main



!----------------------------------------------------------------------------
SUBROUTINE stop_code( flag )
  !----------------------------------------------------------------------------
  !
  ! ... Synchronize processes before stopping.
  !
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_end, mp_barrier
  USE parallel_include
  IMPLICIT NONE
  LOGICAL :: flag

  call mp_barrier()
  call mp_end()

#if defined (__t3e)
  call set_d_stream( 0 )
#endif

  if ( flag ) then
     stop
  else
     stop 1
  endif

END SUBROUTINE stop_code
