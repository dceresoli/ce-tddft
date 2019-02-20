!
! Copyright (C) 2001-2010 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! TODO: nsave, restart_mode

!-----------------------------------------------------------------------
MODULE tddft_module
  !-----------------------------------------------------------------------
  !
  ! ... This module contains the variables used for TDDFT calculations
  !
  USE kinds,               ONLY : dp
  USE fixed_occ,           ONLY : tfixed_occ, f_inp ! occupations from input
  USE dynamics_module,     ONLY : dt                ! timestep 
  IMPLICIT NONE
  SAVE
  
  character(80) :: job             ! 'optical'
  integer  :: e_direction          ! impulse electric field direction: 1-x 2-y 3-z
  real(dp) :: e_strength           ! impulse electron field strength
  integer  :: nstep                ! number of timesteps for real-time tddft
  real(dp) :: conv_threshold       ! cg convergence threshold
  integer  :: nupdate_Dnm          ! update USPP Dnm matrix every n steps
  logical  :: l_circular_dichroism ! calculate circular dichroism
  logical  :: l_tddft_restart      ! restart propagation from the last step
  integer  :: iverbosity           ! verbosity level (default = 1)
  logical  :: molecule             ! use molecular routuines
  logical  :: ehrenfest            ! .true. if moving atoms

  complex(dp), parameter :: i_complex = (0.0_dp,1.0_dp)

  real(dp), allocatable :: r_pos(:,:)     ! position operator in real space
  real(dp), allocatable :: r_pos_s(:,:)   ! position operator in real space (smooth grid)
  integer, allocatable :: nbnd_occ(:)     ! occupied bands for each k-point
  integer :: nbnd_occ_max                 ! max number of occupied bands
  integer, parameter :: iuntdwfc = 51     ! to save TDDFT intermediate wfcs
  integer :: nwordtdwfc 
  integer, parameter :: iunevcn = 52      ! evc for restart
  real(dp) :: alpha_pv                    ! shift of conduction levels

  real(dp) :: max_seconds           ! max CPU time, in s
  integer :: isave_rho              ! output rho every .. steps

  ! wavepacket variables (wavepacket alwais moving in the -z direction)
  logical :: wavepacket = .false.   ! .true. to simulate a wavepacket
  real(dp) :: wp_pos(3)             ! wavepacket position in angstrom
  real(dp) :: wp_d(3)               ! wavepacket spread in ansgtrom
  real(dp) :: wp_ekin               ! wavepacket energy in eV
  integer :: wp_ibnd                ! (empty) band to accommodate wp

END MODULE tddft_module

