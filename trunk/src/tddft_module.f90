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
  USE kinds, ONLY : DP
  
  IMPLICIT NONE
  SAVE
  
  character(80) :: job             ! 'optical'
  integer  :: e_direction          ! impulse electric field direction: 1-x 2-y 3-z
  real(dp) :: e_strength           ! impulse electron field strength
  real(dp) :: dt                   ! timestep
  integer  :: nstep                ! number of timesteps for real-time tddft
  real(dp) :: conv_threshold       ! cg convergence threshold
  integer  :: nupdate_Dnm          ! update USPP Dnm matrix every n steps
  logical  :: l_circular_dichroism ! calculate circular dichroism
  logical  :: l_tddft_restart      ! restart propagation from the last step
  integer  :: iverbosity           ! verbosity level (default = 1)
  logical  :: molecule             ! use molecular routuines

  complex(dp), parameter :: i_complex = (0.0_dp,1.0_dp)

  real(dp), allocatable :: r_pos(:,:)     ! position operator in real space
  real(dp), allocatable :: r_pos_s(:,:)   ! position operator in real space (smooth grid)
  integer, allocatable :: nbnd_occ(:)     ! occupied bands for each k-point
  integer :: nbnd_occ_max                 ! max number of occupied bands
  integer, parameter :: iuntdwfc = 51     ! to save TDDFT intermediate wfcs
  integer :: nwordtdwfc 
  integer, parameter :: iunevcn = 52      ! evc for restart
  real(dp) :: alpha_pv                    ! shift of conduction levels

  integer :: tddft_exit_code = 0

END MODULE tddft_module

