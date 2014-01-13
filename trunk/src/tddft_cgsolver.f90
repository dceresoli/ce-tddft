!
! Copyright (C) 2001-2014 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Fortran implementation of the Conjugate Gradient Square solver
! author: Xiaofeng Qian, MIT (2008)
!

!-----------------------------------------------------------------------
MODULE tddft_cgsolver_module
  !-----------------------------------------------------------------------
  !
  ! ... this module contains variables and temporaries for the CGS solver
  !
  USE kinds, ONLY : dp
  IMPLICIT NONE
  SAVE

  integer     :: flag
  real(dp)    :: tolb, normr, relres, n2b, normrmin
  complex(dp) :: rho, rho1, alpha, beta, rtvh
  complex(dp), allocatable :: r(:), Ax(:), rt(:), vh(:), u(:), &
                              uh(:), q(:), qh(:), p(:)

END MODULE tddft_cgsolver_module


!-----------------------------------------------------------------------
SUBROUTINE tddft_cgsolver_initialize(ndmx, nbnd)
  !-----------------------------------------------------------------------
  !
  ! ... allocate memory for the solver
  !
  USE tddft_cgsolver_module
  implicit none
  integer, intent(in) :: ndmx
  integer, intent(in) :: nbnd
    
  allocate (r(ndmx), rt(ndmx), Ax(ndmx), u(ndmx), p(ndmx), q(ndmx), &
            qh(ndmx), uh(ndmx), vh(ndmx) )
    
END SUBROUTINE tddft_cgsolver_initialize
  
  
!-----------------------------------------------------------------------
SUBROUTINE tddft_cgsolver_finalize()
  !-----------------------------------------------------------------------
  !
  ! ... deallocate memory
  !
  USE tddft_cgsolver_module
  implicit none
  deallocate( r, rt, Ax, u, p, q, qh, uh, vh )
    
END SUBROUTINE tddft_cgsolver_finalize
  

!----------------------------------------------------------------------
SUBROUTINE tddft_cgsolver (A, b, x, ndmx, ndim, tol, ik, iter, flag_global,  &
                           anorm, nbnd, ee)
  !----------------------------------------------------------------------
  !
  ! ... Conjugate-Gradient Square method for solving:   A * x = b
  ! ... where: A*x is evaluated by subroutine 'A', and 'A' is implicit
  ! ... general square-matrix.
  !                                            Xiaofeng Qian, MIT (2008)
  USE kinds,     ONLY : dp
  USE mp_pools,  ONLY : intra_pool_comm
  USE mp,        ONLY : mp_sum
  USE tddft_cgsolver_module
  !----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: ndmx               ! the maximum dimension of the vectors
  integer, intent(in) :: ndim               ! the actual dimension of the vectors
  integer, intent(in) :: nbnd               ! the number of bands
  integer, intent(in) :: ik                 ! the k point
  integer, intent(out) :: iter              ! counter on iterations
  integer, intent(out) :: flag_global
  real(dp), intent(in) :: tol               ! threshold for convergence
  real(dp), intent(out) :: anorm            ! the norm of the error in the solution
  complex(dp), intent(in) :: ee             ! i*dt/2
  complex(dp), intent(in) :: b(ndmx,nbnd)   ! input: the known term
  complex(dp), intent(out) :: x(ndmx,nbnd)  ! the solution of the linear system
  external A                                ! the subroutine computing A*x
  !----------------------------------------------------------------------
  integer, parameter :: maxit = 200          ! the maximum number of iterations
  complex(dp), external :: zdotc
  real(dp), external    :: ddot
  integer :: imin, stag, i, ibnd

  if (.not. allocated(r)) call errore('tddft_cgsolver', 'cgsolver not initialized', 1)
  
  call start_clock ('cgsolver')

  ! initialize module variables
  tolb = 0.d0
  n2b  = 0.d0
  relres = 0.d0
  normr= 0.d0
  rho  = 0.d0
  rho1 = 0.d0
  r    = (0.d0, 0.d0)
  Ax   = (0.d0, 0.d0)
  
  flag           = 1
  imin           = 0
  
  !----------------------------------------------------------------------
  ! loop over bands
  !----------------------------------------------------------------------
  do ibnd = 1, nbnd

     n2b = dble(zdotc(ndim, b(1,ibnd), 1, b(1,ibnd), 1))
#ifdef __MPI
     call mp_sum(n2b, intra_pool_comm)
#endif
     n2b = dsqrt(n2b)
     tolb = tol * n2b

     call A(ndim, x(1,ibnd), Ax, ee, ik, 1)
     
     r = b(:,ibnd) - Ax
     normr = dble(zdotc( ndim, r, 1, r, 1))
#ifdef __MPI
     call mp_sum(normr, intra_pool_comm)
#endif  
     normr = dsqrt(normr)
     
     if (normr < tolb) then
        flag = 0
        cycle
     endif
     
     rt = r
     normrmin = normr
     stag = 0
     rho  = cmplx(1.d0, 0.d0)
     
     ! CG iteration
     do i = 1, maxit
        rho1 = rho
        rho = zdotc(ndim, rt, 1, r, 1)
#ifdef __MPI
        call mp_sum(rho, intra_pool_comm)
#endif
        
        if (rho == (0.d0, 0.d0)) then  ! TODO: safe FP comparison
           flag  = 2
           return
        endif

        if (i == 1) then
           u = r
           p = u
        else
           beta = rho / rho1
           if ( beta == (0.d0, 0.d0)) then  ! TODO: safe FP comparison
              flag  = 3
              return
           endif
           u = r + beta * q
           p = u + beta * (q + beta * p)
        endif
        
        call A(ndim, p, vh, ee, ik, 1)
        
        rtvh = zdotc(ndim, rt, 1, vh, 1)
#ifdef __MPI
        call mp_sum(rtvh, intra_pool_comm)
#endif
        
        if (rtvh == (0.d0, 0.d0)) then  ! TODO: safe FP comparison
           flag = 4
           return
        else
           alpha = rho / rtvh
        endif
        
        if (alpha == (0.d0, 0.d0)) then  ! TODO: safe FP comparison
           flag = 5
           stag = 1
           return
        endif
        
        q  = u - alpha * vh
        
        uh = u + q
        
        x(:,ibnd)  = x(:,ibnd) + alpha * uh
        
        call A(ndim, x(1,ibnd), Ax, ee, ik, 1)
        
        normr = sum(conjg(b(:,ibnd) - Ax) * (b(:,ibnd) - Ax))
#ifdef __MPI
        call mp_sum(normr, intra_pool_comm)
#endif
        
        if (normr <= tolb) then
           flag = 0
           iter = i
           exit
        endif
        
        if (stag == 1) then
           flag  = 5
           return
        endif
        
        if (normr < normrmin)  then
           normrmin = normr
           imin = i
        endif
        
        call A(ndim, uh, qh, ee, ik, 1)
        
        r  = r - alpha * qh
        
     enddo ! i

  end do ! ibnd
  !----------------------------------------------------------------------
  ! end of the loop over bands
  !----------------------------------------------------------------------
  
  if (flag > 0) then
     call errore('tddft_cgsolver', 'cgsolver cannot achieve convergence', flag)
     stop
  end if
  call stop_clock ('cgsolver')
  
  return
  
END SUBROUTINE tddft_cgsolver


