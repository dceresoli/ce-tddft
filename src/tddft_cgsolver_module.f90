!
! Copyright (C) 2001-2010 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

module tddft_cgsolver_module
  
  USE kinds, ONLY : DP
  
  integer     :: flag
  real(dp)    :: tolb, normr, relres, n2b, normrmin
  complex(dp) :: rho, rho1, alpha, beta, rtvh
  complex(dp), allocatable  :: r(:), Ax(:), rt(:), & 
       vh(:), u(:),  uh(:), q(:), qh(:), p(:)
  
contains
  
  subroutine tddft_cgsolver_initialize(ndmx, nbnd)
    integer :: ndmx
    integer :: nbnd
    
    allocate( r(ndmx) )
    allocate( rt(ndmx) )
    allocate( Ax(ndmx) )
    allocate( u(ndmx) )
    allocate( p(ndmx) )
    allocate( q(ndmx) )
    allocate( qh(ndmx) )
    allocate( uh(ndmx) )
    allocate( vh(ndmx) )
    
  end subroutine tddft_cgsolver_initialize
  
  
  
  subroutine tddft_cgsolver_finalize()
    
    deallocate( r )
    deallocate( rt )
    deallocate( Ax )
    deallocate( u )
    deallocate( p )
    deallocate( q )
    deallocate( qh )
    deallocate( uh )
    deallocate( vh )
    
  end subroutine tddft_cgsolver_finalize
  
!----------------------------------------------------------------------
subroutine tddft_cgsolver (A, b, x, &
     ndmx, ndim, tol, ik, iter, flag_global, anorm, nbnd, ee)
  !----------------------------------------------------------------------
  !  conjugate-gradient square method for solving:   A * x = b
  !  where: A*x is evaluated by subroutine 'A', 
  !         and 'A' is implicit general square-matrix.
  !                                            Xiaofeng Qian, MIT (2008)
  !----------------------------------------------------------------------
  use kinds, only : dp
  use mp_global,  only : intra_pool_comm
  use mp,         only : mp_sum
  
  implicit none
  !
  !   first the I/O variables
  !
  integer :: ndmx, & ! input: the maximum dimension of the vectors
             ndim, & ! input: the actual dimension of the vectors
             iter, & ! output: counter on iterations
             nbnd, & ! input: the number of bands
             ik      ! input: the k point
  
  real(dp) :: tol
  real(dp) :: anorm   ! output: the norm of the error in the solution

  
  complex(dp) :: ee  ! input: i*dt/2
  complex(dp) :: &
             x (ndmx, nbnd), & ! output: the solution of the linear system
             b (ndmx, nbnd)    ! input: the known term

  external A    ! input: the subroutine computing A*x
  !
  complex(dp), external :: ZDOTC
  real(dp), external    :: DDOT
  !
  !  here the local variables
  !
  integer, parameter :: maxit = 200          ! the maximum number of iterations
  
  ! TDDFT internal variables
  integer :: flag_global, imin, stag, i, ibnd
  
  call start_clock ('cgsolver')

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
  
  do ibnd = 1, nbnd
     !
     n2b = dble( ZDOTC (ndim, b(1,ibnd), 1, b(1,ibnd), 1) )
#if defined (__PARA)
     call mp_sum( n2b, intra_pool_comm )
#endif
     n2b = dsqrt( n2b )
     tolb = tol * n2b
     
     call A(ndim, x(1, ibnd), Ax, ee, ik, 1)
     
     r = b(:, ibnd) - Ax
     normr = dble(ZDOTC( ndim, r, 1, r, 1))
#if defined (__PARA)
     call mp_sum( normr, intra_pool_comm )
#endif  
     normr = dsqrt(normr)
     
     if ( normr < tolb ) then
        flag = 0
        cycle
     end if
     
     rt = r
     normrmin = normr
     stag = 0
     rho  = cmplx(1.d0, 0.d0)
     
     do i = 1, maxit
        
        rho1 = rho
        rho  = ZDOTC( ndim, rt, 1, r, 1)
#if defined (__PARA)
        call mp_sum( rho, intra_pool_comm )
#endif
        
        if ( rho == (0.d0, 0.d0) ) then
           flag  = 2
           return
        end if
        if (i == 1) then
           u = r
           p = u
        else
           beta = rho / rho1
           if ( beta == (0.d0, 0.d0) ) then
              flag  = 3
              return
           end if
           u = r + beta * q
           p = u + beta * (q + beta * p)
        end if
        
        call A(ndim, p, vh, ee, ik, 1)
        
        rtvh = ZDOTC(ndim, rt, 1, vh, 1 )
#if defined (__PARA)
        call mp_sum( rtvh, intra_pool_comm )
#endif
        
        if (rtvh == (0.d0, 0.d0) ) then
           flag = 4
           return
        else
           alpha = rho / rtvh
        end if
        
        if (alpha == (0.d0, 0.d0) ) then
           flag = 5
           stag = 1
           return
        end if
        
        q  = u - alpha * vh
        
        uh = u + q
        
        x(:,ibnd)  = x(:,ibnd) + alpha * uh
        
        call A(ndim, x(1, ibnd), Ax, ee, ik, 1)
        
        normr = sum( conjg( b(:,ibnd) - Ax ) * ( b(:,ibnd) - Ax ) )
#if defined (__PARA)
        call mp_sum( normr, intra_pool_comm )
#endif
        
        if (normr <= tolb) then
           flag = 0
           iter = i
           exit
        end if
        
        if (stag == 1) then
           flag  = 5
           return
        end if
        
        if (normr < normrmin)  then
           normrmin = normr
           imin = i
        end if
        
        call A(ndim, uh, qh, ee, ik, 1)
        
        r  = r - alpha * qh
        
     end do                    ! nbnd
     
  end do                       ! maxit
  
  if (flag > 0) then
     call errore('tddft_cgsolver', 'cgsolver cannot achieve convergence', flag)
     stop
  end if
  call stop_clock ('cgsolver')
  
  return
  
end subroutine tddft_cgsolver


end module tddft_cgsolver_module
