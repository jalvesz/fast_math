!
! SPDX-FileCopyrightText: 2023 Transvalor S.A.
!
! SPDX-License-Identifier: MIT
!
module fast_dotp
    !! A faster and more accurate implementation of the dot_product intrinsic. 
    !! It uses the same principle as fsum_chunk but considering local multiplications that can be vectorized for faster summation.
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    private
    
    public :: fprod
    integer, parameter :: chunk64 = 32
    integer, parameter :: chunk32 = 64
    
    interface fprod
        module procedure fprod_sp
        module procedure fprod_sp_weighted
        module procedure fprod_dp
        module procedure fprod_dp_weighted
    end interface
    
    contains

    pure function fprod_sp(a,b) result(p)
        integer, parameter :: wp = sp
        integer, parameter :: chunk = chunk32
        real(wp), intent(in) :: a(:)
        real(wp), intent(in) :: b(:)
        real(wp) :: p
        ! --
        real(wp) :: abatch(chunk)
        integer :: i, n, dr, rr
        ! -----------------------------
        n  = size(a)
        dr = n/chunk
        rr = n - dr*chunk
        
        abatch(:) = 0.0_wp
        do concurrent( i = 1:dr )
          abatch(1:chunk) = abatch(1:chunk) + a(chunk*i-chunk+1:chunk*i)*b(chunk*i-chunk+1:chunk*i)
        end do
        abatch(1:rr) = abatch(1:rr) + a(n-rr+1:n)*b(n-rr+1:n)
        
        p = 0.0_wp
        do concurrent( i = chunk/2:1:-1 )
          p = p + abatch(i)+abatch(chunk/2+i)
        end do
    end function
  
    pure function fprod_dp(a,b) result(p)
        integer, parameter :: wp = dp
        integer, parameter :: chunk = chunk64
        real(wp), intent(in) :: a(:)
        real(wp), intent(in) :: b(:)
        real(wp) :: p
        ! --
        real(wp) :: abatch(chunk)
        integer :: i, n, dr, rr
        ! -----------------------------
        n  = size(a)
        dr = n/chunk
        rr = n - dr*chunk
        
        abatch(:) = 0.0_wp
        do concurrent( i = 1:dr )
          abatch(1:chunk) = abatch(1:chunk) + a(chunk*i-chunk+1:chunk*i)*b(chunk*i-chunk+1:chunk*i)
        end do
        abatch(1:rr) = abatch(1:rr) + a(n-rr+1:n)*b(n-rr+1:n)
        
        p = 0.0_wp
        do concurrent( i = chunk/2:1:-1 )
          p = p + abatch(i)+abatch(chunk/2+i)
        end do
    end function
    
    pure function fprod_sp_weighted(a,b,w) result(p)
        integer, parameter :: wp = sp
        integer, parameter :: chunk = chunk32
        real(wp), intent(in) :: a(:)
        real(wp), intent(in) :: b(:)
        real(wp), intent(in) :: w(:)
        real(wp) :: p
        ! --
        real(wp) :: abatch(chunk)
        integer :: i, n, dr, rr
        ! -----------------------------
        n  = size(a)
        dr = n/chunk
        rr = n - dr*chunk
        
        abatch(:) = 0.0_wp
        do concurrent( i = 1:dr )
          abatch(1:chunk) = abatch(1:chunk) + a(chunk*i-chunk+1:chunk*i)*b(chunk*i-chunk+1:chunk*i)*w(chunk*i-chunk+1:chunk*i)
        end do
        abatch(1:rr) = abatch(1:rr) + a(n-rr+1:n)*b(n-rr+1:n)*w(n-rr+1:n)
        
        p = 0.0_wp
        do concurrent( i = chunk/2:1:-1 )
          p = p + abatch(i)+abatch(chunk/2+i)
        end do
    end function
  
    pure function fprod_dp_weighted(a,b,w) result(p)
        integer, parameter :: wp = dp
        integer, parameter :: chunk = chunk64
        real(wp), intent(in) :: a(:)
        real(wp), intent(in) :: b(:)
        real(wp), intent(in) :: w(:)
        real(wp) :: p
        ! --
        real(wp) :: abatch(chunk)
        integer :: i, n, dr, rr
        ! -----------------------------
        n  = size(a)
        dr = n/chunk
        rr = n - dr*chunk
        
        abatch(:) = 0.0_wp
        do concurrent( i = 1:dr )
          abatch(1:chunk) = abatch(1:chunk) + a(chunk*i-chunk+1:chunk*i)*b(chunk*i-chunk+1:chunk*i)*w(chunk*i-chunk+1:chunk*i)
        end do
        abatch(1:rr) = abatch(1:rr) + a(n-rr+1:n)*b(n-rr+1:n)*w(n-rr+1:n)
        
        p = 0.0_wp
        do concurrent( i = chunk/2:1:-1 )
          p = p + abatch(i)+abatch(chunk/2+i)
        end do
    end function
  end module fast_dotp