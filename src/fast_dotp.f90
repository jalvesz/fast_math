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
    
    public :: fprod, fprod_kahan
    integer, parameter :: chunk64 = 64
    integer, parameter :: chunk32 = 64
    
    interface fprod
        module procedure fprod_sp
        module procedure fprod_sp_weighted
        module procedure fprod_dp
        module procedure fprod_dp_weighted
    end interface

    interface fprod_kahan
        module procedure fprod_kahan_sp
        module procedure fprod_kahan_sp_weighted
        module procedure fprod_kahan_dp
        module procedure fprod_kahan_dp_weighted
    end interface

    interface vkahans
      module procedure vkahans_sp
      module procedure vkahans_dp
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
        integer :: i, n, r
        ! -----------------------------
        n = size(a)
        r = mod(n,chunk)

        abatch(1:r)       = a(1:r)*b(1:r)
        abatch(r+1:chunk) = 0._wp
        do i = r+1, n-r, chunk
         abatch(1:chunk) = abatch(1:chunk) + a(i:i+chunk-1)*b(i:i+chunk-1)
        end do
        
        p = 0.0_wp
        do i = 1, chunk/2
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
        integer :: i, n, r
        ! -----------------------------
        n = size(a)
        r = mod(n,chunk)

        abatch(1:r)       = a(1:r)*b(1:r)
        abatch(r+1:chunk) = 0._wp
        do i = r+1, n-r, chunk
         abatch(1:chunk) = abatch(1:chunk) + a(i:i+chunk-1)*b(i:i+chunk-1)
        end do
        
        p = 0.0_wp
        do i = 1, chunk/2
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
        integer :: i, n, r
        ! -----------------------------
        n = size(a)
        r = mod(n,chunk)

        abatch(1:r)       = a(1:r)*b(1:r)*w(1:r)
        abatch(r+1:chunk) = 0._wp
        do i = r+1, n-r, chunk
         abatch(1:chunk) = abatch(1:chunk) + a(i:i+chunk-1)*b(i:i+chunk-1)*w(i:i+chunk-1)
        end do
        
        p = 0.0_wp
        do i = 1, chunk/2
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
        integer :: i, n, r
        ! -----------------------------
        n = size(a)
        r = mod(n,chunk)

        abatch(1:r)       = a(1:r)*b(1:r)*w(1:r)
        abatch(r+1:chunk) = 0._wp
        do i = r+1, n-r, chunk
         abatch(1:chunk) = abatch(1:chunk) + a(i:i+chunk-1)*b(i:i+chunk-1)*w(i:i+chunk-1)
        end do
        
        p = 0.0_wp
        do i = 1, chunk/2
          p = p + abatch(i)+abatch(chunk/2+i)
        end do
    end function

    pure function fprod_kahan_sp(a,b) result(p)
        integer, parameter :: wp = sp
        integer, parameter :: chunk = chunk32
        real(wp), intent(in) :: a(:)
        real(wp), intent(in) :: b(:)
        real(wp) :: p
        ! --
        real(wp) :: sbatch(chunk)
        real(wp) :: cbatch(chunk)
        integer :: i, n, r
        ! -----------------------------
        n = size(a)
        r = mod(n,chunk)
        
        sbatch(1:r) = a(1:r) * b(1:r)
        sbatch(r+1:chunk)  = 0.0_wp
        cbatch = 0.0_wp 
        do i = r+1, n-r, chunk
          call vkahans( a(i:i+chunk-1) * b(i:i+chunk-1) , sbatch(1:chunk) , cbatch(1:chunk) )
        end do
        
        p = 0.0_wp
        do i = 1,chunk
          call vkahans( sbatch(i) , p , cbatch(i) )
        end do
    end function

    pure function fprod_kahan_dp(a,b) result(p)
        integer, parameter :: wp = dp
        integer, parameter :: chunk = chunk64
        real(wp), intent(in) :: a(:)
        real(wp), intent(in) :: b(:)
        real(wp) :: p
        ! --
        real(wp) :: sbatch(chunk)
        real(wp) :: cbatch(chunk)
        integer :: i, n, r
        ! -----------------------------
        n = size(a)
        r = mod(n,chunk)
        
        sbatch(1:r) = a(1:r) * b(1:r)
        sbatch(r+1:chunk)  = 0.0_wp
        cbatch = 0.0_wp 
        do i = r+1, n-r, chunk
          call vkahans( a(i:i+chunk-1) * b(i:i+chunk-1) , sbatch(1:chunk) , cbatch(1:chunk) )
        end do
        
        p = 0.0_wp
        do i = 1,chunk
          call vkahans( sbatch(i) , p , cbatch(i) )
        end do
    end function

    pure function fprod_kahan_sp_weighted(a,b,w) result(p)
        integer, parameter :: wp = sp
        integer, parameter :: chunk = chunk32
        real(wp), intent(in) :: a(:)
        real(wp), intent(in) :: b(:)
        real(wp), intent(in) :: w(:)
        real(wp) :: p
        ! --
        real(wp) :: sbatch(chunk)
        real(wp) :: cbatch(chunk)
        integer :: i, n, r
        ! -----------------------------
        n = size(a)
        r = mod(n,chunk)
        
        sbatch(1:r) = a(1:r) * b(1:r) * w(1:r)
        sbatch(r+1:chunk)  = 0.0_wp
        cbatch = 0.0_wp
        do i = r+1, n-r, chunk
          call vkahans( a(i:i+chunk-1) * b(i:i+chunk-1) * w(i:i+chunk-1) , sbatch(1:chunk) , cbatch(1:chunk) )
        end do
        
        p = 0.0_wp
        do i = 1,chunk
          call vkahans( sbatch(i) , p , cbatch(i) )
        end do
    end function

    pure function fprod_kahan_dp_weighted(a,b,w) result(p)
        integer, parameter :: wp = dp
        integer, parameter :: chunk = chunk64
        real(wp), intent(in) :: a(:)
        real(wp), intent(in) :: b(:)
        real(wp), intent(in) :: w(:)
        real(wp) :: p
        ! --
        real(wp) :: sbatch(chunk)
        real(wp) :: cbatch(chunk)
        integer :: i, n, r
        ! -----------------------------
        n = size(a)
        r = mod(n,chunk)
        
        sbatch(1:r) = a(1:r) * b(1:r) * w(1:r)
        sbatch(r+1:chunk)  = 0.0_wp
        cbatch = 0.0_wp 
        do i = r+1, n-r, chunk
          call vkahans( a(i:i+chunk-1) * b(i:i+chunk-1) * w(i:i+chunk-1) , sbatch(1:chunk) , cbatch(1:chunk) )
        end do
        
        p = 0.0_wp
        do i = 1,chunk
          call vkahans( sbatch(i) , p , cbatch(i) )
        end do
    end function

    elemental subroutine vkahans_sp(a,s,c)
    integer, parameter :: wp = sp
    include 'utilities/vkahans.inc'
    end subroutine  

    elemental subroutine vkahans_dp(a,s,c)
    integer, parameter :: wp = dp
    include 'utilities/vkahans.inc'
    end subroutine  

  end module fast_dotp