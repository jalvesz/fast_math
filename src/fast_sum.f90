!
! SPDX-FileCopyrightText: 2023 Transvalor S.A.
!
! SPDX-License-Identifier: MIT
!
module fast_sum
  !! Two fast & accurate sum are proposed for 1D arrays:
  !! By default, "fsum" will use the fsum_chunk approach. This method is at worst, one order of magnitud more accurate that "sum" and between 1.5 to 10 times faster
  !! A second approach is also proposed, "fsum_pair" which is the most accurate approach. cpu time can vary between x2 times slower or sometimes faster than intrinsic sum.
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  implicit none
  private

  public :: fsum, fsum_kahan
  integer, parameter :: chunk64 = 64
  integer, parameter :: chunk32 = 64
  
  interface fsum
    !! Source: to the best of knowledge: Alves J. but heavily inspired by this paper https://epubs.siam.org/doi/10.1137/19M1257780
      module procedure fsum_chunk_1d_sp
      module procedure fsum_chunk_1d_sp_mask
      module procedure fsum_chunk_1d_dp
      module procedure fsum_chunk_1d_dp_mask
  end interface

  interface fsum_kahan
      module procedure fsum_kahan_1d_sp
      module procedure fsum_kahan_1d_sp_mask
      module procedure fsum_kahan_1d_dp
      module procedure fsum_kahan_1d_dp_mask
  end interface

  interface vkahans
      module procedure vkahans_sp
      module procedure vkahans_dp
  end interface
  interface vkahans_m
      module procedure vkahans_m_sp
      module procedure vkahans_m_dp
  end interface

  contains

  pure function fsum_chunk_1d_sp(a) result(sout)
      integer, parameter :: wp = sp
      integer, parameter :: chunk = chunk32
      real(wp), intent(in) :: a(:)
      real(wp) :: sout
      ! --
      real(wp) :: abatch(chunk)
      integer :: i, n, r
      ! -----------------------------
      n  = size(a)
      r = mod(n,chunk)
      
      abatch(1:r)       = a(1:r)
      abatch(r+1:chunk) = 0._wp
      do i = r+1, n-r, chunk
       abatch(1:chunk) = abatch(1:chunk) + a(i:i+chunk-1)
      end do
      
      sout = 0.0_wp
      do i = 1, chunk/2
        sout = sout + abatch(i)+abatch(chunk/2+i)
      end do
  end function

  pure function fsum_chunk_1d_dp(a) result(sout)
      integer, parameter :: wp = dp
      integer, parameter :: chunk = chunk64
      real(wp), intent(in) :: a(:)
      real(wp) :: sout
      ! --
      real(wp) :: abatch(chunk)
      integer :: i, n, r
      ! -----------------------------
      n  = size(a)
      r = mod(n,chunk)
      
      abatch(1:r)       = a(1:r)
      abatch(r+1:chunk) = 0._wp
      do i = r+1, n-r, chunk
       abatch(1:chunk) = abatch(1:chunk) + a(i:i+chunk-1)
      end do

      sout = 0.0_wp
      do i = 1, chunk/2
        sout = sout + abatch(i)+abatch(chunk/2+i)
      end do
  end function
  
  pure function fsum_chunk_1d_sp_mask(a,mask) result(sout)
      integer, parameter :: wp = sp
      integer, parameter :: chunk = chunk32
      real(wp), intent(in) :: a(:)
      logical, intent(in) :: mask(:)
      real(wp) :: sout
      ! --
      real(wp) :: abatch(chunk)
      integer :: i, n, r
      ! -----------------------------
      n  = size(a)
      r = mod(n,chunk)
      
      abatch(1:r)       = merge( 0.0_wp , a(1:r) , mask(1:r) )
      abatch(r+1:chunk) = 0._wp
      do i = r+1, n-r, chunk
       abatch(1:chunk) = abatch(1:chunk) + merge( 0.0_wp , a(i:i+chunk-1), mask(i:i+chunk-1) )
      end do
      
      sout = 0.0_wp
      do i = 1, chunk/2
        sout = sout + abatch(i)+abatch(chunk/2+i)
      end do
  end function

  pure function fsum_chunk_1d_dp_mask(a,mask) result(sout)
      integer, parameter :: wp = dp
      integer, parameter :: chunk = chunk64
      real(wp), intent(in) :: a(:)
      logical, intent(in) :: mask(:)
      real(wp) :: sout
      ! --
      real(wp) :: abatch(chunk)
      integer :: i, n, r
      ! -----------------------------
      n  = size(a)
      r = mod(n,chunk)
      
      abatch(1:r)       = merge( 0.0_wp , a(1:r) , mask(1:r) )
      abatch(r+1:chunk) = 0._wp
      do i = r+1, n-r, chunk
       abatch(1:chunk) = abatch(1:chunk) + merge( 0.0_wp , a(i:i+chunk-1), mask(i:i+chunk-1) )
      end do

      sout = 0.0_wp
      do i = 1, chunk/2
        sout = sout + abatch(i)+abatch(chunk/2+i)
      end do
  end function

  pure function fsum_kahan_1d_sp(a) result(sout)
      integer, parameter :: wp = sp
      integer, parameter :: chunk = chunk32
      real(wp), intent(in) :: a(:)
      real(wp) :: sout
      ! --
      real(wp) :: sbatch(chunk)
      real(wp) :: cbatch(chunk)
      integer :: i, n, r
      ! -----------------------------
      n  = size(a)
      r = mod(n,chunk)

      sbatch(1:r) = a(1:r)
      sbatch(r+1:chunk)  = 0.0_wp
      cbatch = 0.0_wp
      do i = r+1, n-r, chunk
        call vkahans( a(i:i+chunk-1) , sbatch(1:chunk) , cbatch(1:chunk) )
      end do  
      
      sout = 0.0_wp
      do i = 1,chunk
        call vkahans( sbatch(i) , sout , cbatch(i) )
      end do
  end function

  pure function fsum_kahan_1d_dp(a) result(sout)
      integer, parameter :: wp = dp
      integer, parameter :: chunk = chunk64
      real(wp), intent(in) :: a(:)
      real(wp) :: sout
      ! --
      real(wp) :: sbatch(chunk)
      real(wp) :: cbatch(chunk)
      integer :: i, n, r
      ! -----------------------------
      n  = size(a)
      r = mod(n,chunk)

      sbatch(1:r) = a(1:r)
      sbatch(r+1:chunk)  = 0.0_wp
      cbatch = 0.0_wp
      do i = r+1, n-r, chunk
        call vkahans( a(i:i+chunk-1) , sbatch(1:chunk) , cbatch(1:chunk) )
      end do    
      
      sout = 0.0_wp
      do i = 1,chunk
        call vkahans( sbatch(i) , sout , cbatch(i) )
      end do
  end function

  pure function fsum_kahan_1d_sp_mask(a,mask) result(sout)
      integer, parameter :: wp = sp
      integer, parameter :: chunk = chunk32
      real(wp), intent(in) :: a(:)
      logical, intent(in) :: mask(:)
      real(wp) :: sout
      ! --
      real(wp) :: sbatch(chunk)
      real(wp) :: cbatch(chunk)
      integer :: i, n, r
      ! -----------------------------
      n  = size(a)
      r = mod(n,chunk)

      sbatch(1:r) = merge( 0.0_wp , a(1:r) , mask(1:r) )
      sbatch(r+1:chunk)  = 0.0_wp
      cbatch = 0.0_wp
      do i = r+1, n-r, chunk
        call vkahans_m( a(i:i+chunk-1) , sbatch(1:chunk) , cbatch(1:chunk), mask(i:i+chunk-1) )
      end do 

      sout = 0.0_wp
      do i = 1,chunk
        call vkahans( sbatch(i) , sout , cbatch(i) )
      end do
  end function

  pure function fsum_kahan_1d_dp_mask(a,mask) result(sout)
      integer, parameter :: wp = dp
      integer, parameter :: chunk = chunk64
      real(wp), intent(in) :: a(:)
      logical, intent(in) :: mask(:)
      real(wp) :: sout
      ! --
      real(wp) :: sbatch(chunk)
      real(wp) :: cbatch(chunk)
      integer :: i, n, r
      ! -----------------------------
      n  = size(a)
      r = mod(n,chunk)

      sbatch(1:r) = merge( 0.0_wp , a(1:r) , mask(1:r) )
      sbatch(r+1:chunk)  = 0.0_wp
      cbatch = 0.0_wp
      do i = r+1, n-r, chunk
        call vkahans_m( a(i:i+chunk-1) , sbatch(1:chunk) , cbatch(1:chunk), mask(i:i+chunk-1) )
      end do 

      sout = 0.0_wp
      do i = 1,chunk
        call vkahans( sbatch(i) , sout , cbatch(i) )
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

  elemental subroutine vkahans_m_sp(a,s,c,m)
  integer, parameter :: wp = sp
  include 'utilities/vkahans_m.inc'
  end subroutine  

  elemental subroutine vkahans_m_dp(a,s,c,m)
  integer, parameter :: wp = dp
  include 'utilities/vkahans_m.inc'
  end subroutine  

end module fast_sum