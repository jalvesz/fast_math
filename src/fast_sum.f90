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

  public :: fsum, fsum_pair
  integer, parameter :: chunk64 = 32
  integer, parameter :: chunk32 = 64
  
  interface fsum
    !! Source: to the best of knowledge: Alves J. but heavily inspired by this paper https://epubs.siam.org/doi/10.1137/19M1257780
      module procedure fsum_chunk_1d_sp
      module procedure fsum_chunk_1d_sp_mask
      module procedure fsum_chunk_1d_dp
      module procedure fsum_chunk_1d_dp_mask
  end interface

  interface fsum_pair
    !! Source: https://fortran-lang.discourse.group/t/some-intrinsic-sums/5760/26
      module procedure fsum_pair_1d_sp
      module procedure fsum_pair_1d_sp_mask
      module procedure fsum_pair_1d_dp
      module procedure fsum_pair_1d_dp_mask
  end interface
  
  contains

  pure recursive function fsum_pair_1d_sp(a) result(bout)
      integer, parameter :: wp=sp
      integer, parameter :: chunk = chunk32
      real(wp), intent(in) :: a(:)
      real(wp) :: bout
      integer :: n, i
      ! -----------------------------
      n = size(a)
      if (n.le.chunk) then
        bout = 0.0_wp
        do concurrent(i = 1:n)
           bout = bout + a(i)
        end do
      else
         bout = fsum_pair(a(1:n/2)) + fsum_pair(a(n/2+1:n))
      end if
  end function

  pure recursive function fsum_pair_1d_dp(a) result(bout)
      integer, parameter :: wp=dp
      integer, parameter :: chunk = chunk64
      real(wp), intent(in) :: a(:)
      real(wp) :: bout
      integer :: n, i
      ! -----------------------------
      n = size(a)
      if (n.le.chunk) then
        bout = 0.0_wp
        do concurrent(i = 1:n)
           bout = bout + a(i)
        end do
      else
         bout = fsum_pair(a(1:n/2)) + fsum_pair(a(n/2+1:n))
      end if
  end function
  
  pure recursive function fsum_pair_1d_sp_mask(a,mask) result(bout)
      integer, parameter :: wp=sp
      integer, parameter :: chunk = chunk32
      real(wp), intent(in) :: a(:)
      logical, intent(in) :: mask(:)
      real(wp) :: bout
      integer :: n, i
      ! -----------------------------
      n = size(a)
      if (n.le.chunk) then
         bout = 0.0_wp
         do concurrent(i = 1:n)
            if(mask(i)) bout = bout + a(i)
         end do
      else
         bout = fsum_pair(a(1:n/2),mask(1:n/2)) + fsum_pair(a(n/2+1:n),mask(n/2+1:n))
      end if
  end function

  pure recursive function fsum_pair_1d_dp_mask(a,mask) result(bout)
      integer, parameter :: wp=dp
      integer, parameter :: chunk = chunk64
      real(wp), intent(in) :: a(:)
      logical, intent(in) :: mask(:)
      real(wp) :: bout
      integer :: n, i
      ! -----------------------------
      n = size(a)
      if (n.le.chunk) then
         bout = 0.0_wp
         do concurrent(i = 1:n)
            if(mask(i)) bout = bout + a(i)
         end do
      else
         bout = fsum_pair(a(1:n/2),mask(1:n/2)) + fsum_pair(a(n/2+1:n),mask(n/2+1:n))
      end if
  end function  
  
  pure function fsum_chunk_1d_sp(a) result(sout)
      integer, parameter :: wp = sp
      integer, parameter :: chunk = chunk32
      real(wp), intent(in) :: a(:)
      real(wp) :: sout
      ! --
      real(wp) :: abatch(chunk)
      integer :: i, dr, rr
      ! -----------------------------
      dr = size(a)/chunk
      rr = size(a) - dr*chunk
      
      abatch(:) = 0.0_wp
      do concurrent( i = 1:dr )
        abatch(1:chunk) = abatch(1:chunk) + a(chunk*i-chunk+1:chunk*i)
      end do
      abatch(1:rr) = abatch(1:rr) + a(size(a)-rr+1:size(a))
      
      sout = 0.0_wp
      do concurrent( i = chunk/2:1:-1 )
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
      integer :: i, dr, rr
      ! -----------------------------
      dr = size(a)/chunk
      rr = size(a) - dr*chunk
      
      abatch(:) = 0.0_wp
      do concurrent( i = 1:dr )
        abatch(1:chunk) = abatch(1:chunk) + a(chunk*i-chunk+1:chunk*i)
      end do
      abatch(1:rr) = abatch(1:rr) + a(size(a)-rr+1:size(a))

      sout = 0.0_wp
      do concurrent( i = chunk/2:1:-1 )
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
      integer :: i, dr, rr
      ! -----------------------------
      dr = size(a)/chunk
      rr = size(a) - dr*chunk
      
      abatch(:) = 0.0_wp
      do concurrent( i = 1:dr )
        where(mask(chunk*i-chunk+1:chunk*i)) abatch(1:chunk) = abatch(1:chunk) + a(chunk*i-chunk+1:chunk*i)
      end do
      where(mask(size(a)-rr+1:size(a))) abatch(1:rr) = abatch(1:rr) + a(size(a)-rr+1:size(a))

      sout = 0.0_wp
      do concurrent( i = chunk/2:1:-1)
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
      integer :: i, dr, rr
      ! -----------------------------
      dr = size(a)/chunk
      rr = size(a) - dr*chunk
      
      abatch(:) = 0.0_wp
      do concurrent( i = 1:dr )
        where(mask(chunk*i-chunk+1:chunk*i)) abatch(1:chunk) = abatch(1:chunk) + a(chunk*i-chunk+1:chunk*i)
      end do
      where(mask(size(a)-rr+1:size(a))) abatch(1:rr) = abatch(1:rr) + a(size(a)-rr+1:size(a))

      sout = 0.0_wp
      do concurrent( i = chunk/2:1:-1)
        sout = sout + abatch(i)+abatch(chunk/2+i)
      end do
  end function
end module fast_sum