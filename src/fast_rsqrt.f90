!
! SPDX-FileCopyrightText: 2016-2022 Federico Perini <perini@wisc.edu>
!
! SPDX-License-Identifier: MIT
!
!   ***********************************************************************************************
!> @brief A FAST reciprocal of a square root, 1/sqrt(x), based on Perini and Reitz, "Fast        **
!>        approximations of exponential and logarithm functions combined with efficient          **
!>        storage/retrieval for combustion kinetics calculations" Comb Flame 194(2018), 37-51.   **
!   ***********************************************************************************************
module fast_rsqrt
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    private
    
    public :: frsqrt

    interface frsqrt
        !! Retranscript of the original Quake III Arena, see https://en.wikipedia.org/wiki/Fast_inverse_square_root
        !! for pure reference
        module procedure frsqrt_dp
        module procedure frsqrt_sp
    end interface
    
contains

    elemental function frsqrt_sp(x) result(y)
        integer, parameter :: wp = sp
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp)    :: x2,y2
        integer(wp) :: i
        integer(wp), parameter :: magic = int(Z'5f3759df',kind=wp)
        !-------------------------------------------------
        x2 = 0.5_wp*x
        i  = transfer(x,i)
        i  = magic - shiftr(i,1)
        y2 = transfer(i,y)

        ! Perform one Newton-Raphson step
        y  = y2*(1.5_wp-x2*y2*y2)

    end function frsqrt_sp

    elemental function frsqrt_dp(x) result(y)
        integer, parameter :: wp = dp
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp)    :: x2,y2
        integer(wp) :: i
        integer(wp), parameter :: magic = 6910469410427058089_wp
        !-------------------------------------------------
        x2 = 0.5_wp*x
        i  = transfer(x,i)
        i  = magic - shiftr(i,1)
        y2 = transfer(i,y)

        ! Perform one Newton-Raphson step
        y  = y2*(1.5_wp-x2*y2*y2)

    end function frsqrt_dp
    
end module fast_rsqrt