!
! SPDX-FileCopyrightText: 2016-2022 Federico Perini <perini@wisc.edu>
!
! SPDX-License-Identifier: MIT
!
!   ***********************************************************************************************
!> @brief A module to compute FAST logarithm functions, based on Perini and Reitz, "Fast         **
!>        approximations of exponential and logarithm functions combined with efficient          **
!>        storage/retrieval for combustion kinetics calculations" Comb Flame 194(2018), 37-51.   **
!   ***********************************************************************************************
module fast_log
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    private
    
    public :: flog_p3, flog_p5

    interface flog_p3
        module procedure flog_p3_dp
    end interface
    interface flog_p5
        module procedure flog_p5_dp
    end interface

#ifdef __NVCOMPILER
    include 'utilities/nvidia_shift_interface.inc'
#endif
    
contains

    elemental function flog_p3_dp(x) result(y)
        integer, parameter :: wp = dp
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp) :: xi,xf
        integer(wp) :: iwp
        integer(wp), parameter :: mantissa_left  = 2_wp**52
        integer(wp), parameter :: mantissa       = -9218868437227405313_wp ! not(shiftl(2047_wp,52))
        integer(wp), parameter :: bias           = 1023_wp
        integer(wp), parameter :: ishift         = mantissa_left*bias

        real(wp), parameter :: log2         = log(2._wp)
        real(wp), parameter :: rlog2        = 1._wp/log2
        real(wp), parameter :: sqrt2        = sqrt(2._wp)
        real(wp), parameter :: s(3)= [rlog2,3.0_wp-2.5_wp*rlog2,1.5_wp*rlog2-2.0_wp]
        !-------------------------------------------------
        iwp = transfer(x,iwp)
        xi = shiftr(iwp,52)-bias

        ! Take mantissa part only
        xf = transfer(iand(iwp,mantissa)+ishift,xf)-1._wp

        ! Apply cubic polynomial
        xf = xf*(s(1)+xf*(s(2)+xf*s(3)))

        ! Compute log and Change of basis: log_2(x) -> log_e(x) = log2*log_2(x)
        y = (xf+xi)*log2

    end function flog_p3_dp

    elemental function flog_p5_dp(x) result(y)
        integer, parameter :: wp = dp
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp) :: xi,xf
        integer(wp) :: iwp
        integer(wp), parameter :: mantissa_left  = 2_wp**52
        integer(wp), parameter :: mantissa       = -9218868437227405313_wp ! not(shiftl(2047_wp,52))
        integer(wp), parameter :: bias           = 1023_wp
        integer(wp), parameter :: ishift         = mantissa_left*bias

        real(wp), parameter :: log2         = log(2._wp)
        real(wp), parameter :: rlog2        = 1._wp/log2
        real(wp), parameter :: sqrt2        = sqrt(2._wp)
        real(wp), parameter :: s(5)= [ 1.44269504088896e+0_wp,&
                                      -7.21347520444482e-1_wp,&
                                       4.42145354110618e-1_wp,&
                                      -2.12375830888126e-1_wp,&
                                       4.88829563330264e-2_wp]
        !-------------------------------------------------
        iwp = transfer(x,iwp)
        xi = shiftr(iwp,52)-bias

        ! Take mantissa part only
        xf = transfer(iand(iwp,mantissa)+ishift,xf)-1._wp

        ! Apply quintic polynomial
        xf = xf*(s(1)+xf*(s(2)+xf*(s(3)+xf*(s(4)+xf*s(5)))))

        ! Compute log and Change of basis: log_2(x) -> log_e(x) = log2*log_2(x)
        y = (xf+xi)*log2

    end function flog_p5_dp

#ifdef __NVCOMPILER
    include 'utilities/nvidia_shift.inc'
#endif
end module fast_log