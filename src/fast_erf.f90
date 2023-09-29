!
! SPDX-FileCopyrightText: 2023 Transvalor S.A.
!
! SPDX-License-Identifier: MIT
!
module fast_erf
    !! Source: https://fortran-lang.discourse.group/t/fastgpt-faster-than-pytorch-in-300-lines-of-fortran/5385/31
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    private
    
    public :: ferf
    
    interface ferf
        module procedure ferf_sp
        module procedure ferf_dp
    end interface
    
    contains
    
    elemental function ferf_sp( x ) result( y )
        integer, parameter :: wp = sp
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp) :: abs_x
        !-------------------------------------------------
        abs_x = abs(x)
        y = 1 - 1 / (1+ 0.278393_wp*abs_x + 0.230389_wp*abs_x**2 + 0.000972_wp*abs_x**3 + 0.078108_wp*abs_x**4)**4
        y = y * sign(1.0_wp,x)
    end function

    elemental function ferf_dp( x ) result( y )
        integer, parameter :: wp = dp
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp) :: abs_x
        !-------------------------------------------------
        abs_x = abs(x)
        y = 1 - 1 / (1+ 0.278393_wp*abs_x + 0.230389_wp*abs_x**2 + 0.000972_wp*abs_x**3 + 0.078108_wp*abs_x**4)**4
        y = y * sign(1.0_wp,x)
    end function
    
end module