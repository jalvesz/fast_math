!
! SPDX-FileCopyrightText: 2023 Transvalor S.A.
!
! SPDX-License-Identifier: MIT
!
module fast_tanh
    !! Source: https://fortran-lang.discourse.group/t/fastgpt-faster-than-pytorch-in-300-lines-of-fortran/5385/31
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    private
    
    public :: ftanh
    
    interface ftanh
        module procedure ftanh_sp
        module procedure ftanh_dp
    end interface
    
    contains
    
    elemental function ftanh_sp( x ) result( y )
        integer, parameter :: wp = sp
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp) :: x2, a, b
        !---------------------------------------------
        if (x > 5_wp) then
            y = 1_wp
        elseif (x < -5_wp) then
            y = -1_wp
        else
            x2 = x*x
            a = x * (135135.0_wp + x2 * (17325.0_wp + x2 * (378.0_wp + x2)))
            b = 135135.0_wp + x2 * (62370.0_wp + x2 * (3150.0_wp + x2 * 28.0_wp))
            y = a / b
        end if
    end function

    elemental function ftanh_dp( x ) result( y )
        integer, parameter :: wp = dp
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp) :: x2, a, b
        !---------------------------------------------
        if (x > 5_wp) then
            y = 1_wp
        elseif (x < -5_wp) then
            y = -1_wp
        else
            x2 = x*x
            a = x * (135135.0_wp + x2 * (17325.0_wp + x2 * (378.0_wp + x2)))
            b = 135135.0_wp + x2 * (62370.0_wp + x2 * (3150.0_wp + x2 * 28.0_wp))
            y = a / b
        end if
    end function
    
end module