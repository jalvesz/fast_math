module fast_erf
    !! Source: https://fortran-lang.discourse.group/t/fastgpt-faster-than-pytorch-in-300-lines-of-fortran/5385/31
    use iso_fortran_env
    implicit none
    private
    
    public :: ferf
    
    interface ferf
        module procedure ferf_r32
        module procedure ferf_r64
    end interface
    
    contains
    
    elemental function ferf_r32( x ) result( y )
        integer, parameter :: wp = real32
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp) :: abs_x
        !-------------------------------------------------
        abs_x = abs(x)
        y = 1 - 1 / (1+ 0.278393_wp*abs_x + 0.230389_wp*abs_x**2 + 0.000972_wp*abs_x**3 + 0.078108_wp*abs_x**4)**4
        y = y * sign(1.0_wp,x)
    end function

    elemental function ferf_r64( x ) result( y )
        integer, parameter :: wp = real64
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