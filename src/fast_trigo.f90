!
! SPDX-FileCopyrightText: 2023 Transvalor S.A.
!
! SPDX-License-Identifier: MIT
!
module fast_trigo
    !! Source for fast sine cosine: http://web.archive.org/web/20141220225551/http://forum.devmaster.net/t/fast-and-accurate-sine-cosine/9648
    use iso_fortran_env
    implicit none
    private
    
    public :: fsin, fcos, ftan
    public :: facos, facos_nvidia
    
    interface fcos
        module procedure fcos_r32
        module procedure fcos_r64
    end interface

    interface fsin
        module procedure fsin_r32
        module procedure fsin_r64
    end interface

    interface ftan
        module procedure ftan_r32
        module procedure ftan_r64
    end interface

    interface facos
        module procedure facos_r32
        module procedure facos_r64
    end interface
    
    interface facos_nvidia
    !! Source : https://developer.download.nvidia.com/cg/acos.html
        module procedure facos_nvidia_r32
        module procedure facos_nvidia_r64
    end interface
    
    contains
    
    elemental function fcos_r32( x ) result( y )
        integer, parameter :: wp = real32
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp), parameter :: half_pi =acos(-1.0_wp)/2
        !---------------------------------------------
        y = fsin_r32( half_pi - x )
    end function

    elemental function fcos_r64( x ) result( y )
        integer, parameter :: wp = real64
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp), parameter :: half_pi =acos(-1.0_wp)/2
        !---------------------------------------------
        y = fsin_r64( half_pi - x )
    end function

    elemental function fsin_r32( x ) result( y )
        integer, parameter :: wp = real32
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp), parameter :: twopi = 2*acos(-1.0_wp)
        real(wp), parameter :: invtwopi = 1.0_wp/twopi
        real(wp), parameter :: c1=4.0_wp/acos(-1.0_wp)
        real(wp), parameter :: c2=-4.0_wp/acos(-1.0_wp)**2
        real(wp), parameter :: c3=0.225_wp
        real(wp) :: x0
        !---------------------------------------------
        x0 = x - (int(x*invtwopi,kind=1) * twopi) 
        y = c1*x0+c2*x0*abs(x0)
        y = c3*(y*abs(y)-y)+y
    end function

    elemental function fsin_r64( x ) result( y )
        integer, parameter :: wp = real64
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp), parameter :: twopi = 2*acos(-1.0_wp)
        real(wp), parameter :: invtwopi = 1.0_wp/twopi
        real(wp), parameter :: c1=4.0_wp/acos(-1.0_wp)
        real(wp), parameter :: c2=-4.0_wp/acos(-1.0_wp)**2
        real(wp), parameter :: c3=0.225_wp
        real(wp) :: x0
        !---------------------------------------------
        x0 = x - (int(x*invtwopi,kind=1) * twopi) 
        y = c1*x0+c2*x0*abs(x0)
        y = c3*(y*abs(y)-y)+y
    end function

    elemental function ftan_r32( x ) result( y )
        integer, parameter :: wp = real32
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp), parameter :: pi = acos(-1.0_wp)
        real(wp), parameter :: invpi = 1.0_wp/acos(-1.0_wp)
        real(wp) :: x0, xsq
        !-------------------------------------------------
        x0 = x - (int(x*invpi,kind=1) * pi) 
        xsq = x0 * x0
        y = x0 * (2.471688400562703_wp - 0.189759681063053_wp * xsq) / &
                 (2.4674011002723397_wp - xsq)
    end function

    elemental function ftan_r64( x ) result( y )
        integer, parameter :: wp = real64
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        real(wp), parameter :: pi = acos(-1.0_wp)
        real(wp), parameter :: invpi = 1.0_wp/acos(-1.0_wp)
        real(wp) :: x0, xsq
        !-------------------------------------------------
        x0 = x - (int(x*invpi,kind=1) * pi) 
        xsq = x0 * x0
        y = x * (2.471688400562703_wp - 0.189759681063053_wp * xsq) / &
                (2.4674011002723397_wp - xsq)
    end function

    !====================================================
    ! Inverse
    !====================================================

    elemental function facos_r32( x ) result( y )
      integer, parameter :: wp = real32
      real(wp), intent(in) :: x
      real(wp) :: y
      !---------------------------------------------
      y =  (-0.69813170079773212_wp * x * x - 0.87266462599716477_wp) * x + 1.5707963267948966_wp;
    end function

    elemental function facos_r64( x ) result( y )
      integer, parameter :: wp = real64
      real(wp), intent(in) :: x
      real(wp) :: y
      !---------------------------------------------
      y =  (-0.69813170079773212_wp * x * x - 0.87266462599716477_wp) * x + 1.5707963267948966_wp;
    end function

    elemental function  facos_nvidia_r32( x ) result( y )
        integer, parameter :: wp = real32
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        integer(1) :: negate
        real(wp) :: xp
        !---------------------------------------------
        negate = merge( 1_1 , 0_1 , x < 0_wp )
        xp = abs(x)
        y = -0.0187293_wp * xp + 0.0742610_wp
        y = y * xp - 0.2121144_wp
        y = y * xp + 1.5707288_wp
        y = y * sqrt(1_wp-xp)
        y = y + negate * (- 2.0_wp * y + 3.14159265358979_wp)
    end function

    elemental function  facos_nvidia_r64( x ) result( y )
        integer, parameter :: wp = real64
        real(wp), intent(in) :: x
        real(wp) :: y
        !-- Internal Variables
        integer(1) :: negate
        real(wp) :: xp
        !---------------------------------------------
        negate = merge( 1_1 , 0_1 , x < 0_wp )
        xp = abs(x)
        y = -0.0187293_wp * xp + 0.0742610_wp
        y = y * xp - 0.2121144_wp
        y = y * xp + 1.5707288_wp
        y = y * sqrt(1_wp-xp)
        y = y + negate * (- 2.0_wp * y + 3.14159265358979_wp)
    end function
    
end module