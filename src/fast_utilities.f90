module fast_utilities
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    private
    
    interface vkahans
      module procedure vkahans_sp
      module procedure vkahans_dp
    end interface
    public :: vkahans
    
#ifdef __NVCOMPILER
    interface shiftl
      module procedure shiftl_sp
      module procedure shiftl_dp
    end interface
    interface shiftr
      module procedure shiftr_sp
      module procedure shiftr_dp
    end interface
    public :: shiftl, shiftr
#endif
    
contains
    
elemental subroutine vkahans_sp(a,s,c)
  integer, parameter :: wp = sp
  real(wp), intent(in) :: a
  real(wp), intent(inout) :: s
  real(wp), intent(inout) :: c
  ! -- internal variables
  real(wp) :: t, y    
  y = a - c
  t = s + y
  c = (t - s) - y
  s = t
end subroutine  

elemental subroutine vkahans_dp(a,s,c)
  integer, parameter :: wp = dp
  real(wp), intent(in) :: a
  real(wp), intent(inout) :: s
  real(wp), intent(inout) :: c
  ! -- internal variables
  real(wp) :: t, y    
  y = a - c
  t = s + y
  c = (t - s) - y
  s = t
end subroutine

#ifdef __NVCOMPILER
elemental integer(sp) function shiftr_sp( I , shift )
  integer(sp), intent(in) :: I 
  integer, intent(in) :: shift
  shiftr_sp = rshift( I, shift )
end function

elemental integer(dp) function shiftr_dp( I , shift )
  integer(dp), intent(in) :: I 
  integer, intent(in) :: shift
  shiftr_dp = rshift( I, shift )
end function

elemental integer(sp) function shiftl_sp( I , shift )
  integer(sp), intent(in) :: I 
  integer, intent(in) :: shift
  shiftl_sp = lshift( I, shift )
end function

elemental integer(dp) function shiftl_dp( I , shift )
  integer(dp), intent(in) :: I 
  integer, intent(in) :: shift
  shiftl_dp = lshift( I, shift )
end function
#endif
end module fast_utilities