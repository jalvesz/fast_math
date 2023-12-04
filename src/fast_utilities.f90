module fast_utilities
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    implicit none
    private

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