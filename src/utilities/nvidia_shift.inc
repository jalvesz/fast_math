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