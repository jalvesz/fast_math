  real(wp), intent(in) :: a
  real(wp), intent(inout) :: s
  real(wp), intent(inout) :: c
  logical, intent(in) :: m
  ! -- internal variables
  real(wp) :: t, y 
  y = a - c 
  t = s + y 
  c = (t - s) - y 
  s = merge( s , t , m )