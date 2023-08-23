!
! SPDX-FileCopyrightText: 2023 Transvalor S.A.
!
! SPDX-License-Identifier: MIT
!
module fast_math
  !! User API: All modules can be referenced from this module as entry point
  !-------------------------
  ! Basics
  use fast_sum
  use fast_dotp
  !-------------------------
  ! Trigonometric
  use fast_trigo
  !-------------------------
  ! Hyperbolic
  use fast_tanh
  use fast_erf

end module fast_math
