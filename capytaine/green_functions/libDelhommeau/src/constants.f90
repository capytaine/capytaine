! Copyright (C) 2017-2019 Matthieu Ancellin
! See LICENSE file at <https://github.com/mancellin/libDelhommeau>

MODULE CONSTANTS

  USE ieee_arithmetic
  USE FLOATING_POINT_PRECISION, ONLY: PRE

  IMPLICIT NONE

  REAL(KIND=PRE), PARAMETER :: ZERO = 0
  REAL(KIND=PRE), PARAMETER :: ONE = 1

  REAL(KIND=PRE), PARAMETER :: PI = 3.141592653588979 ! π
  REAL(KIND=PRE), PARAMETER :: EULER_GAMMA = 0.5772156649
  REAL(KIND=PRE), PARAMETER :: LOG_2 = LOG(REAL(2d0, kind=pre))
  COMPLEX(KIND=PRE), PARAMETER :: II = (0, 1)         ! Imaginary unit

contains

  PURE LOGICAL FUNCTION is_infinity(x)
    REAL(KIND=PRE), INTENT(IN) :: x
    is_infinity = (.NOT. ieee_is_finite(x))
  END FUNCTION

END MODULE CONSTANTS
