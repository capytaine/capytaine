! Copyright (C) 2017-2024 Matthieu Ancellin
! See LICENSE file at <https://github.com/capytaine/libDelhommeau>

MODULE CONSTANTS

  USE ieee_arithmetic
  USE FLOATING_POINT_PRECISION, ONLY: PRE

  IMPLICIT NONE

  REAL(KIND=PRE), PARAMETER    :: ZERO = 0
  COMPLEX(KIND=PRE), PARAMETER :: CZERO = CMPLX(ZERO, ZERO, KIND=PRE)
  REAL(KIND=PRE), PARAMETER    :: ONE = 1

  REAL(KIND=PRE), PARAMETER    :: PI = 3.141592653588979 ! π
  REAL(KIND=PRE), PARAMETER    :: EULER_GAMMA = 0.5772156649
  REAL(KIND=PRE), PARAMETER    :: LOG_2 = LOG(REAL(2d0, kind=pre))
  COMPLEX(KIND=PRE), PARAMETER :: II = (ZERO, ONE) ! Imaginary unit

  integer, parameter :: nb_tabulated_values = 5

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Parameters for different variants of the Green function.
  ! The actual number should be irrelevant as long as they are different.

  ! Values for "tabulation_grid_shape"
  integer, parameter :: LEGACY_GRID = 0  ! Nemoh 2
  integer, parameter :: SCALED_NEMOH3_GRID = 1

  ! Values for "gf_singularities"
  integer, parameter :: HIGH_FREQ = 0  ! legacy from Nemoh
  integer, parameter :: LOW_FREQ = 1  ! aka XieDelhommeau
  integer, parameter :: LOW_FREQ_WITH_RANKINE_PART = 2  ! like LOW_FREQ but with a term integrated exactly

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

  PURE LOGICAL FUNCTION is_infinity(x)
    REAL(KIND=PRE), INTENT(IN) :: x
    is_infinity = (.NOT. ieee_is_finite(x))
  END FUNCTION

END MODULE CONSTANTS
