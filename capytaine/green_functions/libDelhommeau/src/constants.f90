! Copyright (C) 2017-2019 Matthieu Ancellin
! See LICENSE file at <https://github.com/mancellin/libDelhommeau>

MODULE CONSTANTS

  USE FLOATING_POINT_PRECISION, ONLY: PRE

  IMPLICIT NONE

  REAL(KIND=PRE), PARAMETER :: ZERO = 0
  REAL(KIND=PRE), PARAMETER :: ONE = 1

  REAL(KIND=PRE), PARAMETER :: PI = 3.141592653588979 ! π
  COMPLEX(KIND=PRE), PARAMETER :: II = (0, 1)         ! Imaginary unit

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Parameters for different variants of the Green function.
  ! The actual number should be irrelevant as long as they are different.

  ! Values for "tabulation_grid_shape"
  integer, parameter :: LEGACY_GRID = 0  ! Nemoh 2
  integer, parameter :: SCALED_NEMOH3_GRID = 1

  ! Values for "gf_singularities"
  integer, parameter :: HIGH_FREQ = 0  ! legacy from Nemoh
  integer, parameter :: LOW_FREQ = 1  ! aka XieDelhommeau


END MODULE CONSTANTS
