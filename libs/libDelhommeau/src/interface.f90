! Copyright 2026 Capytaine developers
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
module interface
  ! This module is meant to keep functions that are only meant to be used for Python interfacing.

  use floating_point_precision, only: pre
  use constants
  use green_wave, only: wave_part_infinite_depth

  implicit none

  public

contains

  ! =====================================================================

  SUBROUTINE VECTORIZED_WAVE_PART_INFINITE_DEPTH                 &
      (nb_points_pairs, x0i, x0j, wavenumber,                    &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      g, nablag)
    ! Vectorized version of the WAVE_PART_INFINITE_DEPTH, taking as input two lists of points and computing the Green
    ! function in a Fortran loop.
    ! Meant for testing and debugging, not currently used otherwise.

    ! Inputs
    integer,                                          intent(in) :: nb_points_pairs
    real(kind=pre), dimension(nb_points_pairs, 3),    intent(in) :: x0i, x0j
    real(kind=pre),                                   intent(in) :: wavenumber
    integer,                                          intent(in) :: gf_singularities

    ! tabulated data
    integer,                                          intent(in) :: tabulation_nb_integration_points
    integer,                                          intent(in) :: tabulation_grid_shape
    real(kind=pre), dimension(:),                     intent(in) :: tabulated_r_range
    real(kind=pre), dimension(:),                     intent(in) :: tabulated_z_range
    real(kind=pre), dimension(:, :, :),               intent(in) :: tabulated_integrals

    ! outputs
    complex(kind=pre), dimension(nb_points_pairs),    intent(out) :: g
    complex(kind=pre), dimension(nb_points_pairs, 3), intent(out) :: nablag

    ! local variables
    integer :: i

    do i = 1, nb_points_pairs
      call wave_part_infinite_depth(                               &
        x0i(i, :), x0j(i, :), wavenumber,                          &
        tabulation_nb_integration_points, tabulation_grid_shape,   &
        tabulated_r_range, tabulated_z_range, tabulated_integrals, &
        gf_singularities,                                          &
        g(i), nablag(i, :))
    enddo
  END SUBROUTINE VECTORIZED_WAVE_PART_INFINITE_DEPTH

  ! =====================================================================


end module
