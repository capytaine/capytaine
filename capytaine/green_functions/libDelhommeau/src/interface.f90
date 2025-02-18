module interface
  ! This module is meant to keep functions that are only meant to be used for Python interfacing.

  use floating_point_precision, only: pre
  use constants
  use green_wave, only: wave_part_infinite_depth, wave_part_finite_depth

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

  SUBROUTINE VECTORIZED_WAVE_PART_FINITE_DEPTH                   &
      (nb_points_pairs, x0i, x0j, wavenumber, depth,             &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      nexp, ambda, ar,                                           &
      g, nablag_sym, nablag_antisym)
    ! Vectorized version of the WAVE_PART_FINITE_DEPTH, taking as input two lists of points and computing the Green
    ! function in a Fortran loop.
    ! Meant for testing and debugging, not currently used otherwise.

    ! Inputs
    integer,                                          intent(in) :: nb_points_pairs
    real(kind=pre), dimension(nb_points_pairs, 3),    intent(in) :: x0i, x0j
    real(kind=pre),                                   intent(in) :: wavenumber, depth

    ! tabulated data
    integer,                                          intent(in) :: tabulation_nb_integration_points
    integer,                                          intent(in) :: tabulation_grid_shape
    real(kind=pre), dimension(:),                     intent(in) :: tabulated_r_range
    real(kind=pre), dimension(:),                     intent(in) :: tabulated_z_range
    real(kind=pre), dimension(:, :, :),               intent(in) :: tabulated_integrals

    ! Prony decomposition for finite depth
    integer,                                          intent(in) :: nexp
    real(kind=pre), dimension(nexp),                  intent(in) :: ambda, ar

    ! outputs
    complex(kind=pre), dimension(nb_points_pairs),    intent(out) :: g
    complex(kind=pre), dimension(nb_points_pairs, 3), intent(out) :: nablag_sym, nablag_antisym

    ! local variables
    integer :: i

    do i = 1, nb_points_pairs
      call wave_part_finite_depth( &
        x0i(i, :), x0j(i, :), wavenumber, depth,                   &
        tabulation_nb_integration_points, tabulation_grid_shape,   &
        tabulated_r_range, tabulated_z_range, tabulated_integrals, &
        nexp, ambda, ar,                                           &
        g(i), nablag_sym(i, :), nablag_antisym(i, :))
    enddo
  END SUBROUTINE VECTORIZED_WAVE_PART_FINITE_DEPTH

  ! =====================================================================

end module
