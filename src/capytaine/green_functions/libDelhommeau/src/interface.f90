module interface
  ! This module is meant to keep functions that are only meant to be used for Python interfacing.

  use floating_point_precision, only: pre
  use constants
  use mesh_types
  use green_wave, only: wave_part_infinite_depth
  use green_Rankine, only: integral_of_Rankine

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


  ! =====================================================================

  subroutine integral_of_Rankine_array( &
      x,                                &
      vertices, center, normal,         &
      area, radius,                     &
      quad_points, quad_weights,        &
      derivative_with_respect_to_first_variable, &
      int_G, int_nablaG)

    real(kind=pre), dimension(3),             intent(in) :: x
    real(kind=pre), dimension(:, :),          intent(in) :: vertices
    real(kind=pre), dimension(3),             intent(in) :: center, normal
    real(kind=pre),                           intent(in) :: area, radius
    real(kind=pre), dimension(:, :),          intent(in) :: quad_points
    real(kind=pre), dimension(:),             intent(in) :: quad_weights
    logical,                                  intent(in) :: derivative_with_respect_to_first_variable

    real(kind=pre),                           intent(out) :: int_G
    real(kind=pre), dimension(3),             intent(out) :: int_nablaG

    ! Local variables
    type(Face) :: face

    ! Convert array inputs to Face type
    face = create_face(vertices, center, normal, area, radius, quad_points, quad_weights)

    ! Call the Face-based version
    call integral_of_Rankine(x, face, derivative_with_respect_to_first_variable, int_G, int_nablaG)

    ! Clean up
    call destroy_face(face)

  end subroutine integral_of_Rankine_array

end module
