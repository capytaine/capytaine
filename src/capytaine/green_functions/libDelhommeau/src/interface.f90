module interface
  ! This module is meant to keep functions that are only meant to be used for Python interfacing.

  use floating_point_precision, only: pre
  use constants
  use mesh_types, only: face_type, create_face, deallocate_face
  use green_wave, only: &
    wave_part_infinite_depth, &
    integral_of_wave_part_infinite_depth_typed => integral_of_wave_part_infinite_depth, &
    integral_of_wave_parts_finite_depth_typed => integral_of_wave_parts_finite_depth, &
    integral_of_prony_decomp_finite_depth_typed => integral_of_prony_decomp_finite_depth, &
    compute_dispersion_roots_typed => compute_dispersion_roots
  use green_Rankine, only: integral_of_Rankine_typed => integral_of_Rankine

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

  subroutine integral_of_Rankine(                &
      x,                                         &
      vertices, center, normal,                  &
      area, radius,                              &
      quadrature_points, quadrature_weights,                 &
      derivative_with_respect_to_first_variable, &
      int_G, int_nablaG)

    real(kind=pre), dimension(3),             intent(in) :: x
    real(kind=pre), dimension(:, :),          intent(in) :: vertices
    real(kind=pre), dimension(3),             intent(in) :: center, normal
    real(kind=pre),                           intent(in) :: area, radius
    real(kind=pre), dimension(:, :),          intent(in) :: quadrature_points
    real(kind=pre), dimension(:),             intent(in) :: quadrature_weights
    logical,                                  intent(in) :: derivative_with_respect_to_first_variable

    real(kind=pre),                           intent(out) :: int_G
    real(kind=pre), dimension(3),             intent(out) :: int_nablaG

    ! Local variables
    type(face_type) :: face

    ! Convert array inputs to Face type
    face = create_face(vertices, center, normal, area, radius, quadrature_points, quadrature_weights)

    ! Call the face_type version
    call integral_of_Rankine_typed(              &
      x,                                         &
      face,                                      &
      derivative_with_respect_to_first_variable, &
      int_G,                                     &
      int_nablaG                                 &
    )

    ! Clean up
    call deallocate_face(face)

  end subroutine integral_of_Rankine

  ! =====================================================================

  subroutine integral_of_wave_part_infinite_depth(               &
      x,                                                         &
      vertices, center, normal,                                  &
      area, radius,                                              &
      quadrature_points, quadrature_weights,                     &
      wavenumber,                                                &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      derivative_with_respect_to_first_variable,                 &
      int_G, int_nablaG)

    real(kind=pre), dimension(3),          intent(in) :: x
    real(kind=pre), dimension(:, :),       intent(in) :: vertices
    real(kind=pre), dimension(3),          intent(in) :: center, normal
    real(kind=pre),                        intent(in) :: area, radius
    real(kind=pre), dimension(:, :),       intent(in) :: quadrature_points
    real(kind=pre), dimension(:),          intent(in) :: quadrature_weights

    real(kind=pre),                        intent(in) :: wavenumber
    integer,                               intent(in) :: gf_singularities
    integer,                               intent(in) :: tabulation_nb_integration_points
    integer,                               intent(in) :: tabulation_grid_shape
    real(kind=pre), dimension(:),          intent(in) :: tabulated_r_range
    real(kind=pre), dimension(:),          intent(in) :: tabulated_z_range
    real(kind=pre), dimension(:, :, :),    intent(in) :: tabulated_integrals
    logical,                               intent(in) :: derivative_with_respect_to_first_variable

    complex(kind=pre),                     intent(out) :: int_G
    complex(kind=pre), dimension(3),       intent(out) :: int_nablaG

    ! Local variables
    type(face_type) :: face

    ! Convert array inputs to Face type
    face = create_face(vertices, center, normal, area, radius, quadrature_points, quadrature_weights)

    ! Call the face_type version
    call integral_of_wave_part_infinite_depth_typed(             &
      x,                                                         &
      face,                                                      &
      wavenumber,                                                &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      derivative_with_respect_to_first_variable,                 &
      int_G, int_nablaG                                          &
    )

    ! Clean up
    call deallocate_face(face)

  end subroutine integral_of_wave_part_infinite_depth

  ! =====================================================================

  subroutine integral_of_wave_parts_finite_depth                 &
      (x,                                                        &
      vertices, center, normal,                                  &
      area, radius,                                              &
      quadrature_points, quadrature_weights,                     &
      wavenumber, depth,                                         &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      derivative_with_respect_to_first_variable,                 &
      int_G, int_nablaG                                          &
      )

    ! Inputs
    real(kind=pre), dimension(3),          intent(in) :: x
    real(kind=pre), dimension(:, :),       intent(in) :: vertices
    real(kind=pre), dimension(3),          intent(in) :: center, normal
    real(kind=pre),                        intent(in) :: area, radius
    real(kind=pre), dimension(:, :),       intent(in) :: quadrature_points
    real(kind=pre), dimension(:),          intent(in) :: quadrature_weights
    real(kind=pre),                        intent(in) :: wavenumber, depth
    logical,                               intent(in) :: derivative_with_respect_to_first_variable

    ! Tabulated data
    integer,                               intent(in) :: tabulation_nb_integration_points
    integer,                               intent(in) :: tabulation_grid_shape
    real(kind=pre), dimension(:),          intent(in) :: tabulated_r_range
    real(kind=pre), dimension(:),          intent(in) :: tabulated_z_range
    real(kind=pre), dimension(:, :, :),    intent(in) :: tabulated_integrals

    integer,                               intent(in) :: gf_singularities

    ! Outputs
    complex(kind=pre),               intent(out) :: int_G  ! integral of the Green function over the panel.
    complex(kind=pre), dimension(3), intent(out) :: int_nablaG ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    type(face_type) :: face

    ! Convert array inputs to Face type
    face = create_face(vertices, center, normal, area, radius, quadrature_points, quadrature_weights)

    ! Call the face_type version
    call integral_of_wave_parts_finite_depth_typed(              &
      x,                                                         &
      face,                                                      &
      wavenumber, depth,                                         &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      derivative_with_respect_to_first_variable,                 &
      int_G, int_nablaG                                          &
    )

    ! Clean up
    call deallocate_face(face)

  end subroutine integral_of_wave_parts_finite_depth

  ! =====================================================================

  subroutine integral_of_prony_decomp_finite_depth &
      (x,                                          &
      vertices, center, normal,                    &
      area, radius,                                &
      quadrature_points, quadrature_weights,       &
      depth,                                       &
      prony_decomposition,                         &
      derivative_with_respect_to_first_variable,   &
      int_G, int_nablaG                            &
      )

    ! Inputs
    real(kind=pre), dimension(3),          intent(in) :: x
    real(kind=pre), dimension(:, :),       intent(in) :: vertices
    real(kind=pre), dimension(3),          intent(in) :: center, normal
    real(kind=pre),                        intent(in) :: area, radius
    real(kind=pre), dimension(:, :),       intent(in) :: quadrature_points
    real(kind=pre), dimension(:),          intent(in) :: quadrature_weights
    real(kind=pre),                        intent(in) :: depth
    logical,                               intent(in) :: derivative_with_respect_to_first_variable

    real(kind=pre), dimension(:, :),       intent(in) :: prony_decomposition

    ! Outputs
    complex(kind=pre),               intent(out) :: int_G
    complex(kind=pre), dimension(3), intent(out) :: int_nablaG

    ! Local variables
    type(face_type) :: face

    ! Convert array inputs to Face type
    face = create_face(vertices, center, normal, area, radius, quadrature_points, quadrature_weights)

    ! Call the face_type version
    call integral_of_prony_decomp_finite_depth_typed(            &
      x,                                                         &
      face,                                                      &
      depth,                                                     &
      prony_decomposition,                                       &
      derivative_with_respect_to_first_variable,                 &
      int_G, int_nablaG                                          &
    )

    ! Clean up
    call deallocate_face(face)

  end subroutine integral_of_prony_decomp_finite_depth

  ! =====================================================================

  function compute_dispersion_roots(nk, omega2_over_g, depth) &
      result(roots_of_dispersion_relationship)
    ! Just because the function is in a file that is not interfaceable to Python.
    integer, intent(in) :: nk
    real(kind=pre), intent(in) :: omega2_over_g, depth
    real(kind=pre), dimension(nk) :: roots_of_dispersion_relationship

    roots_of_dispersion_relationship = compute_dispersion_roots_typed(nk, omega2_over_g, depth)

  end function compute_dispersion_roots

  ! =====================================================================

end module interface
