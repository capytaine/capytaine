! Copyright (C) 2017-2024 Matthieu Ancellin
! See LICENSE file at <https://github.com/capytaine/libDelhommeau>

module matrices

  use floating_point_precision, only: pre
  use constants

  use green_Rankine
  use green_wave

  implicit none

contains

  ! =====================================================================

  subroutine build_matrices(                                          &
      nb_collocation_points, collocation_points, dot_product_normals, &
      nb_vertices, nb_faces, vertices, faces,                         &
      centers, normals, areas, radiuses,                              &
      nb_quad_points, quad_points, quad_weights,                      &
      wavenumber, depth,                                              &
      tabulation_nb_integration_points,                               &
      tabulation_grid_shape,                                          &
      tabulated_r_range, tabulated_z_range,                           &
      tabulated_integrals,                                            &
      finite_depth_method, prony_decomposition, dispersion_roots,     &
      gf_singularities, adjoint_double_layer,                         &
      S, K)

    ! Mesh data
    integer,                                              intent(in) :: nb_collocation_points
    real(kind=pre), dimension(nb_collocation_points, 3),  intent(in) :: collocation_points
    real(kind=pre), dimension(:, :),           intent(in) :: dot_product_normals
    ! If adjoint_double_layer:     size(dot_product_normals) == (nb_collocation_points, 3)
    ! If not adjoint_double_layer: size(dot_product_normals) == (nb_faces, 3)

    ! dot_product_normals might be identical to normals, especially when adjoint_double_layer is False.
    ! The former is used when computing the dot product with the normal vector in the double layer or adjoint double layer operator
    ! (D or K matrices). The latter is only used when computing the exact formula to integrate the Rankine part of the Green
    ! function over a face.
    ! Hence, both variables fulfill different role, and may or may not be identical.

    integer, intent(in)                                   :: nb_vertices, nb_faces
    real(kind=pre), dimension(nb_vertices, 3), intent(in) :: vertices
    integer,        dimension(nb_faces, 4),    intent(in) :: faces
    real(kind=pre), dimension(nb_faces, 3),    intent(in) :: centers, normals
    real(kind=pre), dimension(nb_faces),       intent(in) :: areas, radiuses
    integer,                                                intent(in) :: nb_quad_points
    real(kind=pre), dimension(nb_faces, nb_quad_points, 3), intent(in) :: quad_points
    real(kind=pre), dimension(nb_faces, nb_quad_points),    intent(in) :: quad_weights

    ! Solver parameters
    integer,                                  intent(in) :: gf_singularities
    logical,                                  intent(in) :: adjoint_double_layer

    real(kind=pre),                           intent(in) :: wavenumber, depth

    ! Tabulated values for the wave part of the Green function
    integer,                                  intent(in) :: tabulation_grid_shape
    integer,                                  intent(in) :: tabulation_nb_integration_points
    real(kind=pre), dimension(:),             intent(in) :: tabulated_r_range
    real(kind=pre), dimension(:),             intent(in) :: tabulated_z_range
    real(kind=pre), dimension(:, :, :),       intent(in) :: tabulated_integrals

    integer,                                  intent(in) :: finite_depth_method
    real(kind=pre), dimension(:, :),          intent(in) :: prony_decomposition  ! For Delhommeau's finite depth, dummy otherwise
    real(kind=pre), dimension(:),             intent(in) :: dispersion_roots  ! For FinGreen3D, dummy otherwise

    ! Outputs
    complex(kind=pre), dimension(:, :),    intent(inout) :: S
    complex(kind=pre), dimension(:, :, :), intent(inout) :: K

    ! Local variables
    integer                         :: I, J
    real(kind=pre)                  :: sign_reflected_Rankine
    real(kind=pre)                  :: int_G_rankine, diagonal_coef
    real(kind=pre), dimension(3)    :: int_nablaG_rankine
    complex(kind=pre)               :: int_G, int_G_wave
    complex(kind=pre), dimension(3) :: int_nablaG, int_nablaG_wave
    logical :: derivative_with_respect_to_first_variable, finite_depth, finite_wavenumber

    derivative_with_respect_to_first_variable = adjoint_double_layer
    ! When computing the adjoint double layer operator (K), the derivative of the Green function is computed with respect to its
    ! first variable (field point, often written x, or sometimes M in this code).
    ! When computing the double layer operator (D), the derivative of the Green function is computed with respect to its second
    ! variable (source point, often written xi, or sometimes M' in this code).

    finite_depth = (.not. is_infinity(depth))

    finite_wavenumber = ((ZERO < wavenumber) .and. (.not. is_infinity(wavenumber)))

    if (gf_singularities == HIGH_FREQ) then
      sign_reflected_Rankine = -ONE
    else
      sign_reflected_Rankine = +ONE
    endif

    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) &
    !$OMP&  PRIVATE(J, I, int_G, int_nablaG, int_G_Rankine, int_nablaG_Rankine, diagonal_coef, &
    !$OMP&          int_G_wave, int_nablaG_wave)
    do J = 1, nb_faces
      do I = 1, nb_collocation_points

        int_G = CZERO
        int_nablaG = CZERO

        !!!!!!!!!!!!!!!!!!
        !  Rankine part  !
        !!!!!!!!!!!!!!!!!!
        call integral_of_Rankine(                    &
          collocation_points(I, :),                  &
          vertices(faces(J, :), :),                  &
          centers(J, :),                             &
          normals(J, :),                             &
          areas(J),                                  &
          radiuses(J),                               &
          derivative_with_respect_to_first_variable, &
          int_G_Rankine, int_nablaG_Rankine          &
          )

        int_G = int_G + int_G_Rankine
        int_nablaG(:) = int_nablaG(:) + int_nablaG_Rankine(:)
        ! Not sign_reflected_Rankine!
        ! The sign of this term does not change in the high frequency limit.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Reflected Rankine part  !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (finite_depth .and. (finite_depth_method == LEGACY_FINITE_DEPTH)) then
          ! Reproduce behavior of legacy finite depth method
          call one_point_integral_of_reflected_Rankine( &
            collocation_points(I, :),                   &
            centers(J, :),                              &
            areas(J),                                   &
            derivative_with_respect_to_first_variable,  &
            [-ONE, ZERO],                               &
            int_G_Rankine,                              &
            int_nablaG_Rankine                          &
          )
        else
          call integral_of_reflected_Rankine(          &
            collocation_points(I, :),                  &
            vertices(faces(J, :), :),                  &
            centers(J, :),                             &
            normals(J, :),                             &
            areas(J),                                  &
            radiuses(J),                               &
            derivative_with_respect_to_first_variable, &
            [-ONE, ZERO],                              &
            int_G_Rankine,                             &
            int_nablaG_Rankine                         &
            )
        endif
        int_G = int_G + sign_reflected_Rankine * int_G_Rankine
        int_nablaG(:) = int_nablaG(:) + sign_reflected_Rankine * int_nablaG_Rankine(:)

        if (gf_singularities == LOW_FREQ_WITH_RANKINE_PART) then
          int_nablaG(3) = int_nablaG(3) + 2*wavenumber * int_G_Rankine
        endif

        if (.not. finite_depth) then

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !  Infinite depth wave part  !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (finite_wavenumber) then
            call integral_of_wave_part_infinite_depth                    &
              (collocation_points(I, :),                                 &
              centers(J, :), areas(J),                                   &
              quad_points(J, :, :), quad_weights(J, :),                  &
              wavenumber,                                                &
              tabulation_nb_integration_points, tabulation_grid_shape,   &
              tabulated_r_range, tabulated_z_range, tabulated_integrals, &
              gf_singularities,                                          &
              derivative_with_respect_to_first_variable,                 &
              int_G_wave, int_nablaG_wave                                &
              )
              int_G = int_G + int_G_wave
              int_nablaG(:) = int_nablaG(:) + int_nablaG_wave(:)
          ! else zero, assuming gf_singularities has been set up correctly, that is
          ! low_freq for wavenumber = 0 and high_freq for wavenumber = inf
          endif

        else  ! Finite depth

          !!!!!!!!!!!!!!!!!!!!!!!!!!
          !  FinGreen3D wave part  !
          !!!!!!!!!!!!!!!!!!!!!!!!!!
          if (finite_depth_method == FINGREEN3D_METHOD) then
            call integral_of_wave_part_fingreen3D                 &
              (collocation_points(I, :),                          &
              quad_points(J, :, :), quad_weights(J, :),           &
              wavenumber, depth, dispersion_roots,                &
              derivative_with_respect_to_first_variable,          &
              int_G_wave, int_nablaG_wave                         &
              )
              int_G = int_G + int_G_wave
              int_nablaG(:) = int_nablaG(:) + int_nablaG_wave(:)

          else  ! Delhommeau's finite depth method

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !  Supplementary Rankine parts in finite depth  !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 1. Reflection through sea bottom
            call integral_of_reflected_Rankine(          &
              collocation_points(I, :),                  &
              vertices(faces(J, :), :),                  &
              centers(J, :),                             &
              normals(J, :),                             &
              areas(J),                                  &
              radiuses(J),                               &
              derivative_with_respect_to_first_variable, &
              [-ONE, -2*depth],                          &
              int_G_Rankine, int_nablaG_Rankine          &
              )
            int_G = int_G + int_G_Rankine
            int_nablaG(:) = int_nablaG(:) + int_nablaG_Rankine(:)
            ! Not sign_reflected_Rankine!
            ! The sign of this term does not change in the high frequency limit.

            ! 2. Reflection through sea bottom and free surface
            call one_point_integral_of_reflected_Rankine( &
              collocation_points(I, :),                   &
              centers(J, :),                              &
              areas(J),                                   &
              derivative_with_respect_to_first_variable,  &
              [ONE, -2*depth],                            &
              int_G_Rankine, int_nablaG_Rankine           &
              )
            int_G = int_G + sign_reflected_Rankine * int_G_Rankine
            int_nablaG(:) = int_nablaG(:) + sign_reflected_Rankine * int_nablaG_Rankine(:)

            ! 3. Reflection through free surface and sea bottom
            call one_point_integral_of_reflected_Rankine( &
              collocation_points(I, :),                   &
              centers(J, :),                              &
              areas(J),                                   &
              derivative_with_respect_to_first_variable,  &
              [ONE, 2*depth],                             &
              int_G_Rankine, int_nablaG_Rankine           &
              )
            int_G = int_G + sign_reflected_Rankine * int_G_Rankine
            int_nablaG(:) = int_nablaG(:) + sign_reflected_Rankine * int_nablaG_Rankine(:)

            ! 4. Reflection through sea bottom and free surface and sea bottom again
            call one_point_integral_of_reflected_rankine( &
              collocation_points(I, :),                   &
              centers(J, :),                              &
              areas(J),                                   &
              derivative_with_respect_to_first_variable,  &
              [-ONE, -4*depth],                           &
              int_G_Rankine, int_nablaG_Rankine           &
              )
            int_G = int_G + sign_reflected_Rankine * int_G_Rankine
            int_nablaG(:) = int_nablaG(:) + sign_reflected_Rankine * int_nablaG_Rankine(:)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !  Delhommeau's finite depth wave term  !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (finite_wavenumber) then
              call integral_of_wave_parts_finite_depth                     &
                (collocation_points(I, :),                                 &
                centers(J, :), areas(J),                                   &
                quad_points(J, :, :), quad_weights(J, :),                  &
                wavenumber, depth,                                         &
                tabulation_nb_integration_points, tabulation_grid_shape,   &
                tabulated_r_range, tabulated_z_range, tabulated_integrals, &
                gf_singularities,                                          &
                derivative_with_respect_to_first_variable,                 &
                int_G_wave, int_nablaG_wave                                &
                )
              int_G = int_G + int_G_wave
              int_nablaG(:) = int_nablaG(:) + int_nablaG_wave(:)
            ! else zero, assuming gf_singularities has been set up correctly, that is
            ! low_freq for wavenumber = 0 and high_freq for wavenumber = inf
            endif

            call integral_of_prony_decomp_finite_depth        &
              (collocation_points(I, :),                      &
              vertices(faces(J, :), :),                       &
              centers(J, :),                                  &
              normals(J, :),                                  &
              areas(J),                                       &
              radiuses(J),                                    &
              depth,                                          &
              prony_decomposition,                            &
              derivative_with_respect_to_first_variable,      &
              int_G_wave, int_nablaG_wave                     &
              )
            int_G = int_G + int_G_wave
            int_nablaG(:) = int_nablaG(:) + int_nablaG_wave(:)
          endif  ! if Delhommeau's finite depth
        endif  ! if finite depth

        !!!!!!!!!!!!!!!!!!!
        !  Add to matrix  !
        !!!!!!!!!!!!!!!!!!!
        S(I, J) = MINUS_ONE_OVER_FOURPI * int_G

        if (size(K, 3) == 1) then  ! early_dot_product=True
          if (adjoint_double_layer) then
            K(I, J, 1) = MINUS_ONE_OVER_FOURPI * DOT_PRODUCT(dot_product_normals(I, :), int_nablaG(:))
          else
            K(I, J, 1) = MINUS_ONE_OVER_FOURPI * DOT_PRODUCT(dot_product_normals(J, :), int_nablaG(:))
          endif
        else
          K(I, J, :) = MINUS_ONE_OVER_FOURPI * int_nablaG(:)
        endif

      end do  ! loop on I
    end do  ! parallelized loop on J
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_diagonal_term(                               &
      nb_faces, centers, dot_product_normals, free_surface, K &
      )

    integer, intent(in)                                :: nb_faces
    real(kind=pre), dimension(nb_faces, 3), intent(in) :: centers
    real(kind=pre), dimension(nb_faces, 3), intent(in) :: dot_product_normals
    real(kind=pre),                         intent(in) :: free_surface
    complex(kind=pre), dimension(:, :, :),  intent(inout) :: K

    ! Local variables
    integer        :: i
    real(kind=pre) :: diagonal_coef

    !$OMP PARALLEL DO PRIVATE(i, diagonal_coef)
    do i = 1, nb_faces
      if (abs(centers(i, 3) - free_surface) < 1e-8) then  ! Panel on the free surface
        diagonal_coef = ONE
      else
        diagonal_coef = ONE/2
      endif

      if (size(K, 3) == 1) then  ! early_dot_product=True
        K(i, i, 1) = K(i, i, 1) + diagonal_coef
      else
        K(i, i, :) = K(i, i, :) + diagonal_coef * dot_product_normals(i, :)
        ! if (.not. adjoint_double_layer) then we should have used the jth normal instead of the i-th,
        ! such that later the dot product with dot_product_normals(j, :) gives 1.
        ! Except that on the diagonal, i==j, so there is no need to branch based on adjoint_double_layer.
      endif
    enddo

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_rankine_term_only(                                   &
      nb_collocation_points, collocation_points, dot_product_normals, &
      nb_vertices, nb_faces, vertices, faces,                         &
      centers, normals, areas, radiuses,                              &
      nb_quad_points, quad_points, quad_weights,                      &
      adjoint_double_layer,                                           &
      S, K)

    integer,                                                intent(in) :: nb_collocation_points
    real(kind=pre), dimension(nb_collocation_points, 3),    intent(in) :: collocation_points
    real(kind=pre), dimension(:, :),                        intent(in) :: dot_product_normals
    ! If adjoint_double_layer:     size(dot_product_normals) == (nb_collocation_points, 3)
    ! If not adjoint_double_layer: size(dot_product_normals) == (nb_faces, 3)

    integer,                                                intent(in) :: nb_faces, nb_vertices
    real(kind=pre), dimension(nb_vertices, 3),              intent(in) :: vertices
    integer,        dimension(nb_faces, 4),                 intent(in) :: faces
    real(kind=pre), dimension(nb_faces, 3),                 intent(in) :: centers, normals
    real(kind=pre), dimension(nb_faces),                    intent(in) :: areas, radiuses
    integer,                                                intent(in) :: nb_quad_points
    real(kind=pre), dimension(nb_faces, nb_quad_points, 3), intent(in) :: quad_points
    real(kind=pre), dimension(nb_faces, nb_quad_points),    intent(in) :: quad_weights

    logical,                                                intent(in) :: adjoint_double_layer

    complex(kind=pre), dimension(:, :),                     intent(inout) :: S
    complex(kind=pre), dimension(:, :, :),                  intent(inout) :: K

    ! Local variables
    real(kind=pre)               :: int_G_Rankine
    real(kind=pre), dimension(3) :: int_nablaG_Rankine
    integer                      :: I, J
    logical                      :: derivative_with_respect_to_first_variable

    derivative_with_respect_to_first_variable = adjoint_double_layer

    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) &
    !$OMP&  PRIVATE(J, I, int_G_Rankine, int_nablaG_Rankine)
    do J = 1, nb_faces
      do I = 1, nb_collocation_points
        call integral_of_Rankine(                    &
          collocation_points(I, :),                  &
          vertices(faces(J, :), :),                  &
          centers(J, :),                             &
          normals(J, :),                             &
          areas(J),                                  &
          radiuses(J),                               &
          derivative_with_respect_to_first_variable, &
          int_G_Rankine, int_nablaG_Rankine          &
          )

        S(I, J) = S(I, J) + MINUS_ONE_OVER_FOURPI * int_G_Rankine

        if (size(K, 3) == 1) then  ! early_dot_product=True
          if (adjoint_double_layer) then
            K(I, J, 1) = K(I, J, 1) + MINUS_ONE_OVER_FOURPI * dot_product(dot_product_normals(I, :), int_nablaG_Rankine(:))
          else
            K(I, J, 1) = K(I, J, 1) + MINUS_ONE_OVER_FOURPI * dot_product(dot_product_normals(J, :), int_nablaG_Rankine(:))
          endif
        else
          K(I, J, :) = K(I, J, :) + MINUS_ONE_OVER_FOURPI * int_nablaG_Rankine(:)
        endif
      enddo
    enddo
  end subroutine

end module matrices
