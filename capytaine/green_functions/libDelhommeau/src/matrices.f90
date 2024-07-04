! Copyright (C) 2017-2024 Matthieu Ancellin
! See LICENSE file at <https://github.com/capytaine/libDelhommeau>

MODULE MATRICES

  USE FLOATING_POINT_PRECISION, ONLY: PRE
  USE CONSTANTS

  USE GREEN_RANKINE
  USE GREEN_WAVE

  IMPLICIT NONE

CONTAINS

  ! =====================================================================

  SUBROUTINE BUILD_MATRICES(                             &
      nb_faces_1, centers_1, dot_product_normals,        &
      nb_vertices_2, nb_faces_2, vertices_2, faces_2,    &
      centers_2, normals_2, areas_2, radiuses_2,         &
      nb_quad_points, quad_points, quad_weights,         &
      wavenumber, depth,                                 &
      coeffs,                                            &
      tabulation_nb_integration_points,                  &
      tabulation_grid_shape,                             &
      tabulated_r_range, tabulated_z_range,              &
      tabulated_integrals,                               &
      NEXP, AMBDA, AR,                                   &
      same_body, gf_singularities, adjoint_double_layer, &
      S, K)

    ! Mesh data
    INTEGER,                                     INTENT(IN) :: nb_faces_1, nb_faces_2, nb_vertices_2
    REAL(KIND=PRE), DIMENSION(nb_faces_1, 3),    INTENT(IN) :: centers_1
    REAL(KIND=PRE), DIMENSION(:, :),             INTENT(IN) :: dot_product_normals
    REAL(KIND=PRE), DIMENSION(nb_vertices_2, 3), INTENT(IN) :: vertices_2
    INTEGER,        DIMENSION(nb_faces_2, 4),    INTENT(IN) :: faces_2
    REAL(KIND=PRE), DIMENSION(nb_faces_2, 3),    INTENT(IN) :: centers_2, normals_2
    REAL(KIND=PRE), DIMENSION(nb_faces_2),       INTENT(IN) :: areas_2, radiuses_2

    ! dot_product_normals might be identical to normals_2, especially when adjoint_double_layer is False.
    ! The former is used when computing the dot product with the normal vector in the double layer or adjoint double layer operator
    ! (D or K matrices). The latter is only used when computing the exact formula to integrate the Rankine part of the Green
    ! function over a face.
    ! Hence, both variables fulfill different role, and may or may not be identical.

    INTEGER,                                                  INTENT(IN) :: nb_quad_points
    REAL(KIND=PRE), DIMENSION(nb_faces_2, nb_quad_points, 3), INTENT(IN) :: quad_points
    REAL(KIND=PRE), DIMENSION(nb_faces_2, nb_quad_points),    INTENT(IN) :: quad_weights

    ! Solver parameters
    LOGICAL,                                  INTENT(IN) :: same_body
    INTEGER,                                  INTENT(IN) :: gf_singularities
    LOGICAL,                                  INTENT(IN) :: adjoint_double_layer
    REAL(KIND=PRE), DIMENSION(3)                         :: coeffs

    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber, depth

    ! Tabulated values for the wave part of the Green function
    INTEGER,                                  INTENT(IN) :: tabulation_grid_shape
    INTEGER,                                  INTENT(IN) :: tabulation_nb_integration_points
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_z_range
    REAL(KIND=PRE), DIMENSION(:, :, :),       INTENT(IN) :: tabulated_integrals

    ! Prony decomposition for finite depth Green function
    INTEGER,                                  INTENT(IN) :: NEXP
    REAL(KIND=PRE), DIMENSION(NEXP),          INTENT(IN) :: AMBDA, AR

    ! Outputs
    COMPLEX(KIND=PRE), DIMENSION(:, :), INTENT(INOUT) :: S  ! integrals of the Green function
    COMPLEX(KIND=PRE), DIMENSION(:, :, :), INTENT(INOUT) :: K  ! integrals of the gradient of the Green function

    ! Local variables
    INTEGER                         :: I, J, Q
    REAL(KIND=PRE), DIMENSION(3)    :: reflected_centers_1_I, reflected_int_nablaG_Rankine
    REAL(KIND=PRE)                  :: int_G_Rankine, diagonal_coef
    REAL(KIND=PRE), DIMENSION(3)    :: int_nablaG_Rankine
    COMPLEX(KIND=PRE)               :: int_G, int_G_wave
    COMPLEX(KIND=PRE), DIMENSION(3) :: int_nablaG, int_nablaG_wave, int_nablaG_wave_sym, int_nablaG_wave_antisym
    LOGICAL :: use_symmetry_of_wave_part

    ! use_symmetry_of_wave_part = ((SAME_BODY) .AND. (nb_quad_points == 1))
    use_symmetry_of_wave_part = .false.

    coeffs(:) = coeffs(:)/(-4*PI)  ! Factored out coefficient

    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) &
    !$OMP&  PRIVATE(J, I, int_G, int_nablaG, int_G_Rankine, int_nablaG_Rankine, diagonal_coef, &
    !$OMP&          int_G_wave, int_nablaG_wave, int_nablaG_wave_sym, int_nablaG_wave_antisym, &
    !$OMP&          reflected_centers_1_I, reflected_int_nablaG_Rankine)
    DO J = 1, nb_faces_2
      DO I = 1, nb_faces_1

        int_G = CZERO
        int_nablaG(:) = CZERO

        !!!!!!!!!!!!!!!!!!!
        !  Diagonal term  !
        !!!!!!!!!!!!!!!!!!!
        IF (I==J .and. same_body) THEN

          if (abs(centers_2(j, 3)) < 1e-8) then  ! Panel on the free surface
            diagonal_coef = ONE
          else
            diagonal_coef = ONE/2
          endif

          if (adjoint_double_layer) then
            int_nablaG(:) = int_nablaG(:) + diagonal_coef * dot_product_normals(I, :)
          else
            int_nablaG(:) = int_nablaG(:) + diagonal_coef * dot_product_normals(J, :)
          endif
        ENDIF

        !!!!!!!!!!!!!!!!!!
        !  Rankine part  !
        !!!!!!!!!!!!!!!!!!
        IF (coeffs(1) .NE. ZERO) THEN

          CALL COMPUTE_INTEGRAL_OF_RANKINE_SOURCE( &
            centers_1(I, :),                       &
            vertices_2(faces_2(J, :), :),          &
            centers_2(J, :),                       &
            normals_2(J, :),                       &
            areas_2(J),                            &
            radiuses_2(J),                         &
            int_G_Rankine, int_nablaG_Rankine      &
            )

          int_G = int_G + coeffs(1) * int_G_Rankine

          IF (adjoint_double_layer) THEN
            ! The gradient is with respect to the point I.
            int_nablaG(:) = int_nablaG(:) + coeffs(1) * int_nablaG_Rankine(:)
          ELSE
            ! The gradient is with respect to the panel J.
            ! Due to the symmetry of the Green function, it is just the opposite.
            int_nablaG(:) = int_nablaG(:) - coeffs(1) * int_nablaG_Rankine(:)
          END IF

        END IF


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Reflected Rankine part  !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF ((coeffs(2) .NE. ZERO) .or. &
            ((gf_singularities == LOW_FREQ_WITH_RANKINE_PART) .and. (coeffs(3) .NE. ZERO))) then

          IF (is_infinity(depth)) THEN
            ! Reflection through free surface
            reflected_centers_1_I(1:2) = centers_1(I, 1:2)
            reflected_centers_1_I(3)   = -centers_1(I, 3)
          ELSE
            ! Reflection through sea bottom
            reflected_centers_1_I(1:2) = centers_1(I, 1:2)
            reflected_centers_1_I(3)   = -centers_1(I, 3) - 2*depth
          END IF

          CALL COMPUTE_INTEGRAL_OF_RANKINE_SOURCE( &
            reflected_centers_1_I(:),              &
            vertices_2(faces_2(J, :), :),          &
            centers_2(J, :),                       &
            normals_2(J, :),                       &
            areas_2(J),                            &
            radiuses_2(J),                         &
            int_G_Rankine, int_nablaG_Rankine      &
            )

          reflected_int_nablaG_Rankine(1:2) = int_nablaG_Rankine(1:2)
          reflected_int_nablaG_Rankine(3) = -int_nablaG_Rankine(3)

          int_G = int_G + coeffs(2) * int_G_Rankine

          IF (adjoint_double_layer) THEN
            int_nablaG(1:2) = int_nablaG(1:2) + coeffs(2) * reflected_int_nablaG_Rankine(1:2)
          ELSE
            int_nablaG(1:2) = int_nablaG(1:2) - coeffs(2) * reflected_int_nablaG_Rankine(1:2)
          END IF
          int_nablaG(3) = int_nablaG(3) + coeffs(2) * reflected_int_nablaG_Rankine(3)

          if (gf_singularities == LOW_FREQ_WITH_RANKINE_PART) then
            int_nablaG(3) = int_nablaG(3) + coeffs(3) * 2*wavenumber * int_G_Rankine
          endif
        endif

        !!!!!!!!!!!!!!!
        !  Wave part  !
        !!!!!!!!!!!!!!!
        IF ((coeffs(3) .NE. ZERO) .AND. (.NOT. use_symmetry_of_wave_part)) THEN

          call INTEGRAL_OF_WAVE_PART(                                    &
            centers_1(I, :),                                             &
            centers_2(J, :), areas_2(J),                                 &
            quad_points(J, :, :), quad_weights(J, :),                    &
            wavenumber, depth,                                           &
            tabulation_nb_integration_points, tabulation_grid_shape,     &
            tabulated_r_range, tabulated_z_range, tabulated_integrals,   &
            gf_singularities,                                            &
            NEXP, AMBDA, AR,                                             &
            int_G_wave, int_nablaG_wave_sym, int_nablaG_wave_antisym     &
          )

          int_G = int_G + coeffs(3) * int_G_wave

          IF (adjoint_double_layer) THEN
            int_nablaG(:) = int_nablaG(:) + coeffs(3) * (int_nablaG_wave_sym(:) + int_nablaG_wave_antisym(:))
          ELSE
            int_nablaG(:) = int_nablaG(:) + coeffs(3) * (int_nablaG_wave_sym(:) - int_nablaG_wave_antisym(:))
          END IF

        END IF

        !!!!!!!!!!!!!!!!!!!
        !  Add to matrix  !
        !!!!!!!!!!!!!!!!!!!
        S(I, J) = int_G

        if (size(K, 3) == 1) then  ! early_dot_product=True
          if (adjoint_double_layer) then
            K(I, J, 1) = DOT_PRODUCT(dot_product_normals(I, :), int_nablaG(:))
          else
            K(I, J, 1) = DOT_PRODUCT(dot_product_normals(J, :), int_nablaG(:))
          endif
        else
          K(I, J, :) = int_nablaG(:)
        endif

      END DO  ! loop on I
    END DO  ! parallelized loop on J


    IF ((coeffs(3) .NE. ZERO) .AND. use_symmetry_of_wave_part) THEN
      ! If we are computing the influence of some cells upon themselves, the resulting matrices have some symmetries.
      ! This is due to the symmetry of the Green function, and the way the integral on the face is approximated.
      ! (More precisely, the Green function is symmetric and its derivative is the sum of a symmetric part and an anti-symmetric
      ! part.)

      !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(J, I, int_G_wave, int_nablaG_wave, int_nablaG_wave_sym, int_nablaG_wave_antisym)
      DO J = 1, nb_faces_2
        DO I = J, nb_faces_1

          call INTEGRAL_OF_WAVE_PART(                                    &
            centers_1(I, :),                                             &
            centers_2(J, :), areas_2(J),                                 &
            quad_points(J, :, :), quad_weights(J, :),                    &
            wavenumber, depth,                                           &
            tabulation_nb_integration_points, tabulation_grid_shape,     &
            tabulated_r_range, tabulated_z_range, tabulated_integrals,   &
            gf_singularities,                                            &
            NEXP, AMBDA, AR,                                             &
            int_G_wave, int_nablaG_wave_sym, int_nablaG_wave_antisym     &
          )

          S(I, J) = S(I, J) + coeffs(3) * int_G_wave

          IF (adjoint_double_layer) THEN
            int_nablaG_wave(:) = int_nablaG_wave_sym(:) + int_nablaG_wave_antisym(:)
          ELSE
            int_nablaG_wave(:) = int_nablaG_wave_sym(:) - int_nablaG_wave_antisym(:)
          END IF

          if (size(K, 3) == 1) then
            if (.NOT. adjoint_double_layer) then
              K(I, J, 1) = K(I, J, 1) + coeffs(3) * DOT_PRODUCT(dot_product_normals(J, :), int_nablaG_wave(:))
            else
              K(I, J, 1) = K(I, J, 1) + coeffs(3) * DOT_PRODUCT(dot_product_normals(I, :), int_nablaG_wave(:))
            endif
          else
            K(I, J, :) = K(I, J, :) + coeffs(3) * int_nablaG_wave(:)
          endif

          IF (.NOT. I==J) THEN

            IF (.NOT. adjoint_double_layer) THEN
              int_nablaG_wave(:) = int_nablaG_wave_sym(:) + int_nablaG_wave_antisym(:)
            ELSE
              int_nablaG_wave(:) = int_nablaG_wave_sym(:) - int_nablaG_wave_antisym(:)
            END IF

            S(J, I) = S(J, I) + coeffs(3) * int_G_wave * areas_2(I)/areas_2(J)
            if (size(K, 3) == 1) then
              if (.NOT. adjoint_double_layer) then
                K(J, I, 1) = K(J, I, 1) + coeffs(3) * DOT_PRODUCT(dot_product_normals(I, :), int_nablaG_wave(:)) * &
                  areas_2(I)/areas_2(J)
              else
                K(J, I, 1) = K(J, I, 1) + coeffs(3) * DOT_PRODUCT(dot_product_normals(J, :), int_nablaG_wave(:)) * &
                  areas_2(I)/areas_2(J)
              endif
            else
              K(J, I, :) = K(J, I, :) + coeffs(3) * int_nablaG_wave(:) * areas_2(I)/areas_2(J)
            endif
          END IF
        END DO
      END DO
    END IF

  END SUBROUTINE

END MODULE MATRICES
