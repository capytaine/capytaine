MODULE MATRICES

  USE FLOATING_POINT_PRECISION, ONLY: PRE
  USE CONSTANTS

  USE GREEN_RANKINE
  USE GREEN_WAVE

  IMPLICIT NONE

CONTAINS

  ! =====================================================================

  SUBROUTINE BUILD_MATRICES(                          &
      nb_faces_1, centers_1, normals_1,               &
      nb_vertices_2, nb_faces_2, vertices_2, faces_2, &
      centers_2, normals_2, areas_2, radiuses_2,      &
      nb_quad_points, quad_points, quad_weights,      &
      wavenumber, depth,                              &
      coeffs,                                         &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      NEXP, AMBDA, AR,                                &
      same_body,                                      &
      S, K)

    ! Mesh data
    INTEGER,                                     INTENT(IN) :: nb_faces_1, nb_faces_2, nb_vertices_2
    REAL(KIND=PRE), DIMENSION(nb_faces_1, 3),    INTENT(IN) :: centers_1, normals_1
    REAL(KIND=PRE), DIMENSION(nb_vertices_2, 3), INTENT(IN) :: vertices_2
    INTEGER,        DIMENSION(nb_faces_2, 4),    INTENT(IN) :: faces_2
    REAL(KIND=PRE), DIMENSION(nb_faces_2, 3),    INTENT(IN) :: centers_2, normals_2
    REAL(KIND=PRE), DIMENSION(nb_faces_2),       INTENT(IN) :: areas_2, radiuses_2

    INTEGER,                                                  INTENT(IN) :: nb_quad_points
    REAL(KIND=PRE), DIMENSION(nb_faces_2, nb_quad_points, 3), INTENT(IN) :: quad_points
    REAL(KIND=PRE), DIMENSION(nb_faces_2, nb_quad_points),    INTENT(IN) :: quad_weights

    LOGICAL,                                  INTENT(IN) :: same_body

    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber, depth

    REAL(KIND=PRE), DIMENSION(3)                         :: coeffs

    ! Tabulated data
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_z_range
    REAL(KIND=PRE), DIMENSION(:, :, :, :),    INTENT(IN) :: tabulated_integrals

    ! Prony decomposition for finite depth
    INTEGER,                                  INTENT(IN) :: NEXP
    REAL(KIND=PRE), DIMENSION(NEXP),          INTENT(IN) :: AMBDA, AR

    ! Output
    COMPLEX(KIND=PRE), DIMENSION(:, :), INTENT(INOUT) :: S
    COMPLEX(KIND=PRE), DIMENSION(:, :, :), INTENT(INOUT) :: K

    ! Local variables
    INTEGER                         :: I, J, Q
    REAL(KIND=PRE), DIMENSION(3)    :: reflected_centers_1_I, reflected_VSP1
    REAL(KIND=PRE)                  :: SP1
    REAL(KIND=PRE), DIMENSION(3)    :: VSP1
    COMPLEX(KIND=PRE)               :: SP2
    COMPLEX(KIND=PRE), DIMENSION(3) :: VSP2_SYM, VSP2_ANTISYM
    LOGICAL :: use_symmetry_of_wave_part

    use_symmetry_of_wave_part = ((SAME_BODY) .AND. (nb_quad_points == 1))


#ifndef XIE_CORRECTION
    IF (is_infinity(depth)) THEN
      ! In infinite depth, for finite frequency,
      ! the original Delhommeau's method as in Nemoh uses coeffs(2) = -1.
      ! But by default, coeffs(2) is set to 1 in delhommeau.py.
      coeffs(2) = coeffs(2) - 2*coeffs(3)
    ENDIF
#endif

    coeffs(:) = coeffs(:)/(4*PI)


    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) &
    !$OMP&  PRIVATE(J, I, SP1, VSP1, SP2, VSP2_SYM, VSP2_ANTISYM, reflected_centers_1_I, reflected_VSP1)
    DO J = 1, nb_faces_2

      !!!!!!!!!!!!!!!!!!!!
      !  Initialization  !
      !!!!!!!!!!!!!!!!!!!!
      S(:, J) = CMPLX(0.0, 0.0, KIND=PRE)
      K(:, J, :) = CMPLX(0.0, 0.0, KIND=PRE)

      !!!!!!!!!!!!!!!!!!
      !  Rankine part  !
      !!!!!!!!!!!!!!!!!!
      IF (coeffs(1) .NE. ZERO) THEN
        DO I = 1, nb_faces_1

          CALL COMPUTE_INTEGRAL_OF_RANKINE_SOURCE( &
            centers_1(I, :),                       &
            vertices_2(faces_2(J, :), :),          &
            centers_2(J, :),                       &
            normals_2(J, :),                       &
            areas_2(J),                            &
            radiuses_2(J),                         &
            SP1, VSP1                              &
            )

          ! Store into influence matrix
          S(I, J) = S(I, J) - coeffs(1) * SP1                                ! Green function
          if (size(K, 3) == 1) then
            K(I, J, 1) = K(I, J, 1) - coeffs(1) * DOT_PRODUCT(normals_1(I, :), VSP1(:))
          else
            K(I, J, :) = K(I, J, :) - coeffs(1) * VSP1(:)
          endif

        END DO
      END IF


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Reflected Rankine part  !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (coeffs(2) .NE. ZERO) THEN
        DO I = 1, nb_faces_1

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
            reflected_centers_1_I(:),                &
            vertices_2(faces_2(J, :), :),          &
            centers_2(J, :),                       &
            normals_2(J, :),                       &
            areas_2(J),                            &
            radiuses_2(J),                         &
            SP1, VSP1                              &
            )

          reflected_VSP1(1:2) = VSP1(1:2)
          reflected_VSP1(3) = -VSP1(3)

          ! Store into influence matrix
          S(I, J) = S(I, J) - coeffs(2) * SP1                                ! Green function
          if (size(K, 3) == 1) then
            K(I, J, 1) = K(I, J, 1) - coeffs(2) * DOT_PRODUCT(normals_1(I, :), reflected_VSP1(:))
          else
            K(I, J, :) = K(I, J, :) - coeffs(2) * reflected_VSP1(:)
          endif
        END DO
      END IF

      !!!!!!!!!!!!!!!
      !  Wave part  !
      !!!!!!!!!!!!!!!

      IF ((coeffs(3) .NE. ZERO) .AND. (.NOT. use_symmetry_of_wave_part)) THEN
        DO I = 1, nb_faces_1

          call INTEGRAL_OF_WAVE_PART(                                    &
            centers_1(I, :),                                             &
            centers_2(J, :), normals_2(J, :), areas_2(J), radiuses_2(J), &
            quad_points(J, :, :), quad_weights(J, :),                    &
            wavenumber, depth,                                           &
            tabulated_r_range, tabulated_z_range, tabulated_integrals,   &
            NEXP, AMBDA, AR,                                             &
            SP2, VSP2_SYM, VSP2_ANTISYM                                  &
          )

          S(I, J) = S(I, J) - coeffs(3) * SP2

          if (size(K, 3) == 1) then
            K(I, J, 1) = K(I, J, 1) - coeffs(3) * &
              DOT_PRODUCT(normals_1(I, :), VSP2_SYM + VSP2_ANTISYM)
          else
            K(I, J, :) = K(I, J, :) - coeffs(3) * (VSP2_SYM + VSP2_ANTISYM)
          endif

        END DO
      END IF

      IF (SAME_BODY) THEN
        if (abs(centers_1(j, 3)) < 1e-8) then  ! Panel on the free surface
          if (size(K, 3) == 1) then
            K(J, J, 1) = K(J, J, 1) - 0.25
          else
            K(J, J, :) = K(J, J, :) - 0.25 * normals_1(J, :)
          endif
        else
          if (size(K, 3) == 1) then
            K(J, J, 1) = K(J, J, 1) + 0.5
          else
            K(J, J, :) = K(J, J, :) + 0.5 * normals_1(J, :)
          endif
        endif
      END IF

    END DO  ! parallelized loop on J


    IF ((coeffs(3) .NE. ZERO) .AND. use_symmetry_of_wave_part) THEN
      ! If we are computing the influence of some cells upon themselves, the resulting matrices have some symmetries.
      ! This is due to the symmetry of the Green function, and the way the integral on the face is approximated.
      ! (More precisely, the Green function is symmetric and its derivative is the sum of a symmetric part and an anti-symmetric
      ! part.)

      !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(J, I, SP2, VSP2_SYM, VSP2_ANTISYM)
      DO J = 1, nb_faces_2
        DO I = J, nb_faces_1

          call INTEGRAL_OF_WAVE_PART(                                    &
            centers_1(I, :),                                             &
            centers_2(J, :), normals_2(J, :), areas_2(J), radiuses_2(J), &
            quad_points(J, :, :), quad_weights(J, :),                    &
            wavenumber, depth,                                           &
            tabulated_r_range, tabulated_z_range, tabulated_integrals,   &
            NEXP, AMBDA, AR,                                             &
            SP2, VSP2_SYM, VSP2_ANTISYM                                  &
          )

          S(I, J) = S(I, J) - coeffs(3) * SP2
          if (size(K, 3) == 1) then
            K(I, J, 1) = K(I, J, 1) - coeffs(3) * &
              DOT_PRODUCT(normals_1(I, :), VSP2_SYM + VSP2_ANTISYM)
          else
            K(I, J, :) = K(I, J, :) - coeffs(3) * (VSP2_SYM + VSP2_ANTISYM)
          endif

          IF (.NOT. I==J) THEN
            S(J, I) = S(J, I) - coeffs(3) * SP2 * areas_2(I)/areas_2(J)
            if (size(K, 3) == 1) then
              K(J, I, 1) = K(J, I, 1) - coeffs(3) * &
                DOT_PRODUCT(normals_1(J, :), VSP2_SYM - VSP2_ANTISYM) * areas_2(I)/areas_2(J)
            else
              K(J, I, :) = K(J, I, :) - coeffs(3) * (VSP2_SYM - VSP2_ANTISYM) * areas_2(I)/areas_2(J)
            endif
          END IF
        END DO
      END DO
    END IF

  END SUBROUTINE

END MODULE MATRICES
