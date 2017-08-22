SUBROUTINE BUILD_MATRICES(                                          &
    nb_faces_1, centers_1, normals_1,                               &
    nb_vertices_2, nb_faces_2,                                      &
    vertices_2, faces_2, centers_2, normals_2, areas_2, radiuses_2, &
    wavenumber, depth,                                              &
    S, V)

  USE GREEN_1, ONLY: VAV
  USE GREEN_2, ONLY: VNSFD, VNSINFD
  USE INITIALIZE_GREEN_2, ONLY: LISC

  INTEGER,                              INTENT(IN) :: nb_faces_1, nb_faces_2, nb_vertices_2
  REAL,    DIMENSION(nb_faces_1, 3),    INTENT(IN) :: centers_1
  REAL,    DIMENSION(nb_faces_1, 3),    INTENT(IN) :: normals_1
  REAL,    DIMENSION(nb_vertices_2, 3), INTENT(IN) :: vertices_2
  INTEGER, DIMENSION(nb_faces_2, 4),    INTENT(IN) :: faces_2
  REAL,    DIMENSION(nb_faces_2, 3),    INTENT(IN) :: centers_2
  REAL,    DIMENSION(nb_faces_2, 3),    INTENT(IN) :: normals_2
  REAL,    DIMENSION(nb_faces_2),       INTENT(IN) :: areas_2, radiuses_2
  REAL,                                 INTENT(IN) :: wavenumber, depth

  COMPLEX, DIMENSION(nb_faces_1, nb_faces_2), INTENT(OUT) :: S
  COMPLEX, DIMENSION(nb_faces_1, nb_faces_2), INTENT(OUT) :: V

  ! Local variables
  INTEGER :: I, J

  REAL                  :: SP1,  SM1
  REAL, DIMENSION(3)    :: VSP1, VSM1
  COMPLEX               :: SP2,  SM2
  COMPLEX, DIMENSION(3) :: VSP2, VSM2

  DO I = 1, nb_faces_1
    DO J = 1, nb_faces_2

      IF (depth == 0.0) THEN
        CALL VAV                        &
          (centers_1(I, :),             &
          vertices_2(faces_2(J, :), :), &
          centers_2(J, :),              &
          normals_2(J, :),              &
          areas_2(J),                   &
          radiuses_2(J),                &
          depth,                        &
          -1,                           &
          SP1, SM1, VSP1, VSM1          &
          )

        CALL VNSINFD                    &
          (wavenumber,                  &
          centers_1(I, :),              &
          centers_2(J, :),              &
          areas_2(J),                   &
          SP2, SM2, VSP2, VSM2          &
          )
      ELSE
        CALL VAV                        &
          (centers_1(I, :),             &
          vertices_2(faces_2(J, :), :), &
          centers_2(J, :),              &
          normals_2(J, :),              &
          areas_2(J),                   &
          radiuses_2(J),                &
          depth,                        &
          1,                            &
          SP1, SM1, VSP1, VSM1          &
          )

        CALL VNSFD                      &
          (wavenumber,                  &
          centers_1(I, :),              &
          centers_2(J, :),              &
          areas_2(J),                   &
          depth,                        &
          SP2, SM2, VSP2, VSM2          &
          )
      END IF

      ! Store into influence matrix
      S(I, J) = SP1 + SP2                              ! Green function
      V(I, J) = DOT_PRODUCT(normals_1(I, :), VSP1 + VSP2) ! Gradient of the Green function

    END DO
  END DO

END SUBROUTINE BUILD_MATRICES
