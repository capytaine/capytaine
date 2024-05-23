! Copyright (C) 2017-2024 Matthieu Ancellin
! See LICENSE file at <https://github.com/capytaine/libDelhommeau>
MODULE GREEN_WAVE

  USE FLOATING_POINT_PRECISION, ONLY: PRE
  USE CONSTANTS
  USE DELHOMMEAU_INTEGRALS
  USE GREEN_RANKINE, ONLY: COMPUTE_ASYMPTOTIC_RANKINE_SOURCE

  IMPLICIT NONE

  ! Dependencies between the functions of this module:
  ! (from top to bottom: "is called by")
  !
  !                                            (Delhommeau_integrals.f90)
  !                                                       |
  !                                           WAVE_PART_INFINITE_DEPTH           (COMPUTE_ASYMPTOTIC_RANKINE_SOURCE)
  !                                                     /   \                    /
  !          INTEGRAL_OF_SINGULARITY_ON_FREE_SURFACE   |    WAVE_PART_FINITE_DEPTH
  !                                                \   \   /
  !                                              INTEGRAL_OF_WAVE_PART
  !                                                       |
  !                                                 (matrices.f90)
  !                                                       |
  !                                                 (python code)

CONTAINS

  ! =====================================================================

  subroutine integral_of_wave_part                               &
      (x,                                                        &
      face_center, face_area,                                    &
      face_quadrature_points, face_quadrature_weights,           &
      wavenumber, depth,                                         &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      nexp, ambda, ar,                                           &
      int_G, int_nablaG_sym, int_nablaG_antisym                  &
      )
    ! Integral over a panel of the wave part of the Green function.

    real(kind=pre), dimension(3),          intent(in) :: x
    real(kind=pre), dimension(3),          intent(in) :: face_center
    real(kind=pre),                        intent(in) :: face_area
    real(kind=pre), dimension(:),          intent(in) :: face_quadrature_weights
    real(kind=pre), dimension(:, :),       intent(in) :: face_quadrature_points
    real(kind=pre),                        intent(in) :: wavenumber, depth
    integer,                               intent(in) :: gf_singularities
    integer,                               intent(in) :: tabulation_nb_integration_points
    integer,                               intent(in) :: tabulation_grid_shape
    real(kind=pre), dimension(:),          intent(in) :: tabulated_r_range
    real(kind=pre), dimension(:),          intent(in) :: tabulated_z_range
    real(kind=pre), dimension(:, :, :),    intent(in) :: tabulated_integrals
    integer,                               intent(in) :: nexp
    real(kind=pre), dimension(nexp),       intent(in) :: ambda, ar

    complex(kind=pre),                     intent(out) :: int_G
    complex(kind=pre), dimension(3),       intent(out) :: int_nablaG_sym, int_nablaG_antisym

    ! Local variables
    real(kind=pre)                  :: r, z
    complex(kind=pre)               :: G_at_point
    complex(kind=pre), dimension(3) :: nablaG_at_point, nablaG_at_point_sym, nablaG_at_point_antisym
    integer                         :: nb_quad_points, Q

    nb_quad_points = size(face_quadrature_weights)

    r = wavenumber * norm2(x(1:2) - face_center(1:2))
    z = wavenumber * (x(3) + face_center(3))

    if ((abs(r) < 1e-10) .and. (abs(z) < 1e-10)) then
      if (gf_singularities == HIGH_FREQ) then
        print*, "WARNING: support of free surface panels not implemented for the high_freq Green function"
      endif

      ! Interaction of a panel on the free surface with itself
      call integral_of_singularity_on_free_surface( &
        face_area, wavenumber, &
        int_G, int_nablaG_sym, int_nablaG_antisym &
        )

    else
      ! Numerical integration
      int_G = zero
      int_nablaG_sym = zero
      int_nablaG_antisym = zero

      do q = 1, nb_quad_points
        if (is_infinity(depth)) then
          call wave_part_infinite_depth                                &
            (x,                                                        &
            face_quadrature_points(q, :),                              &
            wavenumber,                                                &
            tabulation_nb_integration_points, tabulation_grid_shape,   &
            tabulated_r_range, tabulated_z_range, tabulated_integrals, &
            gf_singularities,                                          &
            G_at_point, nablaG_at_point                                &
            )
            nablaG_at_point_sym(1:2) = cmplx(zero, zero, kind=pre)
            nablaG_at_point_sym(3) = nablaG_at_point(3)
            nablaG_at_point_antisym(1:2) = nablaG_at_point(1:2)
            nablaG_at_point_antisym(3) = cmplx(zero, zero, kind=pre)
        else
          call wave_part_finite_depth                                  &
            (x,                                                        &
            face_quadrature_points(q, :),                              &
            wavenumber,                                                &
            depth,                                                     &
            tabulation_nb_integration_points, tabulation_grid_shape,   &
            tabulated_r_range, tabulated_z_range, tabulated_integrals, &
            nexp, ambda, ar,                                           &
            G_at_point, nablaG_at_point_sym, nablaG_at_point_antisym   &
            )
        end if

        int_G = int_G + G_at_point * face_quadrature_weights(q)
        int_nablaG_sym = int_nablaG_sym + nablaG_at_point_sym * face_quadrature_weights(q)
        int_nablaG_antisym = int_nablaG_antisym + nablaG_at_point_antisym * face_quadrature_weights(q)
      end do
    end if

  end subroutine

  ! =====================================================================

  subroutine integral_of_singularity_on_free_surface(face_area, wavenumber, int_G, int_nablaG_sym, int_nablaG_antisym)
  ! Integrating the wave term by approximating the panel by a circle of same area.
  ! The singularities are integrated analytically. The rest is integrated with a 1-point integral.

  ! G_w^+ = - 2 k log(k r) + rest
  ! ∫_Γ G_w^+ dξ ∼ - 2 k ∫_Γ log(k r) dξ + |Γ| rest(0)
  ! with
  ! ∫_Γ log(k r) dξ = |Γ|/2 (log(k^2 Γ/π) - 1)
  ! and
  ! rest(0) = 2 k ( γ - log(2) ) - 2ikπ
  !
  ! Also
  ! dG_w^+/dz = G_w^+ + 2/r
  ! with
  ! ∫_Γ 1/r dξ ∼ 2 √( π |Γ| )
  ! TODO: replace the latter with the actual integral since we computed it anyway for the Rankine term.

  ! TODO: only low_freq singularities are implemented here
    real(kind=pre),                        intent(in) :: wavenumber
    real(kind=pre),                        intent(in) :: face_area

    complex(kind=pre),                     intent(out) :: int_G
    complex(kind=pre), dimension(3),       intent(out) :: int_nablaG_sym, int_nablaG_antisym

    int_G = wavenumber * face_area * (                      &
                  (1 - log(wavenumber**2 * face_area / pi)) &
                  + 2 * (euler_gamma - log_2) - 2*pi*ii     &
                )
    int_nablaG_sym(1:2) = cmplx(zero, zero, kind=pre)  ! TODO: might be needed for velocity reconstruction.
    int_nablaG_sym(3) = wavenumber * (int_G + 4 * sqrt(pi*face_area))
    int_nablaG_antisym(1:3) = cmplx(zero, zero, kind=pre)  ! Irrelevant anyway because we are on the diagonal
  end subroutine

  ! =====================================================================

  SUBROUTINE WAVE_PART_INFINITE_DEPTH                            &
      (X0I, X0J, wavenumber,                                     &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      G, nablaG)
      ! Returns
      ! (G^-, nabla G^-)                    if gf_singularities == HIGH_FREQ
      ! (G^+, nabla G^+)                    if gf_singularities == LOW_FREQ
      ! (G^+, nabla G^+ - (0, 0, 2*k**2/r1) if gf_singularities == LOW_FREQ_WITH_RANKINE_PART
      ! In the last case, the missing term will be computed with the reflected Rankine term

    ! Inputs
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I, X0J
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber
    integer,                                  intent(in) :: gf_singularities

    ! Tabulated data
    integer,                                  intent(in) :: tabulation_nb_integration_points
    INTEGER,                                  INTENT(IN) :: tabulation_grid_shape
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_z_range
    REAL(KIND=PRE), DIMENSION(:, :, :),       INTENT(IN) :: tabulated_integrals

    ! Outputs
    COMPLEX(KIND=PRE),                        INTENT(OUT) :: G
    COMPLEX(KIND=PRE), DIMENSION(3),          INTENT(OUT) :: nablaG

    ! Local variables
    REAL(KIND=PRE) :: r, r1, z, drdx1, drdx2, dzdx3
    complex(kind=pre) :: dGdr
    REAL(KIND=PRE), dimension(nb_tabulated_values) :: integrals

    r = wavenumber * NORM2(X0I(1:2) - X0J(1:2))
    z = wavenumber * (X0I(3) + X0J(3))
    r1 = hypot(r, z)

    !=======================================================
    ! Evaluate the elementary integrals depending on z and r
    !=======================================================
    IF ((size(tabulated_z_range) <= 1) .or. (size(tabulated_r_range) <= 1)) THEN
      ! No tabulation, fully recompute the Green function each time.
      integrals = numerical_integration(r, z, tabulation_nb_integration_points)
    ELSE
      IF ((abs(z) < abs(tabulated_z_range(size(tabulated_z_range)))) &
          .AND. (r < tabulated_r_range(size(tabulated_r_range)))) THEN
        ! Within the range of tabulated data
        integrals = pick_in_default_tabulation( &
            r, z, tabulation_grid_shape, tabulated_r_range, tabulated_z_range, tabulated_integrals &
        )
      ELSE
        ! Delhommeau's asymptotic expression of Green function for distant panels
        integrals = asymptotic_approximations(MAX(r, 1e-10), z)
      ENDIF
    ENDIF

    !===================================================
    ! Add the elementary integrals to build G and nablaG
    !===================================================

    IF (ABS(r) > 16*EPSILON(r)) THEN
      drdx1 = wavenumber**2 * (X0I(1) - X0J(1))/r
      drdx2 = wavenumber**2 * (X0I(2) - X0J(2))/r
    ELSE
      ! Limit when r->0 is not well defined...
      drdx1 = ZERO
      drdx2 = ZERO
    END IF
    dzdx3 = wavenumber

    if ((gf_singularities == LOW_FREQ) .or. (gf_singularities == LOW_FREQ_WITH_RANKINE_PART)) then
      ! G is G^+, nablaG is nablaG^+
      G = CMPLX(integrals(1), integrals(3), KIND=PRE)
      dGdr = CMPLX(integrals(4), integrals(5), KIND=PRE)
      nablaG(1) = drdx1 * dGdr
      nablaG(2) = drdx2 * dGdr
      if (gf_singularities == LOW_FREQ) then
        nablaG(3) = dzdx3 * (G + 2/r1)
      else
        nablaG(3) = dzdx3 * G
        ! The missing Rankine term is integrated as a Rankine term
      endif
    else if (gf_singularities == HIGH_FREQ) then
      ! G is G^-, nablaG is nablaG^-
      G = CMPLX(integrals(2), integrals(3), KIND=PRE)
      dGdr = CMPLX(integrals(4), integrals(5), KIND=PRE)
      nablaG(1) = drdx1 * dGdr - 2*drdx1*r/r1**3
      nablaG(2) = drdx2 * dGdr - 2*drdx2*r/r1**3
      nablaG(3) = dzdx3 * G    - 2*dzdx3*z/r1**3
    endif

    G = wavenumber * G
    nablaG = wavenumber * nablaG

  END SUBROUTINE WAVE_PART_INFINITE_DEPTH

  ! =====================================================================

  SUBROUTINE WAVE_PART_FINITE_DEPTH                              &
      (X0I, X0J, wavenumber, depth,                              &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      NEXP, AMBDA, AR,                                           &
      SP, VSP_SYM, VSP_ANTISYM)
    ! Compute the frequency-dependent part of the Green function in the finite depth case.

    ! Inputs
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber, depth
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I  ! Coordinates of the source point
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0J  ! Coordinates of the center of the integration panel

    ! Tabulated data
    integer,                                  intent(in) :: tabulation_nb_integration_points
    INTEGER,                                  INTENT(IN) :: tabulation_grid_shape
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_z_range
    REAL(KIND=PRE), DIMENSION(:, :, :),       INTENT(IN) :: tabulated_integrals

    ! Prony decomposition for finite depth
    INTEGER,                                  INTENT(IN) :: NEXP
    REAL(KIND=PRE), DIMENSION(NEXP),          INTENT(IN) :: AMBDA, AR

    ! Outputs
    COMPLEX(KIND=PRE),               INTENT(OUT) :: SP  ! Integral of the Green function over the panel.
    COMPLEX(KIND=PRE), DIMENSION(3), INTENT(OUT) :: VSP_SYM, VSP_ANTISYM ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    INTEGER                              :: KE
    REAL(KIND=PRE)                       :: AMH, AKH, A
    REAL(KIND=PRE)                       :: AQT
    REAL(KIND=PRE),    DIMENSION(3)      :: XI, XJ
    REAL(KIND=PRE),    DIMENSION(4)      :: FTS
    REAL(KIND=PRE),    DIMENSION(3, 4)   :: VTS
    COMPLEX(KIND=PRE), DIMENSION(4)      :: FS
    COMPLEX(KIND=PRE), DIMENSION(3, 4)   :: VS

    !========================================
    ! Part 1: Solve 4 infinite depth problems
    !========================================

    XI(:) = X0I(:)
    XJ(:) = X0J(:)

    ! 1.a First infinite depth problem
    CALL WAVE_PART_INFINITE_DEPTH(                               &
      XI(:), XJ(:), wavenumber,                                  &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      LOW_FREQ,                                                  &
      FS(1), VS(:, 1))

    ! 1.b Shift and reflect XI and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) =  X0J(3)
    CALL WAVE_PART_INFINITE_DEPTH(                               &
      XI(:), XJ(:), wavenumber,                                  &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      LOW_FREQ,                                                  &
      FS(2), VS(:, 2))
    VS(3, 2) = -VS(3, 2) ! Reflection of the output vector

    ! 1.c Shift and reflect XJ and compute another value of the Green function
    XI(3) =  X0I(3)
    XJ(3) = -X0J(3) - 2*depth
    CALL WAVE_PART_INFINITE_DEPTH(                               &
      XI(:), XJ(:), wavenumber,                                  &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      LOW_FREQ,                                                  &
      FS(3), VS(:, 3))

    ! 1.d Shift and reflect both XI and XJ and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) = -X0J(3) - 2*depth
    CALL WAVE_PART_INFINITE_DEPTH(                               &
      XI(:), XJ(:), wavenumber,                                  &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      LOW_FREQ,                                                  &
      FS(4), VS(:, 4))
    VS(3, 4) = -VS(3, 4) ! Reflection of the output vector

    ! Add up the results of the four problems

    SP               = SUM(FS(1:4))
    VSP_SYM(1:2)     = CMPLX(ZERO, ZERO, KIND=PRE)
    VSP_ANTISYM(1:2) = VS(1:2, 1) + VS(1:2, 2) + VS(1:2, 3) + VS(1:2, 4)
    VSP_SYM(3)       = VS(3, 1) + VS(3, 4)
    VSP_ANTISYM(3)   = VS(3, 2) + VS(3, 3)

    ! Multiply by some coefficients
    AMH  = wavenumber*depth
    AKH  = AMH*TANH(AMH)
    A    = (AMH+AKH)**2/(4*AMH*(AMH**2-AKH**2+AKH))

    SP          = A*SP
    VSP_ANTISYM = A*VSP_ANTISYM
    VSP_SYM     = A*VSP_SYM

    !=====================================================
    ! Part 2: Integrate (NEXP+1)×4 terms of the form 1/MM'
    !=====================================================

    DO KE = 1, NEXP
      XI(:) = X0I(:)

      ! 2.a Shift observation point and compute integral
      XI(3) =  X0I(3) + depth*AMBDA(KE) - 2*depth
      CALL COMPUTE_ASYMPTOTIC_RANKINE_SOURCE(XI(:), X0J(:), ONE, FTS(1), VTS(:, 1))

      ! 2.b Shift and reflect observation point and compute integral
      XI(3) = -X0I(3) - depth*AMBDA(KE)
      CALL COMPUTE_ASYMPTOTIC_RANKINE_SOURCE(XI(:), X0J(:), ONE, FTS(2), VTS(:, 2))
      VTS(3, 2) = -VTS(3, 2) ! Reflection of the output vector

      ! 2.c Shift and reflect observation point and compute integral
      XI(3) = -X0I(3) + depth*AMBDA(KE) - 4*depth
      CALL COMPUTE_ASYMPTOTIC_RANKINE_SOURCE(XI(:), X0J(:), ONE, FTS(3), VTS(:, 3))
      VTS(3, 3) = -VTS(3, 3) ! Reflection of the output vector

      ! 2.d Shift observation point and compute integral
      XI(3) =  X0I(3) - depth*AMBDA(KE) + 2*depth
      CALL COMPUTE_ASYMPTOTIC_RANKINE_SOURCE(XI(:), X0J(:), ONE, FTS(4), VTS(:, 4))

      AQT = AR(KE)/2

      ! Add all the contributions
      SP               = SP               + AQT*SUM(FTS(1:4))

      VSP_ANTISYM(1:2) = VSP_ANTISYM(1:2) + AQT*(VTS(1:2, 1) + VTS(1:2, 2) + VTS(1:2, 3) + VTS(1:2, 4))
      VSP_ANTISYM(3)   = VSP_ANTISYM(3) + AQT*(VTS(3, 1) + VTS(3, 4))
      VSP_SYM(3)       = VSP_SYM(3) + AQT*(VTS(3, 2) + VTS(3, 3))

    END DO

    RETURN
  END SUBROUTINE

END MODULE GREEN_WAVE
