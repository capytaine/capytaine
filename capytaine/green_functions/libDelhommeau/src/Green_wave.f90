! Copyright (C) 2017-2025 Matthieu Ancellin
! See LICENSE file at <https://github.com/capytaine/libDelhommeau>
module Green_Wave

  use floating_point_precision, only: pre
  use constants
  use delhommeau_integrals
#ifdef LIANGWUNOBLESSE_OPTIONAL_DEPENDENCY
    use liangwunoblessewaveterm, only: havelockgf
#endif
#ifdef FINGREEN3D_OPTIONAL_DEPENDENCY
    use fingreen3d_module, only: fingreen3d_routine, dispersion
#endif
  use green_rankine, only: integral_of_reflected_rankine

  implicit none

  ! Dependencies between the functions of this module:
  ! (from top to bottom: "is called by")
  !
  !                                            (Delhommeau_integrals.f90 or LiangWuNoblesseWaveTerm.f90)
  !                                                                       |
  !                integral_of_singularity_on_free_surface          wave_part_infinite_depth
  !                                                       \               |
  !                    integral_of_reflected_Rankine      integral_of_wave_part_infinite_depth
  !                                             |         /          |
  !                     integral_of_wave_part_finite_depth           |
  !                                              \                   /
  !                                              integral_of_wave_part
  !                                                       |
  !                                                 (matrices.f90)
  !                                                       |
  !                                                 (python code)

CONTAINS

  ! =====================================================================

  subroutine integral_of_wave_part_infinite_depth                &
      (x,                                                        &
      face_center, face_area,                                    &
      face_quadrature_points, face_quadrature_weights,           &
      wavenumber,                                                &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      derivative_with_respect_to_first_variable,                 &
      int_G, int_nablaG                                          &
      )

    real(kind=pre), dimension(3),          intent(in) :: x
    real(kind=pre), dimension(3),          intent(in) :: face_center
    real(kind=pre),                        intent(in) :: face_area
    real(kind=pre), dimension(:),          intent(in) :: face_quadrature_weights
    real(kind=pre), dimension(:, :),       intent(in) :: face_quadrature_points
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
    real(kind=pre)                  :: r, z
    complex(kind=pre)               :: G_at_point
    complex(kind=pre), dimension(3) :: nablaG_at_point
    integer                         :: nb_quad_points, Q

    r = wavenumber * norm2(x(1:2) - face_center(1:2))
    z = wavenumber * (x(3) + face_center(3))

    if ((abs(r) < 1e-10) .and. (abs(z) < 1e-10)) then
      if (gf_singularities == HIGH_FREQ) then
        print*, "WARNING: support of free surface panels not implemented for the high_freq Green function"
        ! Delhommeau().evaluate() should not permit this to happen.
      endif

      ! Interaction of a panel on the free surface with itself
      call integral_of_singularity_on_free_surface( &
        face_area, wavenumber, int_G, int_nablaG)

    else
      ! Numerical integration
      int_G = zero
      int_nablaG = zero

      nb_quad_points = size(face_quadrature_weights)

      do q = 1, nb_quad_points
        call wave_part_infinite_depth                                &
          (x,                                                        &
          face_quadrature_points(q, :),                              &
          wavenumber,                                                &
          tabulation_nb_integration_points, tabulation_grid_shape,   &
          tabulated_r_range, tabulated_z_range, tabulated_integrals, &
          gf_singularities,                                          &
          G_at_point, nablaG_at_point                                &
          )
        int_G = int_G + G_at_point * face_quadrature_weights(q)
        int_nablaG(:) = int_nablaG(:) + nablaG_at_point(:) * face_quadrature_weights(q)
      end do
    end if

    if (.not. derivative_with_respect_to_first_variable) then
      int_nablaG(1:2) = -int_nablaG(1:2)
    endif

  end subroutine

  ! =====================================================================

  subroutine integral_of_singularity_on_free_surface(face_area, wavenumber, int_G, int_nablaG)
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
    complex(kind=pre), dimension(3),       intent(out) :: int_nablaG

    int_G = wavenumber * face_area * (                      &
                  (1 - log(wavenumber**2 * face_area / pi)) &
                  + 2 * (euler_gamma - log_2) - 2*pi*ii     &
                )
    int_nablaG(1:2) = cmplx(zero, zero, kind=pre)  ! TODO: might be needed for velocity reconstruction.
    int_nablaG(3) = wavenumber * (int_G + 4 * sqrt(pi*face_area))
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
    complex(kind=pre) :: G_, dGdr_, dGdr
    REAL(KIND=PRE), dimension(nb_tabulated_values) :: integrals

    r = wavenumber * NORM2(X0I(1:2) - X0J(1:2))
    z = wavenumber * (X0I(3) + X0J(3))
    r1 = hypot(r, z)

    IF (ABS(r) > 16*EPSILON(r)) THEN
      drdx1 = wavenumber**2 * (X0I(1) - X0J(1))/r
      drdx2 = wavenumber**2 * (X0I(2) - X0J(2))/r
    ELSE
      ! Limit when r->0 is not well defined...
      drdx1 = ZERO
      drdx2 = ZERO
    END IF
    dzdx3 = wavenumber

    IF (tabulation_grid_shape == LIANG_WU_NOBLESSE) THEN
#ifdef LIANGWUNOBLESSE_OPTIONAL_DEPENDENCY
        call HavelockGF(r, z, G_, dGdr_)
        G = - G_
        dGdr = - dGdr_
        nablaG(1) = drdx1 * dGdr
        nablaG(2) = drdx2 * dGdr
        nablaG(3) = dzdx3 * (G + 2/r1)
#endif
    ELSE

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
    ENDIF

    G = wavenumber * G
    nablaG = wavenumber * nablaG

  END SUBROUTINE WAVE_PART_INFINITE_DEPTH

  ! =====================================================================

  pure function symmetric_of_vector(n) result(n_sym)
    complex(kind=pre), dimension(3), intent(in) :: n
    complex(kind=pre), dimension(3) :: n_sym

    n_sym(1:2) = n(1:2)
    n_sym(3)   = -n(3)
  end function

  pure function sea_bottom_symmetric_of_point(x, depth) result(x_sym)
    real(kind=pre), dimension(3), intent(in) :: x
    real(kind=pre),               intent(in) :: depth
    real(kind=pre), dimension(3) :: x_sym

    x_sym(1:2) = x(1:2)
    x_sym(3) = -x(3) - 2*depth
  end function

  pure subroutine sea_bottom_symmetric_of_face( &
    face_center, face_quadrature_points,        &
    depth,                                      &
    face_center_sym, face_quadrature_points_sym &
    )
    real(kind=pre), dimension(3),    intent(in) :: face_center
    real(kind=pre), dimension(:, :), intent(in) :: face_quadrature_points
    real(kind=pre),                  intent(in) :: depth
    real(kind=pre), dimension(3),    intent(out) :: face_center_sym
    real(kind=pre), dimension(:, :), intent(out) :: face_quadrature_points_sym

    integer :: i

    face_center_sym = sea_bottom_symmetric_of_point(face_center, depth)
    do i = 1, size(face_quadrature_points, 1)
      face_quadrature_points_sym(i, :) = sea_bottom_symmetric_of_point(face_quadrature_points(i, :), depth)
    end do
  end subroutine

  subroutine integral_of_wave_parts_finite_depth                 &
      (x,                                                        &
      face_center, face_area,                                    &
      face_quadrature_points, face_quadrature_weights,           &
      wavenumber, depth,                                         &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      derivative_with_respect_to_first_variable,                 &
      int_G, int_nablaG                                          &
      )

    ! Inputs
    real(kind=pre), dimension(3),             intent(in) :: x
    real(kind=pre), dimension(3),             intent(in) :: face_center
    real(kind=pre),                           intent(in) :: face_area
    real(kind=pre), dimension(:),             intent(in) :: face_quadrature_weights
    real(kind=pre), dimension(:, :),          intent(in) :: face_quadrature_points
    real(kind=pre),                           intent(in) :: wavenumber, depth
    logical,                                  intent(in) :: derivative_with_respect_to_first_variable

    ! Tabulated data
    integer,                                  intent(in) :: tabulation_nb_integration_points
    integer,                                  intent(in) :: tabulation_grid_shape
    real(kind=pre), dimension(:),             intent(in) :: tabulated_r_range
    real(kind=pre), dimension(:),             intent(in) :: tabulated_z_range
    real(kind=pre), dimension(:, :, :),       intent(in) :: tabulated_integrals

    integer,                                  intent(in) :: gf_singularities

    ! Outputs
    complex(kind=pre),               intent(out) :: int_G  ! integral of the Green function over the panel.
    complex(kind=pre), dimension(3), intent(out) :: int_nablaG ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    real(kind=pre), dimension(3)    :: x_sym, face_center_sym
    real(kind=pre), dimension(size(face_quadrature_points, 1), size(face_quadrature_points, 2)) :: face_quadrature_points_sym
    real(kind=pre)                  :: amh, akh, a
    complex(kind=pre)               :: int_G_term
    complex(kind=pre), dimension(3) :: int_nablaG_term

    int_G = czero
    int_nablaG = czero

    x_sym = sea_bottom_symmetric_of_point(x, depth)
    call sea_bottom_symmetric_of_face(            &
      face_center, face_quadrature_points, depth, &
      face_center_sym, face_quadrature_points_sym)

    ! 1.a First infinite depth problem
    CALL integral_of_wave_part_infinite_depth                    &
      (x,                                                        &
      face_center, face_area,                                    &
      face_quadrature_points, face_quadrature_weights,           &
      wavenumber,                                                &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      derivative_with_respect_to_first_variable,                 &
      int_G_term, int_nablaG_term                                &
      )
    int_G = int_G + int_G_term
    int_nablaG = int_nablaG + int_nablaG_term

    ! 1.b Reflect X and compute another value of the Green function
    CALL integral_of_wave_part_infinite_depth                    &
      (x_sym,                                                    &
      face_center, face_area,                                    &
      face_quadrature_points, face_quadrature_weights,           &
      wavenumber,                                                &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      derivative_with_respect_to_first_variable,                 &
      int_G_term, int_nablaG_term                                &
      )
    int_G = int_G + int_G_term
    if (derivative_with_respect_to_first_variable) then
      int_nablaG = int_nablaG + symmetric_of_vector(int_nablaG_term)
    else
      int_nablaG = int_nablaG + int_nablaG_term
    endif

    ! 1.c Reflect face and compute another value of the Green function
    CALL integral_of_wave_part_infinite_depth                    &
      (x,                                                        &
      face_center_sym, face_area,                                &
      face_quadrature_points_sym, face_quadrature_weights,       &
      wavenumber,                                                &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      derivative_with_respect_to_first_variable,                 &
      int_G_term, int_nablaG_term                                &
      )
    int_G = int_G + int_G_term
    if (derivative_with_respect_to_first_variable) then
      int_nablaG = int_nablaG + int_nablaG_term
    else
      int_nablaG = int_nablaG + symmetric_of_vector(int_nablaG_term)
    endif

    ! 1.d Reflect both x and face and compute another value of the Green function
    CALL integral_of_wave_part_infinite_depth                    &
      (x_sym,                                                    &
      face_center_sym, face_area,                                &
      face_quadrature_points_sym, face_quadrature_weights,       &
      wavenumber,                                                &
      tabulation_nb_integration_points, tabulation_grid_shape,   &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                          &
      derivative_with_respect_to_first_variable,                 &
      int_G_term, int_nablaG_term                                &
      )
    int_G = int_G + int_G_term
    int_nablaG = int_nablaG + symmetric_of_vector(int_nablaG_term)

    ! Some coefficient
    AMH  = wavenumber*depth
    AKH  = AMH*TANH(AMH)
    A    = (AMH+AKH)**2/(4*AMH*(AMH**2-AKH**2+AKH))

    int_G = A * int_G
    int_nablaG = A * int_nablaG

  end subroutine

  subroutine integral_of_prony_decomp_finite_depth               &
      (x,                                                        &
      face_nodes,                                                &
      face_center, face_normal, face_area, face_radius,          &
      depth,                                                     &
      prony_decomposition,                                       &
      derivative_with_respect_to_first_variable,                 &
      int_G, int_nablaG                                          &
      )

    ! Inputs
    real(kind=pre), dimension(3),             intent(in) :: x
    real(kind=pre), dimension(4, 3),          intent(in) :: face_nodes
    real(kind=pre), dimension(3),             intent(in) :: face_center, face_normal
    real(kind=pre),                           intent(in) :: face_area, face_radius
    real(kind=pre),                           intent(in) :: depth
    logical,                                  intent(in) :: derivative_with_respect_to_first_variable

    real(kind=pre), dimension(:, :),          intent(in) :: prony_decomposition

    ! Outputs
    complex(kind=pre),               intent(out) :: int_G
    complex(kind=pre), dimension(3), intent(out) :: int_nablaG

    ! Local variables
    real(kind=pre)                  :: lambda_k, a_k
    real(kind=pre)                  :: int_G_term_Rankine
    real(kind=pre), dimension(3)    :: int_nablaG_term_Rankine
    integer                         :: ke

    int_G = czero
    int_nablaG = czero

    do ke = 1, size(prony_decomposition, 2)
      lambda_k = prony_decomposition(1, ke)
      a_k = prony_decomposition(2, ke)

      ! 2.a
      call integral_of_reflected_Rankine(          &
        x,                                         &
        face_nodes, face_center, face_normal,      &
        face_area, face_radius,                    &
        derivative_with_respect_to_first_variable, &
        [ONE, depth*lambda_k - 2*depth],           &
        int_G_term_Rankine, int_nablaG_term_Rankine)
      int_G = int_G + a_k/2 * int_G_term_Rankine
      int_nablaG = int_nablaG + a_k/2 * int_nablaG_term_Rankine

      ! 2.b
      call integral_of_reflected_Rankine(          &
        x,                                         &
        face_nodes, face_center, face_normal,      &
        face_area, face_radius,                    &
        derivative_with_respect_to_first_variable, &
        [-ONE, -depth*lambda_k],                   &
        int_G_term_Rankine, int_nablaG_term_Rankine)
      int_G = int_G + a_k/2 * int_G_term_Rankine
      int_nablaG = int_nablaG + a_k/2 * int_nablaG_term_Rankine

      ! 2.c
      call integral_of_reflected_Rankine(          &
        x,                                         &
        face_nodes, face_center, face_normal,      &
        face_area, face_radius,                    &
        derivative_with_respect_to_first_variable, &
        [-ONE, depth*lambda_k - 4*depth],          &
        int_G_term_Rankine, int_nablaG_term_Rankine)
      int_G = int_G + a_k/2 * int_G_term_Rankine
      int_nablaG = int_nablaG + a_k/2 * int_nablaG_term_Rankine

      ! 2.d
      call integral_of_reflected_Rankine(          &
        x,                                         &
        face_nodes, face_center, face_normal,      &
        face_area, face_radius,                    &
        derivative_with_respect_to_first_variable, &
        [ONE, -depth*lambda_k + 2*depth],          &
        int_G_term_Rankine, int_nablaG_term_Rankine)
      int_G = int_G + a_k/2 * int_G_term_Rankine
      int_nablaG = int_nablaG + a_k/2 * int_nablaG_term_Rankine
    end do
  end subroutine

  ! =====================================================================


  function compute_dispersion_roots(nk, omega2_over_g, depth) result(roots_of_dispersion_relationship)
    integer, intent(in) :: nk
    real(kind=pre), intent(in) :: omega2_over_g, depth
    real(kind=pre), dimension(nk) :: roots_of_dispersion_relationship

    real(kind=8) :: omega, depth_f64
    real(kind=8), dimension(nk) :: roots_of_dispersion_relationship_f64
    depth_f64 = depth
    omega = sqrt(omega2_over_g * 9.81d0)

#ifdef FINGREEN3D_OPTIONAL_DEPENDENCY
    ! Calls FinGreen3D.f90
    call dispersion(roots_of_dispersion_relationship_f64, nk, omega, depth_f64)
#else
    print*, "The library has not been compiled with FinGreen3D.f90 optional dependency"
    error stop
#endif

    roots_of_dispersion_relationship = roots_of_dispersion_relationship_f64
  end function

  subroutine integral_of_wave_part_fingreen3D                  &
    (x,                                                        &
    face_quadrature_points, face_quadrature_weights,           &
    wavenumber, depth, dispersion_roots,                       &
    derivative_with_respect_to_first_variable,                 &
    int_G, int_nablaG                                          &
    )

    real(kind=pre), dimension(3),          intent(in) :: x
    real(kind=pre), dimension(:),          intent(in) :: face_quadrature_weights
    real(kind=pre), dimension(:, :),       intent(in) :: face_quadrature_points
    real(kind=pre),                        intent(in) :: wavenumber, depth
    real(kind=pre), dimension(:),          intent(in) :: dispersion_roots
    logical,                               intent(in) :: derivative_with_respect_to_first_variable

    complex(kind=pre),                     intent(out) :: int_G
    complex(kind=pre), dimension(3),       intent(out) :: int_nablaG

    integer, parameter :: nk = 200
    integer :: q, nb_quad_points
    real(kind=pre) :: omega2_over_g, r, drdx1, drdx2
    real(kind=pre), dimension(3) :: xi_q
    real(kind=pre), dimension(nk) :: roots_of_dispersion_relationship
    complex(kind=pre) :: G_at_point
    complex(kind=pre), dimension(3) :: reduced_G_nablaG, nablaG_at_point

    omega2_over_g  = wavenumber*TANH(wavenumber*depth)

    int_G = czero
    int_nablaG = czero

    nb_quad_points = size(face_quadrature_weights)

    do q = 1, nb_quad_points
      xi_q = face_quadrature_points(q, :)
      r = norm2(x(1:2) - xi_q(1:2))
#ifdef FINGREEN3D_OPTIONAL_DEPENDENCY
      if (.not. derivative_with_respect_to_first_variable) then
        ! For direct method as implemented in HAMS
        call fingreen3d_routine(r, x(3), xi_q(3), omega2_over_g, dispersion_roots, nk, depth, reduced_G_nablaG)
      else
        ! Switched inputs
        call fingreen3d_routine(r, xi_q(3), x(3), omega2_over_g, dispersion_roots, nk, depth, reduced_G_nablaG)
      endif
#else
    print*, "The library has not been compiled with FinGreen3D.f90 optional dependency"
    error stop
#endif

      G_at_point = reduced_G_nablaG(1)

      if (abs(r) > 16*epsilon(r)) then
        if (derivative_with_respect_to_first_variable) then
          drdx1 = (x(1) - xi_q(1))/r
          drdx2 = (x(2) - xi_q(2))/r
        else
          drdx1 = -(x(1) - xi_q(1))/r
          drdx2 = -(x(2) - xi_q(2))/r
        endif
      else
        ! Limit when r->0 is not well defined...
        drdx1 = ZERO
        drdx2 = ZERO
      end if
      nablaG_at_point(1) = drdx1 * reduced_G_nablaG(2)
      nablaG_at_point(2) = drdx2 * reduced_G_nablaG(2)
      nablaG_at_point(3) = reduced_G_nablaG(3)

      int_G = int_G + G_at_point * face_quadrature_weights(q)
      int_nablaG(:) = int_nablaG(:) + nablaG_at_point(:) * face_quadrature_weights(q)
    enddo
  end subroutine

end module Green_Wave
