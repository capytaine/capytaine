! Copyright (C) 2017-2019 Matthieu Ancellin
! See LICENSE file at <https://github.com/mancellin/libDelhommeau>
MODULE GREEN_WAVE

  USE FLOATING_POINT_PRECISION, ONLY: PRE
  USE CONSTANTS
  USE DELHOMMEAU_INTEGRALS
  USE GREEN_RANKINE, ONLY: COMPUTE_ASYMPTOTIC_RANKINE_SOURCE

  IMPLICIT NONE

  ! Dependencies between the functions of this module:
  ! (from top to bottom: "is called by")
  !
  !               (Delhommeau_integrals.f90)
  !                          |
  !              COLLECT_DELHOMMEAU_INTEGRALS       (COMPUTE_ASYMPTOTIC_RANKINE_SOURCE)
  !                        /   \                    /
  ! WAVE_PART_INFINITE_DEPTH   WAVE_PART_FINITE_DEPTH
  !                        \   /
  !                 INTEGRAL_OF_WAVE_PART
  !                          |
  !                    (matrices.f90)
  !                          |
  !                    (python code)

CONTAINS

  ! =====================================================================

  subroutine integral_of_wave_part                               &
      (x,                                                        &
      face_center, face_normal, face_area, face_radius,          &
      face_quadrature_points, face_quadrature_weights,           &
      wavenumber, depth,                                         &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      nexp, ambda, ar,                                           &
      int_G, int_grad_G_sym, int_grad_G_antisym                  &
      )
    ! Integral over a panel of the wave part of the Green function.

    real(kind=pre), dimension(3),          intent(in) :: x
    real(kind=pre), dimension(3),          intent(in) :: face_center, face_normal
    real(kind=pre),                        intent(in) :: face_area, face_radius
    real(kind=pre), dimension(:),          intent(in) :: face_quadrature_weights
    real(kind=pre), dimension(:, :),       intent(in) :: face_quadrature_points
    real(kind=pre),                        intent(in) :: wavenumber, depth
    real(kind=pre), dimension(:),          intent(in) :: tabulated_r_range
    real(kind=pre), dimension(:),          intent(in) :: tabulated_z_range
    real(kind=pre), dimension(:, :, :, :), intent(in) :: tabulated_integrals
    integer,                               intent(in) :: nexp
    real(kind=pre), dimension(nexp),       intent(in) :: ambda, ar

    complex(kind=pre),                     intent(out) :: int_G
    complex(kind=pre), dimension(3),       intent(out) :: int_grad_G_sym, int_grad_G_antisym

    ! Local variables
    real(kind=pre)                  :: r, z
    complex(kind=pre)               :: G_at_point
    complex(kind=pre), dimension(3) :: grad_G_at_point_sym, grad_G_at_point_antisym
    integer                         :: nb_quad_points, Q

    nb_quad_points = size(face_quadrature_weights)

    r = wavenumber * norm2(x(1:2) - face_center(1:2))
    z = wavenumber * (x(3) + face_center(3))

    if ((abs(r) < 1e-10) .and. (abs(z) < 1e-10)) then
      ! Interaction of a panel on the free surface with itself
      int_G = wavenumber * face_area * (                      &
                    (1 - log(wavenumber**2 * face_area / pi)) &
                    + 2 * (euler_gamma - log_2) - 2*pi*ii    &
                  )
      int_grad_G_sym(1:2) = cmplx(zero, zero, kind=pre)
      int_grad_G_sym(3) = wavenumber * (int_G + 4 * sqrt(pi*face_area))
      int_grad_G_antisym(1:3) = cmplx(zero, zero, kind=pre)  ! Irrelevant because we are on the diagonal

    else
      ! Numerical integration
      int_G = zero
      int_grad_G_sym = zero
      int_grad_G_antisym = zero

      do q = 1, nb_quad_points
        if (is_infinity(depth)) then
          call wave_part_infinite_depth                                &
            (x,                                                        &
            face_quadrature_points(q, :),                              &
            wavenumber,                                                &
            tabulated_r_range, tabulated_z_range, tabulated_integrals, &
            g_at_point, grad_g_at_point_sym, grad_g_at_point_antisym   &
            )
        else
          call wave_part_finite_depth                                  &
            (x,                                                        &
            face_quadrature_points(q, :),                              &
            wavenumber,                                                &
            depth,                                                     &
            tabulated_r_range, tabulated_z_range, tabulated_integrals, &
            nexp, ambda, ar,                                           &
            g_at_point, grad_g_at_point_sym, grad_g_at_point_antisym   &
            )
        end if

        int_g = int_g + g_at_point * face_quadrature_weights(q)
        int_grad_g_sym = int_grad_g_sym + grad_g_at_point_sym * face_quadrature_weights(q)
        int_grad_g_antisym = int_grad_g_antisym + grad_g_at_point_antisym * face_quadrature_weights(q)
      end do
    end if

  end subroutine


  SUBROUTINE COLLECT_DELHOMMEAU_INTEGRALS                        &
      (X0I, X0J, wavenumber,                                     &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      FS, VS)

    ! Inputs
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I, X0J
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber

    ! Tabulated data
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_z_range
    REAL(KIND=PRE), DIMENSION(size(tabulated_r_range), size(tabulated_z_range), 2, 2), INTENT(IN) :: tabulated_integrals

    ! Outputs
    COMPLEX(KIND=PRE),                        INTENT(OUT) :: FS  ! the integral
    COMPLEX(KIND=PRE), DIMENSION(3),          INTENT(OUT) :: VS  ! its gradient

    ! Local variables
    REAL(KIND=PRE) :: r, z, drdx1, drdx2, dzdx3
    REAL(KIND=PRE), dimension(2, 2) :: integrals

    r = wavenumber * NORM2(X0I(1:2) - X0J(1:2))
    z = wavenumber * (X0I(3) + X0J(3))

    !=======================================================
    ! Evaluate the elementary integrals depending on z and r
    !=======================================================
    IF ((size(tabulated_z_range) <= 1) .or. (size(tabulated_r_range) <= 1)) THEN
      ! No tabulation, fully recompute the Green function each time.
      integrals = numerical_integration(r, z, 500)
    ELSE
      IF ((abs(z) < abs(tabulated_z_range(size(tabulated_z_range)))) .AND. (r < tabulated_r_range(size(tabulated_r_range)))) THEN
        ! Within the range of tabulated data
        integrals = pick_in_default_tabulation(r, z, tabulated_r_range, tabulated_z_range, tabulated_integrals)
      ELSE
        ! Asymptotic expression for distant panels
        integrals = asymptotic_approximations(MAX(r, 1e-10), z)
      ENDIF
    ENDIF

    !================================================
    ! Add the elementary integrals to build FS and VS
    !================================================

    IF (ABS(r) > 16*EPSILON(r)) THEN
      drdx1 = wavenumber**2 * (X0I(1) - X0J(1))/r
      drdx2 = wavenumber**2 * (X0I(2) - X0J(2))/r
    ELSE
      ! Limit when r->0 is not well defined...
      drdx1 = ZERO
      drdx2 = ZERO
    END IF
    dzdx3 = wavenumber

    FS    = CMPLX(integrals(1, 2), integrals(2, 2), KIND=PRE)
    VS(1) = -drdx1 * CMPLX(integrals(1, 1), integrals(2, 1), KIND=PRE)
    VS(2) = -drdx2 * CMPLX(integrals(1, 1), integrals(2, 1), KIND=PRE)
#ifdef XIE_CORRECTION
    VS(3) = dzdx3 * CMPLX(integrals(1, 2) + ONE/hypot(r, z), integrals(2, 2), KIND=PRE)
#else
    VS(3) = dzdx3 * CMPLX(integrals(1, 2), integrals(2, 2), KIND=PRE)
#endif

    RETURN
  END SUBROUTINE COLLECT_DELHOMMEAU_INTEGRALS

  ! =========================

  SUBROUTINE WAVE_PART_INFINITE_DEPTH &
      (X0I, X0J, wavenumber,          &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      SP, VSP_SYM, VSP_ANTISYM)
    ! Compute the wave part of the Green function in the infinite depth case.
    ! This is mostly the integral computed by the subroutine above.

    ! Inputs
    REAL(KIND=PRE),                           INTENT(IN)  :: wavenumber
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN)  :: X0I   ! Coordinates of the source point
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN)  :: X0J   ! Coordinates of the center of the integration panel

    ! Tabulated data
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_z_range
    REAL(KIND=PRE), DIMENSION(size(tabulated_r_range), size(tabulated_z_range), 2, 2), INTENT(IN) :: tabulated_integrals

    ! Outputs
    COMPLEX(KIND=PRE),               INTENT(OUT) :: SP  ! Integral of the Green function over the panel.
    COMPLEX(KIND=PRE), DIMENSION(3), INTENT(OUT) :: VSP_SYM, VSP_ANTISYM ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    REAL(KIND=PRE), DIMENSION(3) :: XJ_REFLECTION
    COMPLEX(KIND=PRE), DIMENSION(3) :: VSP

    ! The integrals
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      X0I, X0J, wavenumber,                                      &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      SP, VSP(:))
    SP  = 2*wavenumber*SP
    VSP = 2*wavenumber*VSP

#ifndef XIE_CORRECTION
    ! In the original Delhommeau method
    XJ_REFLECTION(1:2) = X0J(1:2)
    XJ_REFLECTION(3) = - X0J(3)
    ! Only one singularity is missing in the derivative
    VSP = VSP - 2*(X0I - XJ_REFLECTION)/(NORM2(X0I-XJ_REFLECTION)**3)
#endif

    VSP_SYM(1:2)     = CMPLX(ZERO, ZERO, KIND=PRE)
    VSP_SYM(3)       = VSP(3)
    VSP_ANTISYM(1:2) = VSP(1:2)
    VSP_ANTISYM(3)   = CMPLX(ZERO, ZERO, KIND=PRE)

    RETURN
  END SUBROUTINE WAVE_PART_INFINITE_DEPTH

  ! ======================

  SUBROUTINE WAVE_PART_FINITE_DEPTH &
      (X0I, X0J, wavenumber, depth, &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      NEXP, AMBDA, AR,              &
      SP, VSP_SYM, VSP_ANTISYM)
    ! Compute the frequency-dependent part of the Green function in the finite depth case.

    ! Inputs
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber, depth
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I  ! Coordinates of the source point
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0J  ! Coordinates of the center of the integration panel

    ! Tabulated data
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_z_range
    REAL(KIND=PRE), DIMENSION(size(tabulated_r_range), size(tabulated_z_range), 2, 2), INTENT(IN) :: tabulated_integrals

    ! Prony decomposition for finite depth
    INTEGER,                                  INTENT(IN) :: NEXP
    REAL(KIND=PRE), DIMENSION(NEXP),          INTENT(IN) :: AMBDA, AR

    ! Outputs
    COMPLEX(KIND=PRE),               INTENT(OUT) :: SP  ! Integral of the Green function over the panel.
    COMPLEX(KIND=PRE), DIMENSION(3), INTENT(OUT) :: VSP_SYM, VSP_ANTISYM ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    INTEGER                              :: KE
    REAL(KIND=PRE)                       :: AMH, AKH, A
    REAL(KIND=PRE)                       :: AQT, R
    REAL(KIND=PRE),    DIMENSION(3)      :: XI, XJ
    REAL(KIND=PRE),    DIMENSION(4)      :: FTS, PSR
    REAL(KIND=PRE),    DIMENSION(3, 4)   :: VTS
    COMPLEX(KIND=PRE), DIMENSION(4)      :: FS
    COMPLEX(KIND=PRE), DIMENSION(3, 4)   :: VS

    !========================================
    ! Part 1: Solve 4 infinite depth problems
    !========================================

    XI(:) = X0I(:)
    XJ(:) = X0J(:)

    ! Distance in xOy plane
    R = NORM2(XI(1:2) - XJ(1:2))

    ! 1.a First infinite depth problem
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      XI(:), XJ(:), wavenumber,                                  &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      FS(1), VS(:, 1))

#ifndef XIE_CORRECTION
    ! In the original Delhommeau method, the integrals are Re[ ∫(J(ζ) - 1/ζ)dθ ]/π + i Re[ ∫(e^ζ)dθ ]
    ! whereas we need Re[ ∫(J(ζ))dθ ]/π + i Re[ ∫(e^ζ)dθ ]
    ! So the term PSR is the difference because  Re[ ∫ 1/ζ dθ ] = - π/sqrt(r² + z²)
    !
    ! Note however, that the derivative part of Delhommeau integrals is the derivative of
    ! Re[ ∫(J(ζ))dθ ]/π + i Re[ ∫(e^ζ)dθ ] so no fix is needed for the derivative.
    PSR(1) = ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))
#endif

    ! 1.b Shift and reflect XI and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) =  X0J(3)
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      XI(:), XJ(:), wavenumber,                                  &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      FS(2), VS(:, 2))
    VS(3, 2) = -VS(3, 2) ! Reflection of the output vector

#ifndef XIE_CORRECTION
    PSR(2) = ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))
#endif

    ! 1.c Shift and reflect XJ and compute another value of the Green function
    XI(3) =  X0I(3)
    XJ(3) = -X0J(3) - 2*depth
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      XI(:), XJ(:), wavenumber,                                  &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      FS(3), VS(:, 3))

#ifndef XIE_CORRECTION
    PSR(3) = ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))
#endif

    ! 1.d Shift and reflect both XI and XJ and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) = -X0J(3) - 2*depth
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      XI(:), XJ(:), wavenumber,                                  &
      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      FS(4), VS(:, 4))
    VS(3, 4) = -VS(3, 4) ! Reflection of the output vector

#ifndef XIE_CORRECTION
    PSR(4) = ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))
#endif

    ! Add up the results of the four problems
#ifdef XIE_CORRECTION
    SP               = SUM(FS(1:4))
#else
    SP               = SUM(FS(1:4)) - SUM(PSR(1:4))
#endif
    VSP_SYM(1:2)     = CMPLX(ZERO, ZERO, KIND=PRE)
    VSP_ANTISYM(1:2) = VS(1:2, 1) + VS(1:2, 2) + VS(1:2, 3) + VS(1:2, 4)
    VSP_SYM(3)       = VS(3, 1) + VS(3, 4)
    VSP_ANTISYM(3)   = VS(3, 2) + VS(3, 3)

    ! Multiply by some coefficients
    AMH  = wavenumber*depth
    AKH  = AMH*TANH(AMH)
    A    = (AMH+AKH)**2/(2*depth*(AMH**2-AKH**2+AKH))

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
