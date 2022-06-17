! Copyright (C) 2017-2019 Matthieu Ancellin
! See LICENSE file at <https://github.com/mancellin/capytaine>
MODULE GREEN_WAVE

  USE CONSTANTS
  USE DELHOMMEAU_INTEGRALS
  USE GREEN_RANKINE, ONLY: COMPUTE_ASYMPTOTIC_RANKINE_SOURCE

  IMPLICIT NONE

  ! Dependencies between the functions of this module:
  ! (from top to bottom: "is called by")
  !
  !               (MODULE DELHOMMEAU_INTEGRALS)
  !                          |
  !              COLLECT_DELHOMMEAU_INTEGRALS       (COMPUTE_ASYMPTOTIC_RANKINE_SOURCE)
  !                        /   \                    /
  ! WAVE_PART_INFINITE_DEPTH   WAVE_PART_FINITE_DEPTH
  !                        \   /
  !                  (build_matrices.f90)
  !                          |
  !                    (python code)

CONTAINS

  ! =====================================================================

  SUBROUTINE COLLECT_DELHOMMEAU_INTEGRALS                        &
      (XI, XJ, wavenumber,                                       &
      tabulated_r_range, tabulated_Z_range, tabulated_integrals, &
      FS, VS)

    ! Inputs
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: XI, XJ
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber

    ! Tabulated data
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_Z_range
    REAL(KIND=PRE), DIMENSION(size(tabulated_r_range), size(tabulated_Z_range), 2, 2), INTENT(IN) :: tabulated_integrals

    ! Outputs
    COMPLEX(KIND=PRE),                        INTENT(OUT) :: FS  ! the integral
    COMPLEX(KIND=PRE), DIMENSION(3),          INTENT(OUT) :: VS  ! its gradient

    ! Local variables
    REAL(KIND=PRE) :: r, z, r1
    REAL(KIND=PRE), dimension(2, 2) :: integrals
    REAL(KIND=PRE) :: D1, D2, Z1, Z2

    r = wavenumber * NORM2(XI(1:2) - XJ(1:2))
    z = wavenumber * (XI(3) + XJ(3))
    r1 = hypot(r, z)

    IF (z > -1e-8) THEN
      PRINT*, "Error: Impossible to compute the wave part of the Green function due to panels on the free surface (z=0) or above."
      ERROR STOP
    ENDIF

    !=====================================================================================
    ! Evaluate the elementary integrals PDnX and PDnZ depending on dimless_Z and dimless_r
    !=====================================================================================
    IF ((MINVAL(tabulated_Z_range) < z) .AND. (z < MAXVAL(tabulated_Z_range))) THEN
      IF ((MINVAL(tabulated_r_range) <= r) .AND. (r < MAXVAL(tabulated_r_range))) THEN
        ! Within the range of tabulated data
        ! Note that MINVAL(tabulated_r_range) == 0, so one of the conditions is not actually useful.

        integrals = pick_in_default_tabulation(r, z, tabulated_r_range, tabulated_Z_range, tabulated_integrals)

      ELSE  ! MAXVAL(tabulated_r_range) < dimless_r
        ! Asymptotic expression for (horizontally) distant panels
        integrals = asymptotic_approximations(r, z)
      ENDIF

      !================================================
      ! Add the elementary integrals to build FS and VS
      !================================================

      D1 = integrals(1, 1)
      D2 = integrals(2, 1)
      Z1 = integrals(1, 2)
      Z2 = integrals(2, 2)

#ifdef XIE_CORRECTION
      FS    = CMPLX(Z1/PI + ONE/R1, Z2, KIND=PRE)
      VS(3) = CMPLX(Z1/PI + ONE/R1, Z2, KIND=PRE)
#else
      FS    = CMPLX(Z1/PI, Z2, KIND=PRE)
      VS(3) = CMPLX(Z1/PI, Z2, KIND=PRE)
#endif
      VS(1) = wavenumber*(XJ(1) - XI(1))/r * CMPLX(D1/PI, D2, KIND=PRE)
      VS(2) = wavenumber*(XJ(2) - XI(2))/r * CMPLX(D1/PI, D2, KIND=PRE)

      IF (r < REAL(1e-5, KIND=PRE)) THEN
        ! Limit case r ~ 0 ?
        VS(1:2) = CMPLX(0.0, 0.0, KIND=PRE)
      END IF

    ELSE  ! dimless_Z < MINVAL(tabulated_Z_range) or MAXVAL(tabulated_Z_range) < dimless_Z
      FS      = CMPLX(z/r1**3, 0.0, KIND=PRE)
      VS(1:3) = CMPLX(0.0, 0.0, KIND=PRE)
    ENDIF

    RETURN
  END SUBROUTINE COLLECT_DELHOMMEAU_INTEGRALS

  ! =========================

  SUBROUTINE WAVE_PART_INFINITE_DEPTH &
      (wavenumber, X0I, X0J,          &
      X_AXIS, Z_AXIS, TABULATION,     &
      SP, VSP)
    ! Compute the wave part of the Green function in the infinite depth case.
    ! This is mostly the integral computed by the subroutine above.

    ! Inputs
    REAL(KIND=PRE),                           INTENT(IN)  :: wavenumber
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN)  :: X0I   ! Coordinates of the source point
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN)  :: X0J   ! Coordinates of the center of the integration panel

    ! Tabulated data
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: X_AXIS
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: Z_AXIS
    REAL(KIND=PRE), DIMENSION(size(X_AXIS), size(Z_AXIS), 2, 2), INTENT(IN) :: TABULATION

    ! Outputs
    COMPLEX(KIND=PRE),               INTENT(OUT) :: SP  ! Integral of the Green function over the panel.
    COMPLEX(KIND=PRE), DIMENSION(3), INTENT(OUT) :: VSP ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    REAL(KIND=PRE), DIMENSION(3) :: XJ_REFLECTION

    ! The integrals
    CALL COLLECT_DELHOMMEAU_INTEGRALS(X0I, X0J, wavenumber, X_AXIS, Z_AXIS, TABULATION, SP, VSP(:))
    SP  = 2*wavenumber*SP
    VSP = 2*wavenumber**2*VSP

    ! Only one singularity is missing in the derivative
    XJ_REFLECTION(1:2) = X0J(1:2)
    XJ_REFLECTION(3) = - X0J(3)
    VSP = VSP - 2*(X0I - XJ_REFLECTION)/(NORM2(X0I-XJ_REFLECTION)**3)

    RETURN
  END SUBROUTINE WAVE_PART_INFINITE_DEPTH

  ! ======================

  SUBROUTINE WAVE_PART_FINITE_DEPTH &
      (wavenumber, X0I, X0J, depth, &
      X_AXIS, Z_AXIS, TABULATION,                  &
      NEXP, AMBDA, AR,              &
      SP, VSP_SYM, VSP_ANTISYM)
    ! Compute the frequency-dependent part of the Green function in the finite depth case.

    ! Inputs
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber, depth
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I  ! Coordinates of the source point
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0J  ! Coordinates of the center of the integration panel

    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: X_AXIS
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: Z_AXIS
    REAL(KIND=PRE), DIMENSION(size(X_AXIS), size(Z_AXIS), 2, 2), INTENT(IN) :: TABULATION

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
    CALL COLLECT_DELHOMMEAU_INTEGRALS(XI(:), XJ(:), wavenumber, X_AXIS, Z_AXIS, TABULATION, FS(1), VS(:, 1))

    PSR(1) = ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))

    ! 1.b Shift and reflect XI and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) =  X0J(3)
    CALL COLLECT_DELHOMMEAU_INTEGRALS(XI(:), XJ(:), wavenumber, X_AXIS, Z_AXIS, TABULATION, FS(2), VS(:, 2))
    VS(3, 2) = -VS(3, 2) ! Reflection of the output vector

    PSR(2) = ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))

    ! 1.c Shift and reflect XJ and compute another value of the Green function
    XI(3) =  X0I(3)
    XJ(3) = -X0J(3) - 2*depth
    CALL COLLECT_DELHOMMEAU_INTEGRALS(XI(:), XJ(:), wavenumber, X_AXIS, Z_AXIS, TABULATION, FS(3), VS(:, 3))

    PSR(3) = ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))

    ! 1.d Shift and reflect both XI and XJ and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) = -X0J(3) - 2*depth
    CALL COLLECT_DELHOMMEAU_INTEGRALS(XI(:), XJ(:), wavenumber, X_AXIS, Z_AXIS, TABULATION, FS(4), VS(:, 4))
    VS(3, 4) = -VS(3, 4) ! Reflection of the output vector

    PSR(4) = ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))

    ! Add up the results of the four problems
    SP               = SUM(FS(1:4)) - SUM(PSR(1:4))
    VSP_SYM(1:3)     = VS(1:3, 1) + VS(1:3, 4)
    VSP_ANTISYM(1:3) = VS(1:3, 2) + VS(1:3, 3)

    ! Multiply by some coefficients
    AMH  = wavenumber*depth
    AKH  = AMH*TANH(AMH)
    A    = (AMH+AKH)**2/(2*depth*(AMH**2-AKH**2+AKH))

    SP          = A*SP
    VSP_ANTISYM = A*wavenumber*VSP_ANTISYM
    VSP_SYM     = A*wavenumber*VSP_SYM

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
      VSP_ANTISYM(1:3) = VSP_ANTISYM(1:3) + AQT*(VTS(1:3, 1) + VTS(1:3, 4))
      VSP_SYM(1:3)     = VSP_SYM(1:3)     + AQT*(VTS(1:3, 2) + VTS(1:3, 3))

    END DO

    RETURN
  END SUBROUTINE

END MODULE GREEN_WAVE
