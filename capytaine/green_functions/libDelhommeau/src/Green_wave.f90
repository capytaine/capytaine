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
  !                    (matrices.f90)
  !                          |
  !                    (python code)

CONTAINS

  ! =====================================================================

  SUBROUTINE COLLECT_DELHOMMEAU_INTEGRALS                        &
      ! Returns (G^-, nabla G^+) if legacy_delhommeau and (G^+, nabla G^+) otherwise
      (X0I, X0J, wavenumber,                                     &
      tabulation_method, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      legacy_delhommeau,                                         &
      FS, VS)

    ! Inputs
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I, X0J
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber
    logical,                                  intent(in) :: legacy_delhommeau

    ! Tabulated data
    INTEGER,                                  INTENT(IN) :: tabulation_method
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_z_range
    REAL(KIND=PRE), DIMENSION(size(tabulated_r_range), size(tabulated_z_range), 2, 2), INTENT(IN) :: tabulated_integrals

    ! Outputs
    COMPLEX(KIND=PRE),                        INTENT(OUT) :: FS  ! the integral
    COMPLEX(KIND=PRE), DIMENSION(3),          INTENT(OUT) :: VS  ! its gradient

    ! Local variables
    REAL(KIND=PRE) :: r, z, r1, drdx1, drdx2, dzdx3
    REAL(KIND=PRE), dimension(2, 2) :: integrals

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

    !=======================================================
    ! Evaluate the elementary integrals depending on z and r
    !=======================================================
    IF ((size(tabulated_z_range) <= 1) .or. (size(tabulated_r_range) <= 1)) THEN
      ! No tabulation, fully recompute the Green function each time.
      integrals = numerical_integration(r, z, 500, legacy_delhommeau)
    ELSE
      IF ((abs(z) < abs(tabulated_z_range(size(tabulated_z_range)))) .AND. (r < tabulated_r_range(size(tabulated_r_range)))) THEN
        ! Within the range of tabulated data
        integrals = pick_in_default_tabulation(r, z, tabulation_method, tabulated_r_range, tabulated_z_range, tabulated_integrals)
      ELSE
        ! Asymptotic expression of legacy's Delhommeau Green function for distant panels
        integrals = asymptotic_approximations(MAX(r, 1e-10), z)
        if (.not. legacy_delhommeau) then
          ! asymptotic_approximations returns always G^-
          ! Here is the correction to retrieve G^+
          integrals(1, 2) = integrals(1, 2) - 2/r1
        endif
      ENDIF
    ENDIF

    !================================================
    ! Add the elementary integrals to build FS and VS
    !================================================

    FS    = CMPLX(integrals(1, 2), integrals(2, 2), KIND=PRE)
    VS(1) = -drdx1 * CMPLX(integrals(1, 1), integrals(2, 1), KIND=PRE)
    VS(2) = -drdx2 * CMPLX(integrals(1, 1), integrals(2, 1), KIND=PRE)
    VS(3) = dzdx3 * CMPLX(integrals(1, 2), integrals(2, 2), KIND=PRE)
    if (.not. legacy_delhommeau) then
      ! integrals(:, 2) are G^+, but the formula is that dG^+/dz = G^-
      ! Here is the correction
      VS(3) = VS(3) + 2*dzdx3/r1
    endif

    RETURN
  END SUBROUTINE COLLECT_DELHOMMEAU_INTEGRALS

  ! =========================

  SUBROUTINE WAVE_PART_INFINITE_DEPTH &
      (X0I, X0J, wavenumber,          &
      tabulation_method, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      legacy_delhommeau, &
      SP, VSP_SYM, VSP_ANTISYM)
    ! Compute the wave part of the Green function in the infinite depth case.
    ! This is mostly the integral computed by the subroutine above.

    ! Inputs
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I   ! Coordinates of the source point
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0J   ! Coordinates of the center of the integration panel
    logical,                                  intent(in) :: legacy_delhommeau

    ! Tabulated data
    INTEGER,                                  INTENT(IN) :: tabulation_method
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
      tabulation_method, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      legacy_delhommeau, &
      SP, VSP(:))
    SP  = wavenumber*SP
    VSP = wavenumber*VSP

    if (legacy_delhommeau) then
      ! COLLECT_DELHOMMEAU_INTEGRALS returned nabla G^+, but we actually needed nabla G^-
      ! Here is the correction:
      XJ_REFLECTION(1:2) = X0J(1:2)
      XJ_REFLECTION(3) = - X0J(3)
      ! Only one singularity is missing in the derivative
      VSP = VSP - 2*(X0I - XJ_REFLECTION)/(NORM2(X0I-XJ_REFLECTION)**3)
    endif

    VSP_SYM(1:2)     = CMPLX(ZERO, ZERO, KIND=PRE)
    VSP_SYM(3)       = VSP(3)
    VSP_ANTISYM(1:2) = VSP(1:2)
    VSP_ANTISYM(3)   = CMPLX(ZERO, ZERO, KIND=PRE)

    RETURN
  END SUBROUTINE WAVE_PART_INFINITE_DEPTH

  ! ======================

  SUBROUTINE WAVE_PART_FINITE_DEPTH &
      (X0I, X0J, wavenumber, depth, &
      tabulation_method, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      legacy_delhommeau, &
      NEXP, AMBDA, AR,              &
      SP, VSP_SYM, VSP_ANTISYM)
    ! Compute the frequency-dependent part of the Green function in the finite depth case.

    ! Inputs
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber, depth
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I  ! Coordinates of the source point
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0J  ! Coordinates of the center of the integration panel
    logical,                                  intent(in) :: legacy_delhommeau

    ! Tabulated data
    INTEGER,                                  INTENT(IN) :: tabulation_method
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
      tabulation_method, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      legacy_delhommeau, &
      FS(1), VS(:, 1))

    if (legacy_delhommeau) then
      ! In the original Delhommeau method, the integrals are Re[ ∫(J(ζ) - 1/ζ)dθ ]/π + i Re[ ∫(e^ζ)dθ ]
      ! whereas we need Re[ ∫(J(ζ))dθ ]/π + i Re[ ∫(e^ζ)dθ ]
      ! So the term PSR is the difference because  Re[ ∫ 1/ζ dθ ] = - π/sqrt(r² + z²)
      !
      ! Note however, that the derivative part of Delhommeau integrals is the derivative of
      ! Re[ ∫(J(ζ))dθ ]/π + i Re[ ∫(e^ζ)dθ ] so no fix is needed for the derivative.
      PSR(1) = 2*ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))
    endif

    ! 1.b Shift and reflect XI and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) =  X0J(3)
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      XI(:), XJ(:), wavenumber,                                  &
      tabulation_method, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      legacy_delhommeau, &
      FS(2), VS(:, 2))
    VS(3, 2) = -VS(3, 2) ! Reflection of the output vector

    if (legacy_delhommeau) then
      PSR(2) = 2*ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))
    endif

    ! 1.c Shift and reflect XJ and compute another value of the Green function
    XI(3) =  X0I(3)
    XJ(3) = -X0J(3) - 2*depth
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      XI(:), XJ(:), wavenumber,                                  &
      tabulation_method, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      legacy_delhommeau, &
      FS(3), VS(:, 3))

    if (legacy_delhommeau) then
      PSR(3) = 2*ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))
    endif

    ! 1.d Shift and reflect both XI and XJ and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) = -X0J(3) - 2*depth
    CALL COLLECT_DELHOMMEAU_INTEGRALS(                           &
      XI(:), XJ(:), wavenumber,                                  &
      tabulation_method, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      legacy_delhommeau, &
      FS(4), VS(:, 4))
    VS(3, 4) = -VS(3, 4) ! Reflection of the output vector

    if (legacy_delhommeau) then
      PSR(4) = 2*ONE/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))
    endif

    ! Add up the results of the four problems

    SP               = SUM(FS(1:4))
    if (legacy_delhommeau) then
      SP               = SP - SUM(PSR(1:4))
    endif
    VSP_SYM(1:2)     = CMPLX(ZERO, ZERO, KIND=PRE)
    VSP_ANTISYM(1:2) = VS(1:2, 1) + VS(1:2, 2) + VS(1:2, 3) + VS(1:2, 4)
    VSP_SYM(3)       = VS(3, 1) + VS(3, 4)
    VSP_ANTISYM(3)   = VS(3, 2) + VS(3, 3)

    ! Multiply by some coefficients
    AMH  = wavenumber*depth
    AKH  = AMH*TANH(AMH)
    A    = (AMH+AKH)**2/(4*depth*(AMH**2-AKH**2+AKH))

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
