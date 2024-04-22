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
  !              WAVE_PART_INFINITE_DEPTH           (COMPUTE_ASYMPTOTIC_RANKINE_SOURCE)
  !                        /   \                    /
  !                       |    WAVE_PART_FINITE_DEPTH
  !                       \   /
  !                    (matrices.f90)
  !                          |
  !                    (python code)

CONTAINS

  ! =====================================================================

  SUBROUTINE WAVE_PART_INFINITE_DEPTH                        &
      ! Returns (G^-, nabla G^+) if gf_singularities == HIGH_FREQ
      ! and (G^+, nabla G^+) if gf_singularities == LOW_FREQ
      (X0I, X0J, wavenumber,                                     &
      tabulation_grid_shape, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      gf_singularities,                                         &
      G, nablaG)

    ! Inputs
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I, X0J
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber
    integer,                                  intent(in) :: gf_singularities

    ! Tabulated data
    INTEGER,                                  INTENT(IN) :: tabulation_grid_shape
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(:),             INTENT(IN) :: tabulated_z_range
    REAL(KIND=PRE), DIMENSION(size(tabulated_r_range), size(tabulated_z_range), 2, 2), INTENT(IN) :: tabulated_integrals

    ! Outputs
    COMPLEX(KIND=PRE),                        INTENT(OUT) :: G
    COMPLEX(KIND=PRE), DIMENSION(3),          INTENT(OUT) :: nablaG

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
      integrals = numerical_integration(r, z, 500)
      if (gf_singularities == HIGH_FREQ) then
        ! numerical_integration always computes the low_freq version,
        ! so need a fix to get the high_freq
        integrals(1, 2) = integrals(1, 2) + 2/r1
      endif
    ELSE
      IF ((abs(z) < abs(tabulated_z_range(size(tabulated_z_range)))) &
          .AND. (r < tabulated_r_range(size(tabulated_r_range)))) THEN
        ! Within the range of tabulated data
        integrals = pick_in_default_tabulation( &
            r, z, tabulation_grid_shape, tabulated_r_range, tabulated_z_range, tabulated_integrals &
        )
        if (gf_singularities == HIGH_FREQ) then
          ! numerical_integration always computes the low_freq version,
          ! so need a fix to get the high_freq
          integrals(1, 2) = integrals(1, 2) + 2/r1
        endif
      ELSE
        ! Delhommeau's asymptotic expression of Green function for distant panels
        integrals = asymptotic_approximations(MAX(r, 1e-10), z)
        if (gf_singularities == LOW_FREQ) then
          ! numerical_integration always computes the high_freq version,
          ! so need a fix to get the low_freq
          integrals(1, 2) = integrals(1, 2) - 2/r1
        endif
      ENDIF
    ENDIF

    !================================================
    ! Add the elementary integrals to build FS and VS
    !================================================

    G    = CMPLX(integrals(1, 2), integrals(2, 2), KIND=PRE)
    nablaG(1) = -drdx1 * CMPLX(integrals(1, 1), integrals(2, 1), KIND=PRE)
    nablaG(2) = -drdx2 * CMPLX(integrals(1, 1), integrals(2, 1), KIND=PRE)
    nablaG(3) = dzdx3 * CMPLX(integrals(1, 2), integrals(2, 2), KIND=PRE)
    if (gf_singularities == LOW_FREQ) then
      ! integrals(:, 2) are G^+, but the formula is that dG^+/dz = G^-
      ! Here is the correction
      nablaG(3) = nablaG(3) + 2*dzdx3/r1
    elseif (gf_singularities == HIGH_FREQ) then
      ! we have nabla G^+, but we actually need nabla G^-
      nablaG(1) = nablaG(1) - 2*drdx1*r/r1**3
      nablaG(2) = nablaG(2) - 2*drdx2*r/r1**3
      nablaG(3) = nablaG(3) - 2*dzdx3*z/r1**3
    endif

    G = wavenumber * G
    nablaG = wavenumber * nablaG

  END SUBROUTINE WAVE_PART_INFINITE_DEPTH

  ! ======================

  SUBROUTINE WAVE_PART_FINITE_DEPTH &
      (X0I, X0J, wavenumber, depth, &
      tabulation_grid_shape, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      NEXP, AMBDA, AR,              &
      SP, VSP_SYM, VSP_ANTISYM)
    ! Compute the frequency-dependent part of the Green function in the finite depth case.

    ! Inputs
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber, depth
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I  ! Coordinates of the source point
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0J  ! Coordinates of the center of the integration panel

    ! Tabulated data
    INTEGER,                                  INTENT(IN) :: tabulation_grid_shape
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
    REAL(KIND=PRE),    DIMENSION(4)      :: FTS
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
    CALL WAVE_PART_INFINITE_DEPTH(                               &
      XI(:), XJ(:), wavenumber,                                  &
      tabulation_grid_shape, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      LOW_FREQ, &
      FS(1), VS(:, 1))

    ! 1.b Shift and reflect XI and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) =  X0J(3)
    CALL WAVE_PART_INFINITE_DEPTH(                               &
      XI(:), XJ(:), wavenumber,                                  &
      tabulation_grid_shape, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      LOW_FREQ, &
      FS(2), VS(:, 2))
    VS(3, 2) = -VS(3, 2) ! Reflection of the output vector

    ! 1.c Shift and reflect XJ and compute another value of the Green function
    XI(3) =  X0I(3)
    XJ(3) = -X0J(3) - 2*depth
    CALL WAVE_PART_INFINITE_DEPTH(                               &
      XI(:), XJ(:), wavenumber,                                  &
      tabulation_grid_shape, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      LOW_FREQ, &
      FS(3), VS(:, 3))

    ! 1.d Shift and reflect both XI and XJ and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) = -X0J(3) - 2*depth
    CALL WAVE_PART_INFINITE_DEPTH(                               &
      XI(:), XJ(:), wavenumber,                                  &
      tabulation_grid_shape, tabulated_r_range, tabulated_z_range, tabulated_integrals, &
      LOW_FREQ, &
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
