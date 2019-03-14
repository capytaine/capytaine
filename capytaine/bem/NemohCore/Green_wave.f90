! Copyright (C) 2017-2019 Matthieu Ancellin
! See LICENSE file at <https://github.com/mancellin/capytaine>
MODULE GREEN_WAVE

  USE CONSTANTS
  USE INITIALIZE_GREEN_WAVE
  USE GREEN_RANKINE, ONLY: COMPUTE_ASYMPTOTIC_RANKINE_SOURCE

  IMPLICIT NONE

  ! Dependancies between the functions of this module:
  ! (from top to bottom: "is called by")
  !
  !            LAGRANGE_POLYNOMIAL_INTERPOLATION
  !                          |
  !              COMPUTE_INTEGRALS_WRT_THETA       (COMPUTE_ASYMPTOTIC_RANKINE_SOURCE)
  !                        /   \                    /
  ! WAVE_PART_INFINITE_DEPTH   WAVE_PART_FINITE_DEPTH
  !                        \   /
  !              BUILD_MATRICES_WAVE_SOURCE
  !                          |
  !                    (python code)

CONTAINS

  ! =====================================================================

  SUBROUTINE LAGRANGE_POLYNOMIAL_INTERPOLATION &
    (AKR, AKZ,                                 &
     XR, XZ, APD,                              &
     PD1X, PD2X, PD1Z, PD2Z)
   ! Helper function used in the following subroutine to interpolate between the tabulated integrals.

    ! Inputs
    REAL(KIND=PRE),                        INTENT(IN) :: AKR, AKZ
    REAL(KIND=PRE), DIMENSION(3),          INTENT(IN) :: XR
    REAL(KIND=PRE), DIMENSION(3),          INTENT(IN) :: XZ
    REAL(KIND=PRE), DIMENSION(3, 3, 2, 2), INTENT(IN) :: APD

    ! Output
    REAL(KIND=PRE), INTENT(OUT) :: PD1X, PD2X, PD1Z, PD2Z

    ! Local variable
    REAL(KIND=PRE), DIMENSION(3) :: XL, ZL

    XL(1) = PL2(XR(2), XR(3), XR(1), AKR)
    XL(2) = PL2(XR(3), XR(1), XR(2), AKR)
    XL(3) = PL2(XR(1), XR(2), XR(3), AKR)
    ZL(1) = PL2(XZ(2), XZ(3), XZ(1), AKZ)
    ZL(2) = PL2(XZ(3), XZ(1), XZ(2), AKZ)
    ZL(3) = PL2(XZ(1), XZ(2), XZ(3), AKZ)

    PD1Z = DOT_PRODUCT(XL, MATMUL(APD(:, :, 1, 2), ZL))
    PD2Z = DOT_PRODUCT(XL, MATMUL(APD(:, :, 2, 2), ZL))
    PD1X = DOT_PRODUCT(XL, MATMUL(APD(:, :, 1, 1), ZL))
    PD2X = DOT_PRODUCT(XL, MATMUL(APD(:, :, 2, 1), ZL))

  CONTAINS

    REAL(KIND=PRE) FUNCTION PL2(U1, U2, U3, XU)
      REAL(KIND=PRE) :: U1, U2, U3, XU
      PL2 = ((XU-U1)*(XU-U2))/((U3-U1)*(U3-U2))
      RETURN
    END FUNCTION

  END SUBROUTINE LAGRANGE_POLYNOMIAL_INTERPOLATION

  ! =====================================================================

  SUBROUTINE COMPUTE_INTEGRALS_WRT_THETA                         &
      (XI, XJ, wavenumber,                                       &
      tabulated_r_range, tabulated_Z_range, tabulated_integrals, &
      FS, VS)
    ! Compute the expression FS = Re[ ∫(J(ζ) - 1/ζ)dθ ] + i Re[ ∫(e^ζ)dθ ] and the integrals appearing in its gradient.
    ! For this, this function uses tabulated values of the following integrals:
    ! PD1X = Re[ ∫(-i cosθ)(J(ζ) - 1/ζ)dθ ]
    ! PD2X = Re[ ∫(-i cosθ)(e^ζ)dθ ]
    ! PD1Z = Re[ ∫(J(ζ) - 1/ζ)dθ ]
    ! PD2Z = Re[ ∫(e^ζ)dθ ]
    ! (See theory manual for the definition of the symbols above.)

    ! Inputs
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: XI, XJ
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber

    ! Tabulated data
    REAL(KIND=PRE), DIMENSION(328),           INTENT(IN) :: tabulated_r_range
    REAL(KIND=PRE), DIMENSION(46),            INTENT(IN) :: tabulated_Z_range
    REAL(KIND=PRE), DIMENSION(328, 46, 2, 2), INTENT(IN) :: tabulated_integrals

    ! Outputs
    COMPLEX(KIND=PRE),                        INTENT(OUT) :: FS  ! the integral
    COMPLEX(KIND=PRE), DIMENSION(3),          INTENT(OUT) :: VS  ! its gradient

    ! Local variables
    INTEGER        :: KI, KJ
    REAL(KIND=PRE) :: r, dimless_r, Z, dimless_Z, R1
    REAL(KIND=PRE) :: SIK, CSK, SQ, EPZ
    REAL(KIND=PRE) :: PD1X, PD2X, PD1Z, PD2Z

    r = NORM2(XI(1:2) - XJ(1:2))
    dimless_r = wavenumber*r

    Z = XI(3) + XJ(3)
    dimless_Z = wavenumber*Z

    R1 = SQRT(r**2 + Z**2)

    IF ((R1 < 1e-5) .OR. (wavenumber == 0)) THEN
      PRINT*, "Error: Impossible to compute the wave part of the Green function (division by zero)."
      ERROR STOP
    ENDIF

    !=====================================================================================
    ! Evaluate the elementary integrals PDnX and PDnZ depending on dimless_Z and dimless_r
    !=====================================================================================
    IF ((MINVAL(tabulated_Z_range) < dimless_Z) .AND. (dimless_Z < MAXVAL(tabulated_Z_range))) THEN
      IF ((MINVAL(tabulated_r_range) <= dimless_r) .AND. (dimless_r < MAXVAL(tabulated_r_range))) THEN
        ! Within the range of tabulated data
        ! Note that MINVAL(tabulated_r_range) == 0, so one of the conditions is not actually useful.

        IF (dimless_r < 1) THEN
          KI = INT(5*(LOG10(dimless_r+1e-20)+6)+1)
        ELSE
          KI = INT(3*dimless_r+28)
        ENDIF
        KI = MAX(MIN(KI, 327), 2)

        IF (dimless_Z < -1e-2) THEN
          KJ = INT(8*(LOG10(-dimless_Z)+4.5))
        ELSE
          KJ = INT(5*(LOG10(-dimless_Z)+6))
        ENDIF
        KJ = MAX(MIN(KJ, 45), 2)

        CALL LAGRANGE_POLYNOMIAL_INTERPOLATION                           &
             (dimless_r, dimless_Z,                                      &
             tabulated_r_range(KI-1:KI+1), tabulated_Z_range(KJ-1:KJ+1), &
             tabulated_integrals(KI-1:KI+1, KJ-1:KJ+1, :, :),            &
             PD1X, PD2X, PD1Z, PD2Z)

      ELSE  ! MAXVAL(tabulated_r_range) < dimless_r
        ! Asymptotic expression for (horizontally) distant panels

        EPZ  = EXP(dimless_Z)
        SQ   = SQRT(2*PI/dimless_r)
        CSK  = COS(dimless_r-PI/4)
        SIK  = SIN(dimless_r-PI/4)

        PD1X = PI*EPZ*SQ*(CSK - 0.5*SIK/dimless_r) - PI*r/(wavenumber**2*R1**3)
        PD2X =    EPZ*SQ*(SIK + 0.5*CSK/dimless_r)
        PD1Z = -EPZ*SQ*SIK*PI + PI*dimless_Z/(wavenumber*R1)**3
        PD2Z =  EPZ*SQ*CSK
      ENDIF

      !================================================
      ! Add the elementary integrals to build FS and VS
      !================================================

      FS    = CMPLX(PD1Z, PD2Z, KIND=PRE)
      VS(1) = (XI(1) - XJ(1))/r * CMPLX(-PD1X, -PD2X, KIND=PRE)
      VS(2) = (XI(2) - XJ(2))/r * CMPLX(-PD1X, -PD2X, KIND=PRE)
      VS(3) = CMPLX(PD1Z, PD2Z, KIND=PRE)

      IF (r < REAL(1e-5, KIND=PRE)) THEN
        ! I'm not totally sure why this is here, but it seems to be necessary.
        VS(1:2) = CMPLX(0.0, 0.0, KIND=PRE)
      END IF

    ELSE  ! dimless_Z < MINVAL(tabulated_Z_range) or MAXVAL(tabulated_Z_range) < dimless_Z
      FS      = CMPLX(PI*dimless_Z/(wavenumber*R1)**3, 0.0, KIND=PRE)
      VS(1:3) = CMPLX(0.0, 0.0, KIND=PRE)
    ENDIF

    RETURN
  END SUBROUTINE COMPUTE_INTEGRALS_WRT_THETA

  ! =========================

  SUBROUTINE WAVE_PART_INFINITE_DEPTH &
      (wavenumber, X0I, X0J,          &
      XR, XZ, APD,                    &
      SP, VSP)
    ! Compute the wave part of the Green function in the infinite depth case.
    ! This is mostly the integral computed by the subroutine above.

    ! Inputs
    REAL(KIND=PRE),                           INTENT(IN)  :: wavenumber
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN)  :: X0I   ! Coordinates of the source point
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN)  :: X0J   ! Coordinates of the center of the integration panel

    ! Tabulated data
    REAL(KIND=PRE), DIMENSION(328),           INTENT(IN) :: XR
    REAL(KIND=PRE), DIMENSION(46),            INTENT(IN) :: XZ
    REAL(KIND=PRE), DIMENSION(328, 46, 2, 2), INTENT(IN) :: APD

    ! Outputs
    COMPLEX(KIND=PRE),               INTENT(OUT) :: SP  ! Integral of the Green function over the panel.
    COMPLEX(KIND=PRE), DIMENSION(3), INTENT(OUT) :: VSP ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    REAL(KIND=PRE), DIMENSION(3) :: XJ_REFLECTION

    ! The integrals
    CALL COMPUTE_INTEGRALS_WRT_THETA(X0I, X0J, wavenumber, XR, XZ, APD, SP, VSP(:))
    SP  = 2*CMPLX(wavenumber/PI*REAL(SP),     wavenumber*AIMAG(SP),   KIND=PRE)
    VSP = 2*CMPLX(wavenumber**2/PI*REAL(VSP), wavenumber**2*AIMAG(VSP), KIND=PRE)

    ! Only one singularity is missing in the derivative
    XJ_REFLECTION(1:2) = X0J(1:2)
    XJ_REFLECTION(3) = - X0J(3)
    VSP = VSP - 2*(X0I - XJ_REFLECTION)/(NORM2(X0I-XJ_REFLECTION)**3)

    RETURN
  END SUBROUTINE WAVE_PART_INFINITE_DEPTH

  ! ======================

  SUBROUTINE WAVE_PART_FINITE_DEPTH &
      (wavenumber, X0I, X0J, depth, &
      XR, XZ, APD,                  &
      NEXP, AMBDA, AR,              &
      SP, VSP_SYM, VSP_ANTISYM)
    ! Compute the frequency-dependent part of the Green function in the finite depth case.

    ! Inputs
    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber, depth
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0I  ! Coordinates of the source point
    REAL(KIND=PRE), DIMENSION(3),             INTENT(IN) :: X0J  ! Coordinates of the center of the integration panel

    REAL(KIND=PRE), DIMENSION(328),           INTENT(IN) :: XR
    REAL(KIND=PRE), DIMENSION(46),            INTENT(IN) :: XZ
    REAL(KIND=PRE), DIMENSION(328, 46, 2, 2), INTENT(IN) :: APD

    INTEGER,                                  INTENT(IN) :: NEXP
    REAL(KIND=PRE), DIMENSION(NEXP),          INTENT(IN) :: AMBDA, AR

    ! Outputs
    COMPLEX(KIND=PRE),               INTENT(OUT) :: SP  ! Integral of the Green function over the panel.
    COMPLEX(KIND=PRE), DIMENSION(3), INTENT(OUT) :: VSP_SYM, VSP_ANTISYM ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    INTEGER                              :: KE
    REAL(KIND=PRE)                       :: AMH, AKH, A, COF1, COF2, COF3, COF4
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
    CALL COMPUTE_INTEGRALS_WRT_THETA(XI(:), XJ(:), wavenumber, XR, XZ, APD, FS(1), VS(:, 1))

    PSR(1) = PI/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))

    ! 1.b Shift and reflect XI and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) =  X0J(3)
    CALL COMPUTE_INTEGRALS_WRT_THETA(XI(:), XJ(:), wavenumber, XR, XZ, APD, FS(2), VS(:, 2))
    VS(3, 2) = -VS(3, 2) ! Reflection of the output vector

    PSR(2) = PI/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))

    ! 1.c Shift and reflect XJ and compute another value of the Green function
    XI(3) =  X0I(3)
    XJ(3) = -X0J(3) - 2*depth
    CALL COMPUTE_INTEGRALS_WRT_THETA(XI(:), XJ(:), wavenumber, XR, XZ, APD, FS(3), VS(:, 3))

    PSR(3) = PI/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))

    ! 1.d Shift and reflect both XI and XJ and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) = -X0J(3) - 2*depth
    CALL COMPUTE_INTEGRALS_WRT_THETA(XI(:), XJ(:), wavenumber, XR, XZ, APD, FS(4), VS(:, 4))
    VS(3, 4) = -VS(3, 4) ! Reflection of the output vector

    PSR(4) = PI/(wavenumber*SQRT(R**2+(XI(3)+XJ(3))**2))

    ! Add up the results of the four problems
    SP               = SUM(FS(1:4)) - SUM(PSR(1:4))
    VSP_SYM(1:3)     = VS(1:3, 1) + VS(1:3, 4)
    VSP_ANTISYM(1:3) = VS(1:3, 2) + VS(1:3, 3)

    ! Multiply by some coefficients
    AMH  = wavenumber*depth
    AKH  = AMH*TANH(AMH)
    A    = (AMH+AKH)**2/(depth*(AMH**2-AKH**2+AKH))
    COF1 = A/(2*PI)
    COF2 = A/2
    COF3 = wavenumber*COF1
    COF4 = wavenumber*COF2

    SP          = CMPLX(REAL(SP)*COF1,          AIMAG(SP)*COF2, KIND=PRE)
    VSP_ANTISYM = CMPLX(REAL(VSP_ANTISYM)*COF3, AIMAG(VSP_ANTISYM)*COF4, KIND=PRE)
    VSP_SYM     = CMPLX(REAL(VSP_SYM)*COF3,     AIMAG(VSP_SYM)*COF4, KIND=PRE)

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

  ! =====================================================================

  SUBROUTINE BUILD_MATRICES_WAVE_SOURCE  &
      (nb_faces_1, centers_1, normals_1, &
      nb_faces_2,                        &
      centers_2, areas_2,                &
      wavenumber, depth,                 &
      XR, XZ, APD,                       &
      NEXP, AMBDA, AR,                   &
      same_body,                         &
      S, V)

    ! Mesh data
    INTEGER,                                  INTENT(IN) :: nb_faces_1, nb_faces_2
    REAL(KIND=PRE), DIMENSION(nb_faces_1, 3), INTENT(IN) :: normals_1, centers_1
    REAL(KIND=PRE), DIMENSION(nb_faces_2, 3), INTENT(IN) :: centers_2
    REAL(KIND=PRE), DIMENSION(nb_faces_2),    INTENT(IN) :: areas_2

    REAL(KIND=PRE),                           INTENT(IN) :: wavenumber, depth

    ! Tabulated integrals
    REAL(KIND=PRE), DIMENSION(328),           INTENT(IN) :: XR
    REAL(KIND=PRE), DIMENSION(46),            INTENT(IN) :: XZ
    REAL(KIND=PRE), DIMENSION(328, 46, 2, 2), INTENT(IN) :: APD

    ! Prony decomposition for finite depth
    INTEGER,                                  INTENT(IN) :: NEXP
    REAL(KIND=PRE), DIMENSION(NEXP),          INTENT(IN) :: AMBDA, AR

    ! Trick to save some time
    LOGICAL,                                  INTENT(IN) :: same_body

    ! Output
    COMPLEX(KIND=PRE), DIMENSION(nb_faces_1, nb_faces_2), INTENT(OUT) :: S
    COMPLEX(KIND=PRE), DIMENSION(nb_faces_1, nb_faces_2), INTENT(OUT) :: V

    ! Local variables
    INTEGER                         :: I, J
    COMPLEX(KIND=PRE)               :: SP2
    COMPLEX(KIND=PRE), DIMENSION(3) :: VSP2_SYM, VSP2_ANTISYM

    IF (SAME_BODY) THEN
      ! If we are computing the influence of some cells upon themselves, the resulting matrices have some symmetries.
      ! This is due to the symmetry of the Green function, and the way the integral on the face is approximated.
      ! (More precisely, the Green function is symmetric and its derivative is the sum of a symmetric part and an anti-symmetric
      ! part.)

      DO I = 1, nb_faces_1
        !$OMP PARALLEL DO PRIVATE(J, SP2, VSP2_SYM, VSP2_ANTISYM)
        DO J = I, nb_faces_2

          IF (depth == INFINITE_DEPTH) THEN
            CALL WAVE_PART_INFINITE_DEPTH &
              (wavenumber,                &
              centers_1(I, :),            &
              centers_2(J, :),            &
              XR, XZ, APD,                &
              SP2, VSP2_SYM               &
              )
            VSP2_ANTISYM(:) = ZERO
          ELSE
            CALL WAVE_PART_FINITE_DEPTH   &
              (wavenumber,                &
              centers_1(I, :),            &
              centers_2(J, :),            &
              depth,                      &
              XR, XZ, APD,                &
              NEXP, AMBDA, AR,            &
              SP2, VSP2_SYM, VSP2_ANTISYM &
              )
          END IF

          S(I, J) = -1/(4*PI) * SP2*areas_2(J)
          V(I, J) = -1/(4*PI) * DOT_PRODUCT(normals_1(I, :),         &
                                            VSP2_SYM + VSP2_ANTISYM) &
                                            *areas_2(J)

          IF (.NOT. I==J) THEN
            VSP2_SYM(1:2) = -VSP2_SYM(1:2)
            S(J, I) = -1/(4*PI) * SP2*areas_2(I)
            V(J, I) = -1/(4*PI) * DOT_PRODUCT(normals_1(J, :),         &
                                              VSP2_SYM - VSP2_ANTISYM) &
                                              *areas_2(I)
          END IF

        END DO
        !$OMP END PARALLEL DO
      END DO

    ELSE
      ! General case: if we are computing the influence of a some cells on other cells, we have to compute all the coefficients.

      DO I = 1, nb_faces_1
        !$OMP PARALLEL DO PRIVATE(J, SP2, VSP2_SYM, VSP2_ANTISYM)
        DO J = 1, nb_faces_2

          IF (depth == INFINITE_DEPTH) THEN
            CALL WAVE_PART_INFINITE_DEPTH &
              (wavenumber,                &
              centers_1(I, :),            &
              centers_2(J, :),            &
              XR, XZ, APD,                &
              SP2, VSP2_SYM               &
              )
            VSP2_ANTISYM(:) = ZERO 
          ELSE
            CALL WAVE_PART_FINITE_DEPTH   &
              (wavenumber,                &
              centers_1(I, :),            &
              centers_2(J, :),            &
              depth,                      &
              XR, XZ, APD,                &
              NEXP, AMBDA, AR,            &
              SP2, VSP2_SYM, VSP2_ANTISYM &
              )
          END IF

          S(I, J) = -1/(4*PI) * SP2*areas_2(J)                                ! Green function
          V(I, J) = -1/(4*PI) * DOT_PRODUCT(normals_1(I, :),         &
                                            VSP2_SYM + VSP2_ANTISYM) &
                                            *areas_2(J) ! Gradient of the Green function

        END DO
        !$OMP END PARALLEL DO
      END DO
   END IF

  END SUBROUTINE

  ! =====================================================================

END MODULE GREEN_WAVE
