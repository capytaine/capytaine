MODULE Green_2

  USE Initialize_Green_2
  USE Green_1, ONLY: COMPUTE_ASYMPTOTIC_S0

  IMPLICIT NONE

  REAL, PARAMETER :: DPI  = 6.283185307179586 ! 2π
  REAL, PARAMETER :: DPI2 = 19.73920880217871 ! 2π²

CONTAINS

  REAL FUNCTION PL2(U1,U2,U3,XU)
    REAL::U1,U2,U3,XU
    PL2=((XU-U1)*(XU-U2))/((U3-U1)*(U3-U2))
    RETURN
  END FUNCTION

  SUBROUTINE COMPUTE_S2(XI, XJ, depth, wavenumber, &
                        XR, XZ, APD,               &
                        FS, VS)

    ! Inputs
    REAL, DIMENSION(3),    INTENT(IN)  :: XI, XJ
    REAL,                  INTENT(IN)  :: depth, wavenumber

    REAL, DIMENSION(328), INTENT(IN) :: XR
    REAL, DIMENSION(46),  INTENT(IN) :: XZ
    REAL, DIMENSION(328, 46, 2, 2), INTENT(IN) :: APD

    ! Outputs
    COMPLEX,               INTENT(OUT) :: FS
    COMPLEX, DIMENSION(3), INTENT(OUT) :: VS

    ! Local variables
    INTEGER                            :: KI, KJ
    REAL                               :: RRR, AKR, ZZZ, AKZ, DD, PSURR
    REAL                               :: SIK, CSK, SQ, EPZ
    REAL                               :: PD1X, PD2X, PD1Z, PD2Z
    REAL, DIMENSION(3)                 :: XL, ZL

    RRR = NORM2(XI(1:2) - XJ(1:2))
    AKR = wavenumber*RRR

    ZZZ = XI(3) + XJ(3)
    AKZ = wavenumber*ZZZ

    DD  = SQRT(RRR**2 + ZZZ**2)

    IF ((DD > 1e-5) .AND. (wavenumber > 0)) THEN
      PSURR = PI/(wavenumber*DD)**3
    ELSE
      PSURR = 0.0
    ENDIF

    ! IF (AKZ > -1.5e-6) THEN
    !   WRITE(*,*)'AKZ < -1.5 E-6' ! Not a very explicit warning...
    ! END IF

    IF (AKZ > -16) THEN             !   -16 < AKZ < -1.5e-6

      !================================================
      ! Evaluate PDnX and PDnZ depending on AKZ and AKR
      !================================================

      IF (AKR < 99.7) THEN          !     0 < AKR < 99.7

        IF (AKZ < -1e-2) THEN       !   -16 < AKZ < -1e-2
          KJ = INT(8*(ALOG10(-AKZ)+4.5))
        ELSE                        ! -1e-2 < AKZ < -1.5e-6
          KJ = INT(5*(ALOG10(-AKZ)+6))
        ENDIF
        KJ = MAX(MIN(KJ, 45), 2)

        IF (AKR < 1) THEN           !     0 < AKR < 1
          KI = INT(5*(ALOG10(AKR+1e-20)+6)+1)
        ELSE                        !     1 < AKR < 99.7
          KI = INT(3*AKR+28)
        ENDIF
        KI = MAX(MIN(KI, 327), 2)

        XL(1) = PL2(XR(KI),   XR(KI+1), XR(KI-1), AKR)
        XL(2) = PL2(XR(KI+1), XR(KI-1), XR(KI),   AKR)
        XL(3) = PL2(XR(KI-1), XR(KI),   XR(KI+1), AKR)
        ZL(1) = PL2(XZ(KJ),   XZ(KJ+1), XZ(KJ-1), AKZ)
        ZL(2) = PL2(XZ(KJ+1), XZ(KJ-1), XZ(KJ),   AKZ)
        ZL(3) = PL2(XZ(KJ-1), XZ(KJ),   XZ(KJ+1), AKZ)

        PD1Z = DOT_PRODUCT(XL, MATMUL(APD(KI-1:KI+1, KJ-1:KJ+1, 1, 2), ZL))
        PD2Z = DOT_PRODUCT(XL, MATMUL(APD(KI-1:KI+1, KJ-1:KJ+1, 2, 2), ZL))

        IF (RRR > 1e-5) THEN
          PD1X = DOT_PRODUCT(XL, MATMUL(APD(KI-1:KI+1, KJ-1:KJ+1, 1, 1), ZL))
          PD2X = DOT_PRODUCT(XL, MATMUL(APD(KI-1:KI+1, KJ-1:KJ+1, 2, 1), ZL))
        END IF

      ELSE  ! 99.7 < AKR

        EPZ  = EXP(AKZ)
        SQ   = SQRT(2*PI/AKR)
        CSK  = COS(AKR-PI/4)
        SIK  = SIN(AKR-PI/4)

        PD1Z = PSURR*AKZ - PI*EPZ*SQ*SIK
        PD2Z =                EPZ*SQ*CSK

        IF (RRR > 1e-5) THEN
          ! PD1X=-PSURR*AKR-PI*EPZ*SQ*(CSK-0.5/AKR*SIK) ! correction par GD le 17/09/2010
          PD1X = PI*EPZ*SQ*(CSK - 0.5*SIK/AKR) - PSURR*AKR
          PD2X =    EPZ*SQ*(SIK + 0.5*CSK/AKR)
        END IF

      ENDIF

      !====================================
      ! Deduce FS ans VS from PDnX and PDnZ
      !====================================

      FS    = -CMPLX(PD1Z, PD2Z)
      IF (depth == 0.0) THEN
        VS(3) = -CMPLX(PD1Z-PSURR*AKZ, PD2Z)
      ELSE
        VS(3) = -CMPLX(PD1Z, PD2Z)
      END IF

      IF (RRR > 1e-5) THEN
        IF (depth == 0.0) THEN
          VS(1) = (XI(1) - XJ(1))/RRR * CMPLX(PD1X+PSURR*AKR, PD2X)
          VS(2) = (XI(2) - XJ(2))/RRR * CMPLX(PD1X+PSURR*AKR, PD2X)
        ELSE
          VS(1) = (XI(1) - XJ(1))/RRR * CMPLX(PD1X, PD2X)
          VS(2) = (XI(2) - XJ(2))/RRR * CMPLX(PD1X, PD2X)
        END IF
      ELSE
        VS(1:2) = CMPLX(0.0, 0.0)
      END IF

    ELSE ! AKZ < -16
      FS      = CMPLX(-PSURR*AKZ, 0.0)
      VS(1:3) = CMPLX(0.0, 0.0)
    ENDIF

    RETURN
  END SUBROUTINE COMPUTE_S2

  ! =========================

  SUBROUTINE VNSINFD         &
      (wavenumber, X0I, X0J, &
      XR, XZ, APD,           &
      SP, VSP)
    ! Compute the frequency-dependent part of the Green function in the infinite depth case.

    ! Inputs
    REAL,               INTENT(IN)  :: wavenumber
    REAL, DIMENSION(3), INTENT(IN)  :: X0I   ! Coordinates of the source point
    REAL, DIMENSION(3), INTENT(IN)  :: X0J   ! Coordinates of the center of the integration panel

    REAL, DIMENSION(328), INTENT(IN) :: XR
    REAL, DIMENSION(46),  INTENT(IN) :: XZ
    REAL, DIMENSION(328, 46, 2, 2), INTENT(IN) :: APD


    ! Outputs
    COMPLEX,               INTENT(OUT) :: SP  ! Integral of the Green function over the panel.
    COMPLEX, DIMENSION(3), INTENT(OUT) :: VSP ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    REAL                     :: ADPI, ADPI2, AKDPI, AKDPI2
    REAL, DIMENSION(3)       :: XI

    XI(:) = X0I(:)
    ! XI(3) = MIN(X0I(3), -1e-5*Mesh%xy_diameter)
    CALL COMPUTE_S2(XI, X0J, 0.0, wavenumber, XR, XZ, APD, SP, VSP(:))

    ADPI2  = wavenumber/DPI2
    ADPI   = wavenumber/DPI
    AKDPI2 = wavenumber**2/DPI2
    AKDPI  = wavenumber**2/DPI

    SP  = CMPLX(REAL(SP)*ADPI2,   AIMAG(SP)*ADPI)
    VSP = CMPLX(REAL(VSP)*AKDPI2, AIMAG(VSP)*AKDPI)

    RETURN
  END SUBROUTINE VNSINFD

  ! ======================

  SUBROUTINE VNSFD &
      (wavenumber, X0I, X0J, depth, &
      XR, XZ, APD,  &
      AMBDA, AR, NEXP, &
      SP, VSP_SYM, VSP_ANTISYM)
    ! Compute the frequency-dependent part of the Green function in the finite depth case.

    ! Inputs
    REAL,               INTENT(IN)  :: wavenumber, depth
    REAL, DIMENSION(3), INTENT(IN)  :: X0I   ! Coordinates of the source point
    REAL, DIMENSION(3), INTENT(IN)  :: X0J   ! Coordinates of the center of the integration panel

    REAL, DIMENSION(328), INTENT(IN) :: XR
    REAL, DIMENSION(46),  INTENT(IN) :: XZ
    REAL, DIMENSION(328, 46, 2, 2), INTENT(IN) :: APD

    INTEGER, INTENT(IN) :: NEXP
    REAL, DIMENSION(31), INTENT(INOUT)  :: AMBDA, AR

    ! Outputs
    COMPLEX,               INTENT(OUT) :: SP  ! Integral of the Green function over the panel.
    COMPLEX, DIMENSION(3), INTENT(OUT) :: VSP_SYM, VSP_ANTISYM ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    INTEGER                     :: KE
    REAL                        :: AMH, AKH, A, COF1, COF2, COF3, COF4
    REAL                        :: AQT, RRR
    REAL, DIMENSION(3)          :: XI, XJ
    REAL, DIMENSION(4)          :: FTS, PSR
    REAL, DIMENSION(3, 4)       :: VTS
    COMPLEX, DIMENSION(4)       :: FS
    COMPLEX, DIMENSION(3, 4)    :: VS

    !========================================
    ! Part 1: Solve 4 infinite depth problems
    !========================================

    XI(:) = X0I(:)
    XJ(:) = X0J(:)

    ! Distance in xOy plane
    RRR = NORM2(XI(1:2) - XJ(1:2))

    ! 1.a First infinite depth problem
    CALL COMPUTE_S2(XI(:), XJ(:), depth, wavenumber, XR, XZ, APD, FS(1), VS(:, 1))

    PSR(1) = PI/(wavenumber*SQRT(RRR**2+(XI(3)+XJ(3))**2))

    ! 1.b Shift and reflect XI and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) =  X0J(3)
    CALL COMPUTE_S2(XI(:), XJ(:), depth, wavenumber, XR, XZ, APD, FS(2), VS(:, 2))
    VS(3, 2) = -VS(3, 2) ! Reflection of the output vector

    PSR(2) = PI/(wavenumber*SQRT(RRR**2+(XI(3)+XJ(3))**2))

    ! 1.c Shift and reflect XJ and compute another value of the Green function
    XI(3) =  X0I(3)
    XJ(3) = -X0J(3) - 2*depth
    CALL COMPUTE_S2(XI(:), XJ(:), depth, wavenumber, XR, XZ, APD, FS(3), VS(:, 3))

    PSR(3) = PI/(wavenumber*SQRT(RRR**2+(XI(3)+XJ(3))**2))

    ! 1.d Shift and reflect both XI and XJ and compute another value of the Green function
    XI(3) = -X0I(3) - 2*depth
    XJ(3) = -X0J(3) - 2*depth
    CALL COMPUTE_S2(XI(:), XJ(:), depth, wavenumber, XR, XZ, APD, FS(4), VS(:, 4))
    VS(3, 4) = -VS(3, 4) ! Reflection of the output vector

    PSR(4) = PI/(wavenumber*SQRT(RRR**2+(XI(3)+XJ(3))**2))

    ! Add up the results of the four problems
    SP               = -SUM(FS(1:4)) - SUM(PSR(1:4))
    VSP_SYM(1:3)     = -VS(1:3, 1) - VS(1:3, 4)
    VSP_ANTISYM(1:3) = -VS(1:3, 2) - VS(1:3, 3)

    ! Multiply by some coefficients
    AMH  = wavenumber*depth
    AKH  = AMH*TANH(AMH)
    A    = (AMH+AKH)**2/(depth*(AMH**2-AKH**2+AKH))
    COF1 = -A/(8*PI**2)
    COF2 = -A/(8*PI)
    COF3 = wavenumber*COF1
    COF4 = wavenumber*COF2

    SP          = CMPLX(REAL(SP)*COF1,          AIMAG(SP)*COF2)
    VSP_ANTISYM = CMPLX(REAL(VSP_ANTISYM)*COF3, AIMAG(VSP_ANTISYM)*COF4)
    VSP_SYM     = CMPLX(REAL(VSP_SYM)*COF3,     AIMAG(VSP_SYM)*COF4)

    !=====================================================
    ! Part 2: Integrate (NEXP+1)×4 terms of the form 1/MM'
    !=====================================================

    AMBDA(NEXP+1) = 0
    AR(NEXP+1)    = 2

    DO KE = 1, NEXP+1
      XI(:) = X0I(:)

      ! 2.a Shift observation point and compute integral
      XI(3) =  X0I(3) + depth*AMBDA(KE) - 2*depth
      CALL COMPUTE_ASYMPTOTIC_S0(XI(:), X0J(:), 1.0, FTS(1), VTS(:, 1))

      ! 2.b Shift and reflect observation point and compute integral
      XI(3) = -X0I(3) - depth*AMBDA(KE)
      CALL COMPUTE_ASYMPTOTIC_S0(XI(:), X0J(:), 1.0, FTS(2), VTS(:, 2))
      VTS(3, 2) = -VTS(3, 2) ! Reflection of the output vector

      ! 2.c Shift and reflect observation point and compute integral
      XI(3) = -X0I(3) + depth*AMBDA(KE) - 4*depth
      CALL COMPUTE_ASYMPTOTIC_S0(XI(:), X0J(:), 1.0, FTS(3), VTS(:, 3))
      VTS(3, 3) = -VTS(3, 3) ! Reflection of the output vector

      ! 2.d Shift observation point and compute integral
      XI(3) =  X0I(3) - depth*AMBDA(KE) + 2*depth
      CALL COMPUTE_ASYMPTOTIC_S0(XI(:), X0J(:), 1.0, FTS(4), VTS(:, 4))

      AQT = -AR(KE)/(8*PI)

      ! Add all the contributions
      SP               = SP               + AQT*SUM(FTS(1:4))
      VSP_ANTISYM(1:3) = VSP_ANTISYM(1:3) + AQT*(VTS(1:3, 1) + VTS(1:3, 4))
      VSP_SYM(1:3)     = VSP_SYM(1:3)     + AQT*(VTS(1:3, 2) + VTS(1:3, 3))

    END DO

    RETURN
  END SUBROUTINE

  ! =====================================================================

  SUBROUTINE BUILD_MATRIX_2(            &
      nb_faces_1, centers_1, normals_1, &
      nb_faces_2,                       &
      centers_2, areas_2,               &
      wavenumber, depth,                &
      XR, XZ, APD,                      &
      AMBDA, AR, NEXP,                  &
      same_body,                        &
      S, V)

    INTEGER,                              INTENT(IN) :: nb_faces_1, nb_faces_2
    REAL,    DIMENSION(nb_faces_1, 3),    INTENT(IN) :: normals_1, centers_1
    REAL,    DIMENSION(nb_faces_2, 3),    INTENT(IN) :: centers_2
    REAL,    DIMENSION(nb_faces_2),       INTENT(IN) :: areas_2
    REAL,                                 INTENT(IN) :: wavenumber, depth

    REAL, DIMENSION(328), INTENT(IN) :: XR
    REAL, DIMENSION(46),  INTENT(IN) :: XZ
    REAL, DIMENSION(328, 46, 2, 2), INTENT(IN) :: APD
    INTEGER,                              INTENT(IN) :: NEXP
    REAL,    DIMENSION(31),               INTENT(INOUT) :: AMBDA, AR

    LOGICAL,                              INTENT(IN) :: same_body

    COMPLEX, DIMENSION(nb_faces_1, nb_faces_2), INTENT(OUT) :: S
    COMPLEX, DIMENSION(nb_faces_1, nb_faces_2), INTENT(OUT) :: V

    ! Local variables
    INTEGER               :: I, J
    COMPLEX               :: SP2
    COMPLEX, DIMENSION(3) :: VSP2_SYM, VSP2_ANTISYM

    IF (SAME_BODY) THEN
      ! Use the symmetry of SP2 and VSP2

      DO I = 1, nb_faces_1
        DO J = I, nb_faces_2

          IF (depth == 0.0) THEN
            CALL VNSINFD                  &
              (wavenumber,                &
              centers_1(I, :),            &
              centers_2(J, :),            &
              XR, XZ, APD,                &
              SP2, VSP2_SYM               &
              )
            VSP2_ANTISYM(:) = 0
          ELSE
            CALL VNSFD                    &
              (wavenumber,                &
              centers_1(I, :),            &
              centers_2(J, :),            &
              depth,                      &
              XR, XZ, APD,                &
              AMBDA, AR, NEXP,            &
              SP2, VSP2_SYM, VSP2_ANTISYM &
              )
          END IF

          S(I, J) = SP2*areas_2(J)
          V(I, J) = DOT_PRODUCT(normals_1(I, :),         &
                                VSP2_SYM + VSP2_ANTISYM) &
                                *areas_2(J)

          IF (.NOT. I==J) THEN
            VSP2_SYM(1:2) = -VSP2_SYM(1:2)
            S(J, I) = SP2*areas_2(I)
            V(J, I) = DOT_PRODUCT(normals_1(J, :),         &
                                  VSP2_SYM - VSP2_ANTISYM) &
                                  *areas_2(I)
          END IF

        END DO
      END DO

    ELSE

      DO I = 1, nb_faces_1
        DO J = 1, nb_faces_2

          IF (depth == 0.0) THEN
            CALL VNSINFD                  &
              (wavenumber,                &
              centers_1(I, :),            &
              centers_2(J, :),            &
              XR, XZ, APD,                &
              SP2, VSP2_SYM               &
              )
            VSP2_ANTISYM(:) = 0
          ELSE
            CALL VNSFD                    &
              (wavenumber,                &
              centers_1(I, :),            &
              centers_2(J, :),            &
              depth,                      &
              XR, XZ, APD,                &
              AMBDA, AR, NEXP,            &
              SP2, VSP2_SYM, VSP2_ANTISYM &
              )
          END IF

          S(I, J) = SP2*areas_2(J)                                ! Green function
          V(I, J) = DOT_PRODUCT(normals_1(I, :),         &
                                VSP2_SYM + VSP2_ANTISYM) &
                                *areas_2(J) ! Gradient of the Green function

        END DO
      END DO
   END IF

  END SUBROUTINE

  ! =====================================================================

END MODULE Green_2
