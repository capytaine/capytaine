MODULE Green_2

  IMPLICIT NONE

  REAL, PARAMETER :: PI = 3.141592653588979 ! π
  REAL, PARAMETER :: DPI  = 6.283185307179586 ! 2π
  REAL, PARAMETER :: DPI2 = 19.73920880217871 ! 2π²

  INTEGER, PARAMETER      :: NPINTE = 251
  INTEGER, PARAMETER      :: IR = 328
  INTEGER, PARAMETER      :: JZ = 46
  REAL, DIMENSION(IR)     :: XR
  REAL, DIMENSION(JZ)     :: XZ
  REAL, DIMENSION(IR, JZ) :: APD1X, APD1Z, APD2X, APD2Z

CONTAINS

  COMPLEX FUNCTION GG(Z, CEX)
    ! Estimation of ∫_z^∞ exp(-t)/t dt
    ! See p.367 of G. Delhommeau thesis (referenced as [Del]).

    COMPLEX, INTENT(IN) :: Z, CEX
    COMPLEX             :: Y

    IF (REAL(Z) < -16.0) THEN                                      ! Case 1 p. 368 in [Del]
      Y = 1./Z
      GG = Y*(1.+Y*(-1.+Y*(2.+Y*(-6.+Y*(24.+Y*(-120.))))))
    ELSE IF (ABS(AIMAG(Z)) > 10.0) THEN                            ! Case 3 p. 368 in [Del]
      GG = 0.711093/(Z+0.415775)+0.278518/(Z+2.29428)+0.010389/(Z+6.2900)
      IF (AIMAG(Z) < 0) THEN
        GG = GG-(0., 3.14159265)*CEX
      ELSE
        GG = GG+(0., 3.14159265)*CEX
      END IF
    ELSE IF (REAL(Z) > -0.5) THEN                                  ! Case 2 p. 368 in [Del]
      GG = -(CLOG(Z)*(.1E+01+Z*(0.23721365E+00+Z*(0.206543E-01+Z*(0.763297E-03+    &
        Z*0.97087007E-05))))+0.5772156649E+00*(0.99999207E+00+                    &
        Z*(-0.149545886E+01+Z*(0.41806426E-01+Z*(-0.3000591E-01+                  &
        Z*(0.19387339E-02+Z*(-0.51801555E-03)))))))/(0.1E+01+                     &
        Z*(-0.76273617E+00+Z*(0.28388363E+00+Z*(-0.66786033E-01+Z*(0.12982719E-01 &
        +Z*(-0.8700861E-03+Z*0.2989204E-03))))))
      IF (AIMAG(Z) < 0) THEN
        GG = GG-(0., 3.14159265)*CEX
      ELSE
        GG = GG+(0., 3.14159265)*CEX
      END IF
    ELSE                                                           ! Case 4 p. 369 in [Del]
      IF (AIMAG(Z) < 0) THEN
        GG = ((((((( (1.000000, 1.3935496E-06)*Z+ (15.82958, -20.14222))  &
          *Z+ (-70.52863, -227.9511))*Z+ (-985.4221, -226.6272))*Z        &
          + (-1202.318, 1580.907))*Z+ (953.2441, 1342.447))*Z             &
          + (417.3716, -196.6665))*Z+ (-9.881266, -24.24952))/            &
          (((((((( (1.000000, 0.0000000E+00)*Z+ (16.83184, -20.14481))*Z  &
          + (-55.66969, -248.1167))*Z+ (-1068.640, -434.4707))*Z          &
          + (-2082.250, 1522.471))*Z+ (383.3455, 2730.378))*Z             &
          + (1216.791, 351.7189))*Z+ (115.3926, -161.2647))*Z             &
          + (-3.777369, -4.510900))-(0., 3.14159265)*CEX
      ELSE
        GG = ((((((( (1.000000, -1.3935496E-06)*Z+ (15.82958, 20.14222))  &
          *Z+ (-70.52863, 227.9511))*Z+ (-985.4221, 226.6272))*Z          &
          + (-1202.318, -1580.907))*Z+ (953.2441, -1342.447))*Z           &
          + (417.3716, 196.6665))*Z+ (-9.881266, 24.24952))/              &
          (((((((( (1.000000, 0.0000000E+00)*Z+ (16.83184, 20.14481))*Z   &
          + (-55.66969, 248.1167))*Z+ (-1068.640, 434.4707))*Z            &
          + (-2082.250, -1522.471))*Z+ (383.3455, -2730.378))*Z           &
          + (1216.791, -351.7189))*Z+ (115.3926, 161.2647))*Z             &
          + (-3.777369, 4.510900))+(0., 3.14159265)*CEX
      END IF
    END IF
  END FUNCTION

  REAL FUNCTION PL2(U1,U2,U3,XU)
    REAL::U1,U2,U3,XU
    PL2=((XU-U1)*(XU-U2))/((U3-U1)*(U3-U2))
    RETURN
  END FUNCTION

  SUBROUTINE INITIALIZE_GREEN()
    ! Initialize XR, XZ, APD1X, APD2X, APD1Z, APD2Z
    ! Those parameters are independent of the depth and the frequency.
    ! Thus, they are initialized only once at the beginning of the execution of the code.
    ! Other parameters are initialized in LISC below.

    ! Local variables
    INTEGER :: I, J, K
    REAL :: QQT(NPINTE), CQT(NPINTE)
    REAL :: CT
    COMPLEX :: C1, C2, ZIK, CEX

    ! Initialize XZ
    DO J = 1, JZ
      XZ(J) = -AMIN1(10**(J/5.0-6), 10**(J/8.0-4.5), 16.)
    END DO

    ! Initialize XR
    XR(1) = 0.0
    DO I = 2, IR
      IF (I < 40) THEN
        XR(I) = AMIN1(10**((I-1.0)/5-6), 4.0/3.0 + ABS(I-32)/3.0)
      ELSE
        XR(I) = 4.0/3.0 + ABS(I-32)/3.0
      ENDIF
    END DO

    ! Initialize QQT and CQT
    DO K = 1, NPINTE
      QQT(K) = -PI/2 + (K-1.0)/(NPINTE-1.0)*PI
      IF ((K <= 1) .OR. (K >= NPINTE)) THEN
        CQT(K) = PI/(3*(NPINTE-1))
      ELSEIF (MOD(K,2)==0) THEN
        CQT(K) = 4.0/(3*(NPINTE-1))*PI
      ELSE
        CQT(K) = 2.0/(3*(NPINTE-1))*PI
      ENDIF
    ENDDO

    ! Initialize APD..
    APD1X(:, :) = 0.0
    APD1Z(:, :) = 0.0
    APD2X(:, :) = 0.0
    APD2Z(:, :) = 0.0
    DO J = 1, JZ
      DO I = 1, IR
        DO K = 1, NPINTE
          CT = COS(QQT(K))
          ZIK = CMPLX(XZ(J), XR(I)*CT)
          IF (REAL(ZIK) <= -30.0) THEN
            CEX = (0.0, 0.0)
          ELSE
            CEX = CEXP(ZIK)
          ENDIF
          C1 = CQT(K)*(GG(ZIK, CEX) - 1.0/ZIK)
          C2 = CQT(K)*CEX
          APD1X(I, J) = APD1X(I, J) + CT*AIMAG(C1)
          APD1Z(I, J) = APD1Z(I, J) + REAL(C1)
          APD2X(I, J) = APD2X(I, J) + CT*AIMAG(C2)
          APD2Z(I, J) = APD2Z(I, J) + REAL(C2)
        END DO
      END DO
    END DO

    RETURN

  END SUBROUTINE INITIALIZE_GREEN

  SUBROUTINE COMPUTE_S2(XI, XJ, depth, wavenumber, FS, VS)

    ! Inputs
    REAL, DIMENSION(3),    INTENT(IN)  :: XI, XJ
    REAL,                  INTENT(IN)  :: depth, wavenumber

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

        PD1Z = DOT_PRODUCT(XL, MATMUL(APD1Z(KI-1:KI+1, KJ-1:KJ+1), ZL))
        PD2Z = DOT_PRODUCT(XL, MATMUL(APD2Z(KI-1:KI+1, KJ-1:KJ+1), ZL))

        IF (RRR > 1e-5) THEN
          PD1X = DOT_PRODUCT(XL, MATMUL(APD1X(KI-1:KI+1, KJ-1:KJ+1), ZL))
          PD2X = DOT_PRODUCT(XL, MATMUL(APD2X(KI-1:KI+1, KJ-1:KJ+1), ZL))
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

  SUBROUTINE VNSINFD             &
      (wavenumber, X0I, X0J, AJ, &
      SP, SM, VSP, VSM)
    ! Compute the frequency-dependent part of the Green function in the infinite depth case.

    ! Inputs
    REAL,               INTENT(IN)  :: wavenumber, AJ
    REAL, DIMENSION(3), INTENT(IN)  :: X0I   ! Coordinates of the source point
    REAL, DIMENSION(3), INTENT(IN)  :: X0J   ! Coordinates of the center of the integration panel

    ! Outputs
    COMPLEX,               INTENT(OUT) :: SP, SM   ! Integral of the Green function over the panel.
    COMPLEX, DIMENSION(3), INTENT(OUT) :: VSP, VSM ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    REAL                               :: ADPI, ADPI2, AKDPI, AKDPI2
    REAL, DIMENSION(3)       :: XI
    COMPLEX, DIMENSION(1)    :: FS
    COMPLEX, DIMENSION(3, 1) :: VS

    XI(:) = X0I(:)
    ! XI(3) = MIN(X0I(3), -1e-5*Mesh%xy_diameter)
    CALL COMPUTE_S2(XI, X0J, 0.0, wavenumber, FS(1), VS(:, 1))

    ! IF (Mesh%Isym == 0) THEN
      SP       = FS(1)
      VSP(1:3) = VS(1:3, 1)
      SM       = CMPLX(0.0, 0.0)
      VSM      = CMPLX(0.0, 0.0)

    ! ELSE IF (Mesh%Isym == 1) THEN
    !   ! Reflect the source point across the (xOz) plane and compute another coefficient
    !   XI(2) = -X0I(2)
    !   CALL COMPUTE_S2(XI, X0J, wavenumber, FS(2), VS(:, 2))
    !   VS(2, 2) = -VS(2, 2) ! Reflection of the output vector

    !   ! Assemble the two results
    !   SP       = FS(1)      + FS(2)
    !   VSP(1:3) = VS(1:3, 1) + VS(1:3, 2)
    !   SM       = FS(1)      - FS(2)
    !   VSM(1:3) = VS(1:3, 1) - VS(1:3, 2)
    ! END IF

    ADPI2  = wavenumber*AJ/DPI2
    ADPI   = wavenumber*AJ/DPI
    AKDPI2 = wavenumber**2*AJ/DPI2
    AKDPI  = wavenumber**2*AJ/DPI

    SP  = CMPLX(REAL(SP)*ADPI2,   AIMAG(SP)*ADPI)
    VSP = CMPLX(REAL(VSP)*AKDPI2, AIMAG(VSP)*AKDPI)

    ! IF (Mesh%ISym == Y_SYMMETRY) THEN
    !   SM  = CMPLX(REAL(SM)*ADPI2,   AIMAG(SM)*ADPI)
    !   VSM = CMPLX(REAL(VSM)*AKDPI2, AIMAG(VSM)*AKDPI)
    ! END IF

    RETURN
  END SUBROUTINE VNSINFD
END MODULE Green_2
