MODULE Initialize_Green_2

  IMPLICIT NONE

  PUBLIC :: GG
  PUBLIC :: INITIALIZE_GREEN
  PUBLIC :: FF

  PUBLIC :: LISC   ! Initialization of AMBDA and AR
  PUBLIC :: EXPORS ! Called by LISC
  PUBLIC :: MCAS   ! Called by EXPORS
  PUBLIC :: HOUSRS ! Called by EXPORS and MCAS
  PUBLIC :: SPRBM  ! Called by EXPORS
  PUBLIC :: SPQFB  ! Called by SPRBM

  REAL, PARAMETER :: PI = 3.141592653588979 ! π
  COMPLEX, PARAMETER :: II = (0, 1)         ! Imaginary unit

  ! Independent of Omega
  INTEGER, PARAMETER      :: NPINTE = 251
  INTEGER, PARAMETER      :: IR = 328
  INTEGER, PARAMETER      :: JZ = 46

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
    ELSE IF (REAL(Z) > -0.5) THEN                                  ! Case 2 p. 368 in [Del]
      GG = -(LOG(Z)*(.1E+01+Z*(0.23721365E+00+Z*(0.206543E-01+Z*(0.763297E-03+   &
        Z*0.97087007E-05))))+0.5772156649E+00*(0.99999207E+00+                    &
        Z*(-0.149545886E+01+Z*(0.41806426E-01+Z*(-0.3000591E-01+                  &
        Z*(0.19387339E-02+Z*(-0.51801555E-03)))))))/(0.1E+01+                     &
        Z*(-0.76273617E+00+Z*(0.28388363E+00+Z*(-0.66786033E-01+Z*(0.12982719E-01 &
        +Z*(-0.8700861E-03+Z*0.2989204E-03))))))
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
          + (-3.777369, -4.510900))
      ELSE
        GG = ((((((( (1.000000, -1.3935496E-06)*Z+ (15.82958, 20.14222))  &
          *Z+ (-70.52863, 227.9511))*Z+ (-985.4221, 226.6272))*Z          &
          + (-1202.318, -1580.907))*Z+ (953.2441, -1342.447))*Z           &
          + (417.3716, 196.6665))*Z+ (-9.881266, 24.24952))/              &
          (((((((( (1.000000, 0.0000000E+00)*Z+ (16.83184, 20.14481))*Z   &
          + (-55.66969, 248.1167))*Z+ (-1068.640, 434.4707))*Z            &
          + (-2082.250, -1522.471))*Z+ (383.3455, -2730.378))*Z           &
          + (1216.791, -351.7189))*Z+ (115.3926, 161.2647))*Z             &
          + (-3.777369, 4.510900))
      END IF
    END IF
    IF (AIMAG(Z) < 0) THEN
      GG = GG - II*PI*CEX
    ELSE
      GG = GG + II*PI*CEX
    END IF
  END FUNCTION

!---------------------------------------------------------------------

  SUBROUTINE INITIALIZE_GREEN(XR, XZ, APD)
    ! Initialize XR, XZ and APD
    ! Those parameters are independent of the depth and the frequency.
    ! Thus, they are initialized only once at the beginning of the execution of the code.

    ! References:
    ! [1] Delhommeau, Amélioration des codes de calcul de diffraction-radiation, 2èmes journées de l'hydrodynamique, 1989
    ! [2] Babarit and Delhommeau, Theoretical and numerical aspects of the open source BEM solver NEMOH, EWTEC 2015

    ! Output
    REAL, DIMENSION(328), INTENT(OUT)     :: XR
    REAL, DIMENSION(46),  INTENT(OUT)     :: XZ
    REAL, DIMENSION(328, 46, 2, 2), INTENT(OUT) :: APD

    ! Local variables
    INTEGER :: I, J, K
    REAL :: THETA(NPINTE), CQT(NPINTE)
    REAL :: CT
    COMPLEX :: C1, C2, ZETA, CEX

    ! Initialize XZ (named Z(J) in [1, 2])
    DO J = 1, JZ
      XZ(J) = -AMIN1(10**(J/5.0-6), 10**(J/8.0-4.5), 16.)
    END DO

    ! Initialize XR (named X(I) in [1, 2])
    XR(1) = 0.0
    DO I = 2, IR
      IF (I < 40) THEN
        XR(I) = AMIN1(10**((I-1.0)/5-6), 4.0/3.0 + ABS(I-32)/3.0)
      ELSE
        XR(I) = 4.0/3.0 + ABS(I-32)/3.0
      ENDIF
    END DO

    ! Initialize THETA and CQT for the integration between -pi/2 and pi/2 with the Simpson rule.
    DO K = 1, NPINTE
      THETA(K) = -PI/2 + (K-1)/(NPINTE-1.0)*PI
      IF ((K == 1) .OR. (K == NPINTE)) THEN
        CQT(K) = PI/(3*(NPINTE-1))
      ELSEIF (MOD(K,2)==0) THEN
        CQT(K) = 4.0/(3*(NPINTE-1))*PI
      ELSE
        CQT(K) = 2.0/(3*(NPINTE-1))*PI
      ENDIF
    ENDDO

    ! Initialize APD..
    APD(:, :, :, :) = 0.0
    DO J = 1, JZ
      DO I = 1, IR
        DO K = 1, NPINTE
          CT = COS(THETA(K))
          ZETA = CMPLX(XZ(J), XR(I)*CT)
          IF (REAL(ZETA) <= -30.0) THEN
            CEX = (0.0, 0.0)
          ELSE
            CEX = CEXP(ZETA)
          ENDIF
          C1 = CQT(K)*(GG(ZETA, CEX) - 1.0/ZETA)
          C2 = CQT(K)*CEX
          APD(I, J, 1, 1) = APD(I, J, 1, 1) + CT*AIMAG(C1) ! named D_1(Z, X) in [1, 2]
          APD(I, J, 1, 2) = APD(I, J, 1, 2) + REAL(C1)     ! named D_2(Z, X) in [1, 2]
          APD(I, J, 2, 1) = APD(I, J, 2, 1) + CT*AIMAG(C2) ! named Z_1(Z, X) in [1, 2]
          APD(I, J, 2, 2) = APD(I, J, 2, 2) + REAL(C2)     ! named Z_2(Z, X) in [1, 2]
        END DO
      END DO
    END DO

    RETURN

  END SUBROUTINE INITIALIZE_GREEN

!-------------------------------------------------------------------------------!

  FUNCTION FF(XTT, AK, AM)

    ! Input
    REAL, INTENT(IN) :: XTT, AK, AM

    ! Local variables
    REAL :: COEF, TOL, A, B, C, D, E, F, FF

    COEF = (AM+AK)**2/(AM**2-AK**2+AK)

    TOL = MAX(0.1, 0.1*AM)
    IF (ABS(XTT-AM) > TOL) THEN
      FF = (XTT+AK)*EXP(XTT)/(XTT*SINH(XTT)-AK*COSH(XTT)) - COEF/(XTT-AM) - 2
    ELSE
      A = AM - TOL
      B = AM
      C = AM + TOL
      D = (A+AK)*EXP(A)/(A*SINH(A)-AK*COSH(A)) - COEF/(A-AM) - 2
      E = COEF/(AM+AK)*(AM+AK+1)               - (COEF/(AM+AK))**2*AM - 2
      F = (C+AK)*EXP(C)/(C*SINH(C)-AK*COSH(C)) - COEF/(C-AM) - 2
      FF = D*(XTT-B)*(XTT-C)/((A-B)*(A-C)) + &
           E*(XTT-C)*(XTT-A)/((B-C)*(B-A)) + &
           F*(XTT-A)*(XTT-B)/((C-A)*(C-B))
    ENDIF
    
    RETURN
  END FUNCTION FF

!-------------------------------------------------------------------------------!

! The code below perform the approximation of the function FF above as a sum of
! exponentials. They have been partially commented and refactored. A more
! readable and maintainable Python version has been implemented in Nemoh.py.
! However, as of today (January 2018), the python version is ~100 times slower
! than the Fortran version below. Thus, the code below has been preserved. It
! could become useful if ever this part of the code would become time critical.

  SUBROUTINE LISC(AK0, wavenumber, &
                  AMBDA, AR, NEXP)
    ! Compute AMBDA and AR and test their values

    INTEGER, PARAMETER :: NEXR = 31
    INTEGER, PARAMETER :: NMAX = 4*(31-1)

    ! Range on which the function FF is tabulated
    REAL, PARAMETER :: A = -0.1
    REAL, PARAMETER :: B = 20.0

    REAL, PARAMETER :: PRECI = 1.0e-2

    ! Inputs
    REAL, INTENT(IN)                 :: AK0, wavenumber

    ! Outputs
    INTEGER, INTENT(OUT)             :: NEXP
    REAL, DIMENSION(31), INTENT(OUT) :: AMBDA, AR

    ! Local variables
    LOGICAL :: ISOR
    INTEGER :: I, J 
    INTEGER :: NK, NM, NMO
    REAL :: XX, YY, TT
    REAL :: XT(4*(31-1)+1), YT(4*(31-1)+1) ! Tabulation of FF

    ! Initialize variables to be computed
    AMBDA = 0.0
    AR = 0.0

    ! Number of points in the tabulation (first try)
    NK = 4

62 CONTINUE

    ! Tabulate the function FF
    DO I = 1, (4*NK)+1
      XT(I) = A + (B-A)*(I-1)/(4*NK)
      YT(I) = FF(XT(I), AK0, wavenumber)
    END DO

    ! Compute AMBDA and AR based on this tabulation
    NM = NK
    CALL EXPORS(XT, YT, NM, NMAX, AMBDA, AR)

    ! Test values of AMBDA and AR
    ISOR = .False.
    NMO = 8*NM+1
    DO I = 1, NMO
      XX = (I-1)*B/(NMO-1)

      ! Exact value
      YY = FF(XX, AK0, wavenumber)

      ! Compute the sum of exponentials
      TT = 0.0
      DO J = 1, NM
        TT = TT + AR(J)*EXP(AMBDA(J)*XX)
      END DO

      ! Compare the error to a threshold
      IF (ABS(YY-TT) > PRECI) ISOR = .True.
    END DO

    ! The error was higher than the wanted precision
    IF (ISOR .AND. (NK <= NEXR-2)) THEN
      ! Add more exponentials to the sum and restart
      NK = NK+2
      GOTO 62
    END IF

    ! Final number of exponentials in the sum
    NEXP = NM

    RETURN

  END SUBROUTINE LISC

!---------------------------------------------------------------------

  SUBROUTINE EXPORS(XT, YT, NM, NMAX, AMBDA, AR)
    ! Compute AMBDA and AR based on XT and YT using Prony's method.

    INTEGER, INTENT(IN) :: NMAX
    REAL, INTENT(IN)    :: XT(4*(31-1)+1), YT(4*(31-1)+1)

    INTEGER, INTENT(INOUT) :: NM ! Number of exponentials in the sum.
    ! The number of exponentials will change during the computation, for
    ! instance if some roots of the polynomial are double and thus merged.

    REAL, INTENT(OUT) :: AR(31), AMBDA(31)

    ! Local variables
    INTEGER :: NPP

    REAL    :: SC(31), VR(31), VC(31)
    REAL    :: S(4*(31-1), 31+1)
    COMPLEX :: COM(31)

    INTEGER::I,J,K,JJ,II,IJ,MN,NEXP
    REAL::H,EPS

    ! Assemble an over-determined linear system
    NPP = 4*NM+1
    H = (XT(NPP) - XT(1))/(4*NM)
    K = NPP-NM
    DO I = 1, K
      DO J = 1, NM
        JJ = NM - J + I
        S(I, J) = YT(JJ)
      END DO
      II = NM + I
      S(I, NM+1) = -YT(II)
    END DO

    ! Solve the system
    CALL HOUSRS(S, NMAX, K, NM)

    ! Assemble the coefficients of the polynomial
    DO I = 1, NM
      IJ = NM-I+1
      SC(IJ) = S(I, NM+1)
    END DO
    MN = NM+1
    SC(MN) = 1.0

    ! Get roots of the polynomial
    CALL SPRBM(SC, MN, VR, VC)

    COM(1:NM) = CMPLX(VR(1:NM), VC(1:NM))
    COM(1:NM) = CLOG(COM(1:NM))/H
    VR(1:NM) = REAL(COM(1:NM))
    VC(1:NM) = AIMAG(COM(1:NM))

    ! Keep only the real part of the roots.
    ! In case of a double root, keep only one.
      I=1
      J=0
  100 IF(VC(I))110,111,110
  ! 100 IF (VC(I) == 0) THEN
  111 J=J+1
      AMBDA(J)=VR(I)
      I=I+1
      GO TO 101
    ! ELSE
  110 IF(ABS(VR(I)-VR(I+1))-1.E-5)120,120,121
      ! IF (ABS(VR(I)-VR(I+1)) <= -1.E-5) THEN
  120 J=J+1
      AMBDA(J)=VR(I)
      I=I+2
      GO TO 101
    ! ELSE
  121 J=J+1
      AMBDA(J)=VR(I)
      I=I+1
    ! END IF
    ! END IF
  101 IF(I-NM)100,100,102
  102 NEXP=J

    ! Rewrite AMBDA while keeping only the values between -20.0 and 0.0
    J = 0
    DO I = 1, NEXP
      IF ((AMBDA(I) > -20.0) .AND. (AMBDA(I) <= 0.0)) THEN
        J = J+1
        AMBDA(J) = AMBDA(I)
      END IF
    END DO

    ! Number of relevant values in AMBDA
    NEXP = J
    NM = NEXP

    ! Compute AR
    CALL MCAS(AMBDA, XT, YT, NPP, AR, S, NMAX, NEXP)

    RETURN
  END SUBROUTINE EXPORS

!----------------------------------------------------------------------------

  SUBROUTINE MCAS(AMBDA, XT, YT, NPP, AR, A, NMAX, NEXP)
    ! Compute AR as the least-square solution of an over-determined system.

    REAL, INTENT(IN)    :: AMBDA(31), XT(4*(31-1)+1), YT(4*(31-1)+1)
    INTEGER, INTENT(IN) :: NPP, NMAX, NEXP
    REAL, INTENT(INOUT) :: A(4*(31-1), 31+1)
    REAL, INTENT(OUT)   :: AR(31)

    ! Local variables
    INTEGER :: I, J, L
    REAL :: S, TT, TTT

    ! Assemble the system
    DO I = 1, NEXP
      DO J = 1, NEXP
        S = 0
        DO L = 1, NPP
          TT = (AMBDA(I)+AMBDA(J))*XT(L)
          IF(TT >= -30) THEN
            S = S + EXP(TT)
          END IF
        END DO
        A(I,J) = S
      END DO
    END DO

    DO I = 1, NEXP
      S = 0
      DO L = 1, NPP
        TTT = AMBDA(I)*XT(L)
        IF(TTT >= -30) THEN
          S = S + EXP(TTT)*YT(L)
        END IF
      END DO
      A(I, NEXP+1) = S
    END DO

    ! Solve it
    CALL HOUSRS(A, NMAX, NEXP, NEXP)

    ! Extract results
    DO I = 1, NEXP
      AR(I) = A(I, NEXP+1)
    END DO

    RETURN
  END SUBROUTINE MCAS

!---------------------------------------------------------------------

  SUBROUTINE HOUSRS(A, NMAX, NL, NCC)
    ! Least-square solver for over-determined system.
    ! NCC : Number of columns
    ! NL: Number of lines

    INTEGER, INTENT(IN) :: NMAX, NL, NCC
    REAL, INTENT(INOUT) :: A(NMAX, 31+1)

    ! Local variables
    INTEGER, PARAMETER :: NS = 1
    REAL, PARAMETER :: EPS = 1e-20

    INTEGER :: I, J, K, L
    REAL    :: E, E0, AR, BA, ETA

    IF (NCC > NL) THEN
      PRINT*, ('Error in HOUSRS: number of columns > number of lines')
      STOP
    ENDIF

    DO K = 1, NCC
      E = 0.0
      DO I = K, NL
        E = E + A(I, K)**2
      END DO
      
      E0 = SQRT(E)
      ! IF(E0.LT.EPS)THEN
      !   WRITE(*, 201)EPS
    ! 201 FORMAT(1X, 'NORME INFERIEURE A ', 1PE16.6/)
      !   STOP
      ! ENDIF
      IF (A(K, K) == 0) THEN
        AR = -E0
      ELSE
        AR = -SIGN(1.0, A(K, K)) * E0
      ENDIF

      ETA = AR*(AR - A(K, K))
      DO J = K+1, NCC+NS
        BA = (A(K, K)-AR)*A(K, J)
        DO I = K+1, NL
          BA = BA+A(I, K)*A(I, J)
        END DO
        A(K, J) = A(K, J)+BA/AR
        DO I = K+1, NL
          A(I, J) = A(I, J)-A(I, K)*BA/ETA
        END DO
      END DO

      A(K, K) = AR
      DO I = K+1, NL
        A(I, K) = 0
      END DO
    END DO

    DO J = 1, NS
      A(NCC, NCC+J) = A(NCC, NCC+J)/A(NCC, NCC)
      DO L = 2, NCC
        DO I = NCC+2-L, NCC
          A(NCC+1-L, NCC+J) = A(NCC+1-L, NCC+J) - A(NCC+1-L, I)*A(I, NCC+J)
        END DO
        A(NCC+1-L, NCC+J) = A(NCC+1-L, NCC+J)/A(NCC+1-L, NCC+1-L)
      END DO
    END DO

    RETURN
  END SUBROUTINE HOUSRS

!----------------------------------------------------------------------

  SUBROUTINE SPRBM(C,IC,RR,RC)
    ! Somehow get the roots of a polynomial.

    INTEGER::IC,IR,IER
    REAL :: C(31), RR(31), RC(31)
    REAL :: POL(31)
    INTEGER::I,J,L,N,LIM,IST
    REAL::A,B,H
    REAL::EPS,Q1,Q2,Q(4)

      EPS=1.E-6
      LIM=100
      IR=IC+1
    1 IR=IR-1
      IF(IR-1)42,42,2
    2 IF(C(IR))3,1,3
    3 IER=0
      J=IR
      L=0
      A=C(IR)
      DO 8 I=1,IR
      IF(L)4,4,7
    4 IF(C(I))6,5,6
    5 RR(I)=0.
      RC(I)=0.
      POL(J)=0.
      J=J-1
      GO TO 8
    6 L=1
      IST=I
      J=0
    7 J=J+1
      C(I)=C(I)/A
      POL(J)=C(I)
      IF(ABS(POL(J))-1.E27)8,42,42
    8 CONTINUE
      Q1=0.
      Q2=0.
    9 IF(J-2)33,10,14
   10 A=POL(1)
      RR(IST)=-A
      RC(IST)=0.
      IR=IR-1
      Q2=0.
      IF(IR-1)13,13,11
   11 DO 12 I=2,IR
      Q1=Q2
      Q2=POL(I+1)
   12 POL(I)=A*Q2+Q1
   13 POL(IR+1)=A+Q2
      GO TO 34
   14 DO 22 L=1,10
      N=1
   15 Q(1)=Q1
      Q(2)=Q2
      CALL SPQFB(POL,J,Q,LIM,I)
      IF(I)16,24,23
   16 IF(Q1)18,17,18
   17 IF(Q2)18,21,18
   18 GOTO(19,20,19,21),N
   19 Q1=-Q1
      N=N+1
      GO TO 15
   20 Q2=-Q2
      N=N+1
      GO TO 15
   21 Q1=1.+Q1
   22 Q2=1.-Q2
      IER=3
      IR=IR-J
      GOTO 45
   23 IER=1
   24 Q1=Q(1)
      Q2=Q(2)
      B=0.
      A=0.
      I=J
   25 H=-Q1*B-Q2*A+POL(I)
      POL(I)=B
      B=A
      A=H
      I=I-1
      IF(I-2)26,26,25
   26 POL(2)=B
      POL(1)=A
      L=IR-1
      IF(J-L)27,27,29
   27 DO 28 I=J,L
   28 POL(I-1)=POL(I-1)+POL(I)*Q2+POL(I+1)*Q1
   29 POL(L)=POL(L)+POL(L+1)*Q2+Q1
      POL(IR)=POL(IR)+Q2
      H=-.5*Q2
      A=H*H-Q1
      B=SQRT(ABS(A))
      IF(A)30,30,31
   30 RR(IST)=H
      RC(IST)=B
      IST=IST+1
      RR(IST)=H
      RC(IST)=-B
      GO TO 32
   31 B=H+SIGN(B,H)
      RR(IST)=Q1/B
      RC(IST)=0.
      IST=IST+1
      RR(IST)=B
      RC(IST)=0.
   32 IST=IST+1
      J=J-2
      GO TO 9
   33 IR=IR-1
   34 A=0.
      DO 38 I=1,IR
      Q1=C(I)
      Q2=POL(I+1)
      POL(I)=Q2
      IF(Q1)35,36,35
   35 Q2=(Q1-Q2)/Q1
   36 Q2=ABS(Q2)
      IF(Q2-A)38,38,37
   37 A=Q2
   38 CONTINUE
      I=IR+1
      POL(I)=1.
      RR(I)=A
      RC(I)=0.
      IF(IER)39,39,41
   39 IF(A-EPS)41,41,40
   40 IER=-1
   41 GOTO 45
   42 IER=2
      IR=0
   45 IF(IER-2)46,47,46
   47 WRITE(*,48)IER
   48 FORMAT(/5X,'IER = ',I3,'  ERREUR DANS SPRBM'/)
      STOP
   46 RETURN

  END SUBROUTINE SPRBM

!----------------------------------------------------------------

  SUBROUTINE SPQFB(C,IC,Q,LIM,IER)

    INTEGER::IC,LIM,IER,I,J,L,LL
    REAL:: C(31),Q(4)
    REAL:: H,HH,A,A1,AA,B,BB,B1,C1,CA,CB,CC,CD,DQ1,DQ2,EPS,EPS1
    REAL:: Q1,Q2,QQ1,QQ2,QQQ1,QQQ2

!--------- Value non initialized in previous versions ?!
      H=0.
      HH=0.
!---------
      IER=0
      J=IC+1
    1 J=J-1
      IF(J-1) 40,40,2
    2 IF(C(J)) 3,1,3
    3 A=C(J)
      IF(A-1.) 4,6,4
    4 DO 5 I=1,J
      C(I)=C(I)/A
      IF(ABS(C(I))-1.E27)5,40,40
    5 CONTINUE
    6 IF(J-3) 41,38,7
    7 EPS=1.E-14
      EPS1=1.E-6
      L=0
      LL=0
      Q1=Q(1)
      Q2=Q(2)
      QQ1=0.
      QQ2=0.
      AA=C(1)
      BB=C(2)
      CB=ABS(AA)
      CA=ABS(BB)
      IF(CB-CA) 8,9,10
    8 CC=CB+CB
      CB=CB/CA
      CA=1.
      GO TO 11
    9 CC=CA+CA
      CA=1.
      CB=1.
      GO TO 11
   10 CC=CA+CA
      CA=CA/CB
      CB=1.
   11 CD=CC*.1
   12 A=0.
      B=A
      A1=A
      B1=A
      I=J
      QQQ1=Q1
      QQQ2=Q2
      DQ1=HH
      DQ2=H
   13 H=-Q1*B-Q2*A+C(I)
      IF(ABS(H)-1.E27)14,42,42
   14 B=A
      A=H
      I=I-1
      IF(I-1) 18,15,16
   15 H=0.
   16 H=-Q1*B1-Q2*A1+H
      IF(ABS(H)-1.E27)17,42,42
   17 C1=B1
      B1=A1
      A1=H
      GO TO 13
   18 H=CA*ABS(A)+CB*ABS(B)
      IF(LL) 19,19,39
   19 L=L+1
      IF(ABS(A)-EPS*ABS(C(1))) 20,20,21
   20 IF(ABS(B)-EPS*ABS(C(2))) 39,39,21
   21 IF(H-CC) 22,22,23
   22 AA=A
      BB=B
      CC=H
      QQ1=Q1
      QQ2=Q2
   23 IF(L-LIM) 28,28,24
   24 IF(H-CD) 43,43,25
   25 IF(Q(1)) 27,26,27
   26 IF(Q(2)) 27,42,27
   27 Q(1)=0.
      Q(2)=0.
      GO TO 7
   28 HH=AMAX1(ABS(A1),ABS(B1),ABS(C1))
      IF(HH) 42,42,29
   29 A1=A1/HH
      B1=B1/HH
      C1=C1/HH ! Has C1 been intialized?
      H=A1*C1-B1*B1
      IF(H) 30,42,30
   30 A=A/HH
      B=B/HH
      HH=(B*A1-A*B1)/H
      H=(A*C1-B*B1)/H
      Q1=Q1+HH
      Q2=Q2+H
      IF(ABS(HH)-EPS*ABS(Q1)) 31,31,33
   31 IF(ABS(H)-EPS*ABS(Q2)) 32,32,33
   32 LL=1
      GO TO 12
   33 IF(L-1)12,12,34
   34 IF(ABS(HH)-EPS1*ABS(Q1)) 35,35,12
   35 IF(ABS(H)-EPS1*ABS(Q2)) 36,36,12
   36 IF(ABS(QQQ1*HH)-ABS(Q1*DQ1)) 37,44,44
   37 IF(ABS(QQQ2*H)-ABS(Q2*DQ2)) 12,44,44
   38 Q(1)=C(1)
      Q(2)=C(2)
      Q(3)=0.
      Q(4)=0.
      GOTO 45
   39 Q(1)=Q1
      Q(2)=Q2
      Q(3)=A
      Q(4)=B
      GOTO 45
   40 IER=-1
      GOTO 45
   41 IER=-2
      GOTO 45
   42 IER=-3
      GO TO 44
   43 IER=1
   44 Q(1)=QQ1
      Q(2)=QQ2
      Q(3)=AA
      Q(4)=BB
   45 RETURN

  END SUBROUTINE SPQFB

!----------------------------------------------------------------

END MODULE Initialize_Green_2
