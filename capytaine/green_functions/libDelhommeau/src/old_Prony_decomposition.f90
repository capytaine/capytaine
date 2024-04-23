! Copyright (C) 2017-2024 Matthieu Ancellin
! See LICENSE file at <https://github.com/capytaine/libDelhommeau>

MODULE Old_Prony_decomposition

  USE FLOATING_POINT_PRECISION, ONLY: PRE
  USE CONSTANTS

  IMPLICIT NONE

! The code below performs the approximation of the function FF as a sum of exponentials.
! Unlike the rest of the Fortran code, this module has not been completely refactored.
! A more readable and maintainable Python version has been implemented in libDelhommeau/bem/prony_decomposition.py.
! However, as of today (March 2019), the two versions do not give exactly the same result in all cases.
! For compatibility with Nemoh, the old version has thus been preserved.

! Unlike most of the rest of the Fortran code, these functions are only in single precision.

  PUBLIC :: FF
  PUBLIC :: LISC   ! Initialization of AMBDA and AR
  PUBLIC :: EXPORS ! Called by LISC
  PUBLIC :: MCAS   ! Called by EXPORS
  PUBLIC :: HOUSRS ! Called by EXPORS and MCAS
  PUBLIC :: SPRBM  ! Called by EXPORS
  PUBLIC :: SPQFB  ! Called by SPRBM

CONTAINS

  FUNCTION FF(K, K0, M0)
    ! Function F_2(K) = F_1(K) - 2 from (3.47) of Delhommeau's thesis

    ! Input
    REAL, INTENT(IN) :: K, K0, M0

    ! Local variables
    REAL :: COEF, TOL, A, B, C, D, E, F, FF

    COEF = (M0+K0)**2/(M0**2-K0**2+K0)

    TOL = MAX(0.1, 0.1*M0)
    IF (ABS(K-M0) > TOL) THEN
      FF = (K+K0)*EXP(K)/(K*SINH(K)-K0*COSH(K)) - COEF/(K-M0) - 2
    ELSE  ! The fraction above is undefined 0/0 when K == M0
      A = M0 - TOL
      B = M0
      C = M0 + TOL
      D = (A+K0)*EXP(A)/(A*SINH(A)-K0*COSH(A)) - COEF/(A-M0) - 2
      E = COEF/(M0+K0)*(M0+K0+1)               - (COEF/(M0+K0))**2*M0 - 2
      F = (C+K0)*EXP(C)/(C*SINH(C)-K0*COSH(C)) - COEF/(C-M0) - 2
      FF = D*(K-B)*(K-C)/((A-B)*(A-C)) + &
           E*(K-C)*(K-A)/((B-C)*(B-A)) + &
           F*(K-A)*(K-B)/((C-A)*(C-B))
    ENDIF
  END FUNCTION FF


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
      C1=C1/HH ! Has C1 been initialized?
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

END MODULE Old_Prony_decomposition
