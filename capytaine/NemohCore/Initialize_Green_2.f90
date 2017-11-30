MODULE Initialize_Green_2

  IMPLICIT NONE

  PUBLIC :: INITIALIZE_GREEN
  PUBLIC :: GG

  PUBLIC  :: LISC   ! Initialization of AMBDA and AR
  PUBLIC :: EXPORS ! Called by LISC
  PUBLIC :: MCAS   ! Called by EXPORS
  PUBLIC :: SPRBM  ! Called by EXPORS
  PUBLIC :: SPQFB  ! Called by SPRBM
  PUBLIC :: HOUSRS ! Called by EXPORS and MCAS
  PUBLIC :: FF     ! Called by LISC

  REAL, PARAMETER :: PI = 3.141592653588979 ! π
  COMPLEX, PARAMETER :: II = (0, 1)         ! Imaginary unit

  ! Independent of Omega
  INTEGER, PARAMETER      :: NPINTE = 251
  INTEGER, PARAMETER      :: IR = 328
  INTEGER, PARAMETER      :: JZ = 46
  REAL, DIMENSION(IR)     :: XR
  REAL, DIMENSION(JZ)     :: XZ
  REAL, DIMENSION(IR, JZ) :: APD1X, APD1Z, APD2X, APD2Z

  ! Dependant of Omega
  INTEGER                 :: NEXP
  REAL, DIMENSION(31)     :: AMBDA, AR

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

!-------------------------------------------------------------------------------!

! The subroutine below initialize AMBDA and AR which are dependent on the wave number.
! This part of the code is still in old-fashionned style.
! TODO: clean that up.

  SUBROUTINE LISC(AK0,wavenumber)
    ! Compute AMBDA and AR

    INTEGER :: I,J,NJ,NPP, NM
    REAL:: AK0,wavenumber
    REAL:: POL(31),A,B,depth
    REAL:: S(4*(31-1),31+1),XT(4*(31-1)+1),YT(4*(31-1)+1)
    REAL:: SC(31),VR(31),VC(31)
    INTEGER::ISTIR,NMAX,NK,ISOR,NPI,NMO,NEXR
    REAL:: PRECI,ERMAX,ERMOY,XX,YY,TT,DIF,RT
    COMPLEX:: COM(31)

    AMBDA=0.0

    S=0.
    SC=0.
    AR=0.
    NEXR=31
    PRECI=1.E-02
    ISTIR=0
    NMAX=4*(NEXR-1)
    NK=4
    A=-0.1
    B=20.
62 CONTINUE
    NM=NK
    NJ=4*NM
    NPP=NJ+1
    depth=(B-A)/NJ
    DO I=1,NPP
      XT(I)=A+(I-1)*depth
      YT(I)=FF(XT(I),AK0,wavenumber)
    END DO
    ISOR=0
    CALL EXPORS(XT,YT,NJ,NM,AMBDA,NMAX,S,SC,VR,VC,COM,POL,AR)

    NPI=2
    NMO=NPI*NPP-NPI+1
    ERMAX=0.
    ERMOY=0.
    DO I=1,NMO
      XX=(I-1)*B/(NMO-1)
      YY=FF(XX,AK0,wavenumber)
      TT=0.
      DO J=1,NM
        RT=AMBDA(J)*XX
        IF(RT.GT.-20.)THEN
          TT=TT+AR(J)*EXP(RT)
        ENDIF
      END DO
      DIF=YY-TT
      ERMOY=ERMOY+DIF
      ERMAX=AMAX1(ERMAX,ABS(DIF))
      IF(ABS(DIF).GT.PRECI) ISOR=1
    END DO
    ERMOY=ERMOY/NMO

    ! WRITE(*,1111) NM,ERMAX,ERMOY
    ! 1111 FORMAT(5X,I2,'EXPONENTIELLES  ECART MAXI = ',E10.3,'ECART MOYEN = ',E10.3/)

    IF ((ISTIR .NE. 1) .AND. (ISOR .NE. 0)) THEN
      NK=NK+2
      IF (NK-(NEXR-1)>0) THEN
        ! WRITE(*,6500) PRECI,NM
        ! 6500 FORMAT(/5X,'PRECISION = ',E10.3,'  NON ATTEINTE AVEC ',I2,'  EXPONENTIELLES')
        ! STOP
      ELSE
        GOTO 62
      ENDIF

    ELSE
      DO J=1,NM
        ! WRITE(*,1100) AR(J),AMBDA(J)
        ! 1100 FORMAT(5X,E16.7,'EXP(',E16.7,')')
        IF (AMBDA(J).GT.0.) STOP
      END DO
    END IF

    NEXP=NM

    RETURN

  END SUBROUTINE LISC

!-------------------------------------------------------------------------------!

  SUBROUTINE EXPORS(XT,YT,NJ,NM,VCOM,NMAX,S,SC,VR,VC,COM,POL,AR)

    INTEGER::NJ,NM,NMAX
    REAL:: VCOM(31),POL(31),AR(31),SC(31),VR(31),VC(31)
    REAL:: S(4*(31-1),31+1),XT(4*(31-1)+1),YT(4*(31-1)+1)
    COMPLEX:: COM(31)
    INTEGER::I,J,K,NPP,JJ,II,IJ,MN
    INTEGER::IS,IER
    REAL::H,EPS

    NPP=NJ+1
    H=(XT(NPP)-XT(1))/NJ
    K=NPP-NM
    DO I=1,K
      DO J=1,NM
        JJ=NM-J+I
        S(I,J)=YT(JJ)
      END DO
      II=NM+I
      S(I,NM+1)=-YT(II)
    END DO
    EPS=1.E-20

    CALL HOUSRS(S,NMAX,K,NM,1,EPS)
    DO I=1,NM
      IJ=NM-I+1
      SC(IJ)=S(I,NM+1)
    END DO
    MN=NM+1
    SC(MN)=1.
    CALL SPRBM(SC,MN,VR,VC,POL,IS,IER)
    DO I=1,NM
      COM(I)=CMPLX(VR(I),VC(I))
      COM(I)=CLOG(COM(I))/H
      VR(I)=REAL(COM(I))
      VC(I)=AIMAG(COM(I))
    END DO

      I=1
      J=0
  100 IF(VC(I))110,111,110
  111 J=J+1
      VCOM(J)=VR(I)
      I=I+1
      GO TO 101
  110 IF(ABS(VR(I)-VR(I+1))-1.E-5)120,120,121
  120 J=J+1
      VCOM(J)=VR(I)
      I=I+2
      GO TO 101
  121 J=J+1
      VCOM(J)=VR(I)
      I=I+1
  101 IF(I-NM)100,100,102
  102 NEXP=J
      J=0
      DO 300 I=1,NEXP
      J=J+1
      IF(VCOM(I).GE.0.)GOTO 301
      IF(VCOM(I)+20.)301,301,302
  301 J=J-1
      GO TO 300
  302 VCOM(J)=VCOM(I)
  300 CONTINUE
      NEXP=J
      NM=NEXP
      CALL MCAS(VCOM,XT,YT,NPP,AR,S,NMAX)

    RETURN
  END SUBROUTINE EXPORS
!----------------------------------------------------------------------------

  SUBROUTINE MCAS(TEXP,XT,YT,NPP,AR,A,NMAX)

    INTEGER:: NPP,NMAX
    REAL::XT(4*(31-1)+1),YT(4*(31-1)+1),A(4*(31-1),31+1),AR(31),TEXP(31)
    INTEGER::I,J,L,M,N
    REAL::S,TT,TTT,EPS

      EPS=1.E-20
      DO 1 I=1,NEXP
      DO 1 J=1,NEXP
      S=0
      DO 3 L=1,NPP
      TT=(TEXP(I)+TEXP(J))*XT(L)
      IF(TT+30)3,4,4
    4 S=S+EXP(TT)
    3 CONTINUE
      A(I,J)=S
    1 CONTINUE
      DO 5 I=1,NEXP
      S=0
      DO 6 L=1,NPP
      TTT=TEXP(I)*XT(L)
      IF(TTT+30)6,7,7
    7 S=S+EXP(TTT)*YT(L)
    6 CONTINUE
      A(I,NEXP+1)=S
    5 CONTINUE
      N=NEXP
      M=N+1
      CALL HOUSRS(A,NMAX,N,N,1,EPS)
      DO 10 I=1,NEXP
   10 AR(I)=A(I,NEXP+1)
      RETURN

  END SUBROUTINE MCAS

!----------------------------------------------------------------------

  SUBROUTINE SPRBM(C,IC,RR,RC,POL,IR,IER)

    INTEGER::IC,IR,IER
    REAL:: C(31),RR(31),RC(31),POL(31)
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
!---------------------------------------------------------------------

  SUBROUTINE HOUSRS(A,NMAX,NL,NCC,NS,EPS)

    INTEGER::NMAX,NL,NCC,NS
    REAL:: A(NMAX,31+1),EPS
    INTEGER::I,J,K,L,M,NCJ,NTC,KP1
    REAL::E,E0,AR,BA,ETA
    INTEGER :: I1

      NTC=NCC+NS
      IF(NCC.GT.NL)THEN
        WRITE(*,3010)
  3010 FORMAT(' NBRE DE COLONNES > NBRES DE LIGNES')
        STOP
      ENDIF
      DO 13 K=1,NCC
        E=0
        DO 1101 I=K,NL
          E=E+A(I,K)**2
   1101 CONTINUE
        E0=SQRT(E)
        IF(E0.LT.EPS)THEN
          WRITE(*,201)EPS
      201 FORMAT(1X,'NORME INFERIEURE A ',1PE16.6/)
          STOP
        ENDIF
        IF(A(K,K).EQ.0)THEN
          AR=-E0
        ELSE
          AR=-SIGN(E0,A(K,K))
        ENDIF
        ETA=AR*(AR-A(K,K))
        KP1=K+1
        DO 10 J=KP1,NTC
          BA=(A(K,K)-AR)*A(K,J)
          DO 9 I=KP1,NL
            BA=BA+A(I,K)*A(I,J)
        9 CONTINUE
          A(K,J)=A(K,J)+BA/AR
          DO 11 I=KP1,NL
            A(I,J)=A(I,J)-A(I,K)*BA/ETA
       11 CONTINUE
     10 CONTINUE
        A(K,K)=AR
        DO 12 I=KP1,NL
  12   A(I,K)=0
 13   CONTINUE
      DO 1006 J=1,NS
        NCJ=NCC+J
        A(NCC,NCJ)=A(NCC,NCJ)/A(NCC,NCC)
        DO 1005 L=2,NCC
          I1=NCC+1-L
          M=I1+1
          DO 1004 I=M,NCC
      1004 A(I1,NCJ)=A(I1,NCJ)-A(I1,I)*A(I,NCJ)
    1005 A(I1,NCJ)=A(I1,NCJ)/A(I1,I1)
 1006 CONTINUE
      RETURN

  END SUBROUTINE HOUSRS

!---------------------------------------------------------------------

  FUNCTION FF(XTT,AK,AM)

    REAL::XTT,AK,AM
    REAL::COEF,TOL,H,A,B,C,D,E,F,FF

    COEF=(AM+AK)**2/(AM**2-AK**2+AK)
    H=XTT-AM
    TOL=AMAX1(0.1,0.1*AM)
    IF(ABS(H).GT.TOL)THEN
      FF=(XTT+AK)*EXP(XTT)/(XTT*SINH(XTT)-AK*COSH(XTT))-COEF/(XTT-AM)-2.
    ELSE
      A=AM-TOL
      B=AM
      C=AM+TOL
      D=(A+AK)*EXP(A)/(A*SINH(A)-AK*COSH(A))-COEF/(A-AM)-2
      E=COEF/(AM+AK)*(AM+AK+1)-(COEF/(AM+AK))**2*AM-2
      F=(C+AK)*EXP(C)/(C*SINH(C)-AK*COSH(C))-COEF/(C-AM)-2
      FF=(XTT-B)*(XTT-C)*D/((A-B)*(A-C))+(XTT-C)*(XTT-A)*E/((B-C)*(B-A))+&
        & (XTT-A)*(XTT-B)*F/((C-A)*(C-B))
    ENDIF
    RETURN

  END FUNCTION FF
END MODULE Initialize_Green_2
