module FinGrnExtSubs_module

   implicit none

contains

!!---------------------------------------------------------------------------!!
!        Belows are  third-party (external) Level-3 (low level) subroutines      !
!!---------------------------------------------------------------------------!!
!
!  ============================================================================
!        Purpose: This program computes the error function
!                 erf(x) using subroutine ERROR
!   License:
!
!     This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!     they give permission to incorporate this routine into a user program
!     provided that the copyright is acknowledged.
!
!   Author:
!
!     Shanjie Zhang, Jianming Jin
!
!   Reference:
!
!     Shanjie Zhang, Jianming Jin,
!     Computation of Special Functions,
!     Wiley, 1996,
!     ISBN: 0-471-11963-6,
!     LC: QA351.C45.
!
!   Modified on:
!     25 July 2017, by Yingyi Liu
!
!   Parameters:
!        Input:   x   --- Argument of erf(x)
!        Output:  ERR --- erf(x)
!        Example:
!                   x         erf(x)
!                 ---------------------
!                  1.0       .84270079
!                  2.0       .99532227
!                  3.0       .99997791
!                  4.0       .99999998
!                  5.0      1.00000000
!   ============================================================================

        SUBROUTINE ERROR(X,ERR)

        IMPLICIT NONE
        INTEGER K
        REAL*8 EPS,PI,X,X2,ER,R,C0,ERR

        EPS=1.0D-15
        PI=3.141592653589793D0
        X2=X*X
        IF (DABS(X).LT.3.5D0) THEN
           ER=1.0D0
           R=1.0D0
           DO K=1,50
              R=R*X2/(K+0.5D0)
              ER=ER+R
              IF (DABS(R).LE.DABS(ER)*EPS) EXIT
           ENDDO
           C0=2.0D0/DSQRT(PI)*X*DEXP(-X2)
           ERR=C0*ER
        ELSE
           ER=1.0D0
           R=1.0D0
           DO K=1,12
              R=-R*(K-0.5D0)/X2
              ER=ER+R
           ENDDO
           C0=DEXP(-X2)/(DABS(X)*DSQRT(PI))
           ERR=1.0D0-C0*ER
           IF (X.LT.0.0) ERR=-ERR
        ENDIF
        RETURN
        END

!
!  ============================================================================
!        Purpose: This program computes the exponential integral
!                 En(x) using subroutine ENXB
!   License:
!
!     This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!     they give permission to incorporate this routine into a user program
!     provided that the copyright is acknowledged.
!
!   Author:
!
!     Shanjie Zhang, Jianming Jin
!
!   Reference:
!
!     Shanjie Zhang, Jianming Jin,
!     Computation of Special Functions,
!     Wiley, 1996,
!     ISBN: 0-471-11963-6,
!     LC: QA351.C45.
!
!   Modified on:
!     25 July 2017, by Yingyi Liu
!
!   Parameters:
!        Input :  x --- Argument of En(x)
!                 n --- Order of En(x)  (n = 0,1,2,...)
!        Output:  EN --- En(x)
!        Example: x = 10.0
!                    n         En(x)
!                  ----------------------
!                    0     .45399930D-05
!                    1     .41569689D-05
!                    2     .38302405D-05
!                    3     .35487626D-05
!                    4     .33041014D-05
!                    5     .30897289D-05
!  ============================================================================
!
        SUBROUTINE ENXA(N,X,EN)

        IMPLICIT NONE
        INTEGER N,K
        REAL*8 X,E0,E1,EN

        E0=DEXP(-X)/X

        CALL E1XA(X,E1)

        IF (N.EQ.0) THEN
          EN=E0
        ELSEIF (N.EQ.1) THEN
          EN=E1
        ELSE
          DO K=2,N
             EN=(DEXP(-X)-X*E1)/(K-1.0D0)
             E1=EN
          ENDDO
        ENDIF

        RETURN
        END

!
!  ============================================================================
!        Purpose: This program computes the exponential integral
!                 E1(x) using subroutine E1XA
!   License:
!
!     This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!     they give permission to incorporate this routine into a user program
!     provided that the copyright is acknowledged.
!
!   Author:
!
!     Shanjie Zhang, Jianming Jin
!
!   Reference:
!
!     Shanjie Zhang, Jianming Jin,
!     Computation of Special Functions,
!     Wiley, 1996,
!     ISBN: 0-471-11963-6,
!     LC: QA351.C45.
!
!   Modified on:
!     25 July 2017, by Yingyi Liu
!
!   Parameters:
!        Input :  x  --- Argument of E1(x)  ( x > 0 )
!        Output:  E1 --- E1(x)
!        Example:
!                   x        E1(x)
!                 ----------------------
!                  0.0     .1000000+301
!                  1.0     .2193839E+00
!                  2.0     .4890051E-01
!                  3.0     .1304838E-01
!                  4.0     .3779352E-02
!                  5.0     .1148296E-02
!  ============================================================================
!
        SUBROUTINE E1XA(X,E1)

        IMPLICIT NONE
        REAL*8 X,E1,ES1,ES2

        IF (X.EQ.0.0) THEN
           E1=1.0D+300
        ELSE IF (X.LE.1.0) THEN
           E1=-DLOG(X)+((((1.07857D-3*X-9.76004D-3)*X+5.519968D-2)*X-0.24991055D0)*X+0.99999193D0)*X-0.57721566D0
        ELSE
           ES1=(((X+8.5733287401D0)*X+18.059016973D0)*X+8.6347608925D0)*X+0.2677737343D0
           ES2=(((X+9.5733223454D0)*X+25.6329561486D0)*X+21.0996530827D0)*X+3.9584969228D0
           E1=DEXP(-X)/X*ES1/ES2
        ENDIF
        RETURN
        END

!
!  ============================================================================
!        Purpose: This program computes the exponential integral
!                 Ei(x) using subroutine EIX
!   License:
!
!     This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!     they give permission to incorporate this routine into a user program
!     provided that the copyright is acknowledged.
!
!   Author:
!
!     Shanjie Zhang, Jianming Jin
!
!   Reference:
!
!     Shanjie Zhang, Jianming Jin,
!     Computation of Special Functions,
!     Wiley, 1996,
!     ISBN: 0-471-11963-6,
!     LC: QA351.C45.
!
!   Modified on:
!     25 July 2017, by Yingyi Liu
!
!   Parameters:
!        Input :  x  --- Argument of Ei(x)
!        Output:  EI --- Ei(x) ( x > 0 )
!        Example:
!                   x        Ei(x)
!                 -----------------------
!                   0    -.10000000+301
!                   1     .18951178E+01
!                   2     .49542344E+01
!                   3     .99338326E+01
!                   4     .19630874E+02
!                   5     .40185275E+02
!  ============================================================================

        SUBROUTINE EIX(X,EI)

        IMPLICIT NONE
        INTEGER K
        REAL*8 X,EI,R,GA

        IF (X.EQ.0.0) THEN
           EI=-1.0D+300
        ELSE IF (X.LE.40.0) THEN
           EI=1.0D0
           R=1.0D0
           DO K=1,100
              R=R*K*X/(K+1.0D0)**2
              EI=EI+R
              IF (DABS(R/EI).LE.1.0D-15) EXIT
           ENDDO
           GA=0.5772156649015328D0
           EI=GA+DLOG(X)+X*EI
        ELSE
           EI=1.0D0
           R=1.0D0
           DO K=1,20
              R=R*K/X
              EI=EI+R
           ENDDO
           EI=DEXP(X)/X*EI
        ENDIF
        RETURN
        END

!
!  ============================================================================
!        Purpose: This program computes the gamma function
!                 GAMMA(x) for x > 0 or its logarithm using
!                 subroutine LGAMA
!   License:
!
!     This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!     they give permission to incorporate this routine into a user program
!     provided that the copyright is acknowledged.
!
!   Author:
!
!     Shanjie Zhang, Jianming Jin
!
!   Reference:
!
!     Shanjie Zhang, Jianming Jin,
!     Computation of Special Functions,
!     Wiley, 1996,
!     ISBN: 0-471-11963-6,
!     LC: QA351.C45.
!
!   Modified on:
!     25 July 2017, by Yingyi Liu
!
!   Parameters:
!        Purpose: Compute gamma function GAMMA(x) or ln[GAMMA(x)]
!        Input:   x  --- Argument of GAMMA(x) ( x > 0 )
!                 KF --- Function code
!                        KF=1 for GAMMA(x); KF=0 for ln[GAMMA(x)]
!        Output:  GL --- GAMMA(x) or ln[GAMMA(x)]
!        Examples:
!                   x           GAMMA(x)
!                 -------------------------
!                  0.5     .1772453851D+01
!                  2.5     .1329340388D+01
!                  5.0     .2400000000D+02
!                  7.5     .1871254306D+04
!                 10.0     .3628800000D+06
!  ============================================================================

        SUBROUTINE LGAMA(KF,X,GL)

        IMPLICIT NONE
        INTEGER K,KF,N
        REAL*8 X,X0,GL,X2,XP,GL0

        REAL*8 A(10)
        DATA A/8.333333333333333D-02,-2.777777777777778D-03,     &
               7.936507936507937D-04,-5.952380952380952D-04,     &
               8.417508417508418D-04,-1.917526917526918D-03,     &
               6.410256410256410D-03,-2.955065359477124D-02,     &
               1.796443723688307D-01,-1.39243221690590D+00/

        X0=X
        IF (X.EQ.1.0.OR.X.EQ.2.0) THEN
           GL=0.0D0
           IF (KF.EQ.1) GL=DEXP(GL)
        ELSE IF (X.LE.7.0) THEN
           N=INT(7-X)
           X0=X+N
           X2=1.0D0/(X0*X0)
           XP=6.283185307179586477D0
           GL0=A(10)
           DO K=9,1,-1
            GL0=GL0*X2+A(K)
           ENDDO
           GL=GL0/X0+0.5D0*DLOG(XP)+(X0-.5D0)*DLOG(X0)-X0
           IF (X.LE.7.0) THEN
              DO K=1,N
                 GL=GL-DLOG(X0-1.0D0)
                 X0=X0-1.0D0
              ENDDO
           ENDIF
           IF (KF.EQ.1) GL=DEXP(GL)
        ENDIF

        RETURN
        END




!-----------------------------------------------------------------------------
!+ Finds the limit of a series
      SUBROUTINE Limes      &
      (   N,        &      ! in
          S,        &      ! inout
          rLimes,      &      ! out
          i_Pade,      &      ! out
          k_Pade,     &       ! out
          err  )             ! out

! Description:
!   Finds the limit of a series in the case where only
!   the first N+1 terms are known.
!
! Method:
!   The subroutine operates by applying the epsilon-algorithm
!   to the sequence of partial sums of a seris supplied on input.
!   For description of the algorithm, please see:
!
!   [1] T. Mishonov and E. Penev, Int. J. Mod. Phys. B 14, 3831 (2000)
!
! Owners: Todor Mishonov & Evgeni Penev
!
! History:
! Version   Date         Comment
! =======   ====         =======
! 1.0       01/04/2000   Original code. T. Mishonov & E. Penev
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!                        Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
       IMPLICIT NONE

!* Subroutine arguments
!  Scalar arguments with intent(in):
       INTEGER,   INTENT (IN)    :: N       ! width of the epsilon-table

!  Array arguments with intent(inout):
       REAL*8, INTENT (INOUT) :: S(0:N)   ! sequential row of the epsilon-table

!  Scalar arguments with intent(out)
       REAL*8, INTENT (OUT)   :: rLimes  ! value of the series limes
       INTEGER,   INTENT (OUT)   :: i_Pade  ! power of the numerator
       INTEGER,   INTENT (OUT)   :: k_Pade  ! power of the denominator
       REAL*8, INTENT (OUT)   :: err     ! empirical error

!* End of Subroutine arguments

!  Local parameters                    ! these two need no description ;-)
       REAL*8,    PARAMETER   :: zero = 0.0
       REAL*8,    PARAMETER   :: one  = 1.0

!  Local scalars
       REAL*8                 :: A_max   ! maximum element of A
       INTEGER                   :: i       ! index variable for columns
       INTEGER                   :: k       ! index variable for rows

!  Local arrays
       REAL*8                 :: A(0:N)  ! auxiliary row of the epsilon-table

!- End of header --------------------------------------------------------------

!  Parse input: the algorithm cannot employ more elements than supplied on
!               input, i.e. N <= size(S)
!
      IF ( N > SIZE (S(:)) ) THEN
      WRITE (*, '(A)') '*** Illegal input to Limes: N > size(S)'
       STOP 1
      END IF

!  Algorithm not applicable for N < 2
!
      IF ( N < 2 ) THEN
      WRITE (*, '(A)') '*** Illegal input to Limes: N < 2'
       STOP 2
      END IF

!-----------------------------------------------------------------------------
!  I. Initialize with natural assignments
!-----------------------------------------------------------------------------

      rLimes = S(N)                    ! the N-th partial sum
      err    = ABS ( S(N) - S(N-1) )   ! error -> |S(N) - S(N-1)|
      i_Pade = N                       ! Pade approximant [N/0]
      k_Pade = 0                       !
      A(:)   = zero                ! auxiliary row initially set to zero
      A_max  = zero                ! max. element set to zero
      k      = 1                   ! algorithm starts from the first row

!-----------------------------------------------------------------------------
! II. Main loop: fill in the epsilon table, check for convergence ...
!     (provision against division by zero employs pseudo-inverse numbers)
!-----------------------------------------------------------------------------
      DO
       IF ( N - 2 * k + 1 < 0 ) EXIT

! Update the auxiliary row A(i) of the epsilon-table
! by applying the "cross rule".
!
      DO i=0, N - 2 * k + 1
       IF ( S(i+1) /= S(i) ) THEN
        A(i) = A(i+1) + one/(S(i+1) - S(i))
       ELSE
        A(i) = A(i+1)
       END IF
      END DO
      IF ( N - 2 * k < 0 ) EXIT

!  Update the sequential row S(i) of the epsilon-table
!  by applying the "cross rule".
!
      DO i=0, N - 2 * k
       IF ( A(i+1) /= A(i) ) THEN
        S(i) = S(i+1) + one/(A(i+1) - A(i))
       ELSE
        S(i) = S(i+1)
       END IF

!  Check for convergence, based on A_max; see Ref. [1]
!
       IF ( ABS ( A(i) ) > A_max ) THEN
        A_max  = ABS ( A(i) )
        rLimes = S(i)
        k_Pade = k
        i_Pade = i + k_Pade
        err    = one/A_max
        IF ( S(i+1) == S(i) ) RETURN
       END IF
      END DO
       k = k + 1      ! increment row index
      END DO

      END SUBROUTINE Limes


!
!  ============================================================================
!        Purpose: Compute modified Bessel functions I0(x), I1(1),
!                 K0(x) and K1(x), and their derivatives
!   License:
!
!     This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!     they give permission to incorporate this routine into a user program
!     provided that the copyright is acknowledged.
!
!   Author:
!
!     Shanjie Zhang, Jianming Jin
!
!   Reference:
!
!     Shanjie Zhang, Jianming Jin,
!     Computation of Special Functions,
!     Wiley, 1996,
!     ISBN: 0-471-11963-6,
!     LC: QA351.C45.
!
!   Modified on:
!     25 July 2017, by Yingyi Liu
!
!   Parameters:
!        Input :  x   --- Argument ( x � 0 )
!        Output:  BK0 --- K0(x)
!                 BK1 --- K1(x)
!  ============================================================================
!

        SUBROUTINE IK01A(X,BK0,BK1)

        IMPLICIT NONE
        INTEGER K,K0
        REAL*8 PI,EL,X,X2,BI0,BI1,BK0,BK1,R,XR,CA,CT,W0,WW,CB,XR2
        REAL*8 A(12),B(12),A1(8)

        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        X2=X*X
        IF (X.EQ.0.0D0) THEN
           BI0=1.0D0
           BI1=0.0D0
           BK0=1.0D+300
           BK1=1.0D+300
           RETURN
        ELSE IF (X.LE.18.0D0) THEN
           BI0=1.0D0
           R=1.0D0
           DO K=1,50
              R=0.25D0*R*X2/(K*K)
              BI0=BI0+R
              IF (DABS(R/BI0).LT.1.0D-15) EXIT
           ENDDO
           BI1=1.0D0
           R=1.0D0
           DO K=1,50
              R=0.25D0*R*X2/(K*(K+1))
              BI1=BI1+R
              IF (DABS(R/BI1).LT.1.0D-15) EXIT
           ENDDO
           BI1=0.5D0*X*BI1
        ELSE
           DATA A/0.125D0,7.03125D-2,                 &
                  7.32421875D-2,1.1215209960938D-1,     &
                  2.2710800170898D-1,5.7250142097473D-1,     &
                  1.7277275025845D0,6.0740420012735D0,     &
                  2.4380529699556D01,1.1001714026925D02,     &
                  5.5133589612202D02,3.0380905109224D03/
           DATA B/-0.375D0,-1.171875D-1,               &
                  -1.025390625D-1,-1.4419555664063D-1,     &
                  -2.7757644653320D-1,-6.7659258842468D-1,     &
                  -1.9935317337513D0,-6.8839142681099D0,     &
                  -2.7248827311269D01,-1.2159789187654D02,     &
                  -6.0384407670507D02,-3.3022722944809D03/
           K0=12
           IF (X.GE.35.0) K0=9
           IF (X.GE.50.0) K0=7
           CA=DEXP(X)/DSQRT(2.0D0*PI*X)
           BI0=1.0D0
           XR=1.0D0/X
           DO K=1,K0
              BI0=BI0+A(K)*XR**K
           ENDDO
           BI0=CA*BI0
           BI1=1.0D0
           DO K=1,K0
              BI1=BI1+B(K)*XR**K
           ENDDO
           BI1=CA*BI1
        ENDIF
        IF (X.LE.9.0D0) THEN
           CT=-(DLOG(X/2.0D0)+EL)
           BK0=0.0D0
           W0=0.0D0
           R=1.0D0
           DO K=1,50
              W0=W0+1.0D0/K
              R=0.25D0*R/(K*K)*X2
              BK0=BK0+R*(W0+CT)
              IF (DABS((BK0-WW)/BK0).LT.1.0D-15) EXIT
              WW=BK0
           ENDDO
           BK0=BK0+CT
        ELSE
           DATA A1/0.125D0,0.2109375D0,                &
                   1.0986328125D0,1.1775970458984D01,     &
                   2.1461706161499D02,5.9511522710323D03,     &
                   2.3347645606175D05,1.2312234987631D07/
           CB=0.5D0/X
           XR2=1.0D0/X2
           BK0=1.0D0
           DO K=1,8
            BK0=BK0+A1(K)*XR2**K
           ENDDO
           BK0=CB*BK0/BI0
        ENDIF
        BK1=(1.0D0/X-BI1*BK0)/BI0

        RETURN
        END

!
!  ============================================================================
!        Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
!                 Y1(x), and their derivatives
!   License:
!
!     This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!     they give permission to incorporate this routine into a user program
!     provided that the copyright is acknowledged.
!
!   Author:
!
!     Shanjie Zhang, Jianming Jin
!
!   Reference:
!
!     Shanjie Zhang, Jianming Jin,
!     Computation of Special Functions,
!     Wiley, 1996,
!     ISBN: 0-471-11963-6,
!     LC: QA351.C45.
!
!   Modified on:
!     25 July 2017, by Yingyi Liu
!
!   Parameters:
!        Input :  x   --- Argument of Jn(x) & Yn(x) ( x � 0 )
!        Output:  BJ0 --- J0(x)
!                 BJ1 --- J1(x)
!                 BY0 --- Y0(x)
!                 BY1 --- Y1(x)
!  ============================================================================
!
        SUBROUTINE JY01B(X,BJ0,BJ1,BY0,BY1)

        IMPLICIT NONE
        REAL*8 PI,X,BJ0,BY0,T,T2,A0,P0,Q0,TA0,BJ1,BY1,P1,Q1,TA1

        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           BJ0=1.0D0
           BJ1=0.0D0
           BY0=-1.0D+300
           BY1=-1.0D+300
           RETURN
        ELSE IF (X.LE.4.0D0) THEN
           T=X/4.0D0
           T2=T*T
           BJ0=((((((-.5014415D-3*T2+.76771853D-2)*T2-.0709253492D0)*T2+.4443584263D0)    &
                *T2-1.7777560599D0)*T2+3.9999973021D0)*T2-3.9999998721D0)*T2+1.0D0
           BJ1=T*(((((((-.1289769D-3*T2+.22069155D-2)*T2-.0236616773D0)*T2+.1777582922D0)  &
                *T2-.8888839649D0)*T2+2.6666660544D0)*T2-3.9999999710D0)*T2+1.9999999998D0)
           BY0=(((((((-.567433D-4*T2+.859977D-3)*T2-.94855882D-2)*T2+.0772975809D0)*T2-.4261737419D0)  &
               *T2+1.4216421221D0)*T2-2.3498519931D0)*T2+1.0766115157)*T2+.3674669052D0
           BY0=2.0D0/PI*DLOG(X/2.0D0)*BJ0+BY0
           BY1=((((((((.6535773D-3*T2-.0108175626D0)*T2+.107657606D0)*T2-.7268945577D0)*T2+3.1261399273D0)  &
               *T2-7.3980241381D0)*T2+6.8529236342D0)*T2+.3932562018D0)*T2-.6366197726D0)/X
           BY1=2.0D0/PI*DLOG(X/2.0D0)*BJ1+BY1
        ELSE
           T=4.0D0/X
           T2=T*T
           A0=DSQRT(2.0D0/(PI*X))
           P0=((((-.9285D-5*T2+.43506D-4)*T2-.122226D-3)*T2+.434725D-3)*T2-.4394275D-2)*T2+.999999997D0
           Q0=T*(((((.8099D-5*T2-.35614D-4)*T2+.85844D-4)*T2-.218024D-3)*T2+.1144106D-2)*T2-.031249995D0)
           TA0=X-.25D0*PI
           BJ0=A0*(P0*DCOS(TA0)-Q0*DSIN(TA0))
           BY0=A0*(P0*DSIN(TA0)+Q0*DCOS(TA0))
           P1=((((.10632D-4*T2-.50363D-4)*T2+.145575D-3)*T2-.559487D-3)*T2+.7323931D-2)*T2+1.000000004D0
           Q1=T*(((((-.9173D-5*T2+.40658D-4)*T2-.99941D-4)*T2+.266891D-3)*T2-.1601836D-2)*T2+.093749994D0)
           TA1=X-.75D0*PI
           BJ1=A0*(P1*DCOS(TA1)-Q1*DSIN(TA1))
           BY1=A0*(P1*DSIN(TA1)+Q1*DCOS(TA1))
        ENDIF

        RETURN
        END

!
!  ============================================================================
!        Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
!                 Y1(x), and their derivatives
!   License:
!
!     This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!     they give permission to incorporate this routine into a user program
!     provided that the copyright is acknowledged.
!
!   Author:
!
!     Shanjie Zhang, Jianming Jin
!
!   Reference:
!
!     Shanjie Zhang, Jianming Jin,
!     Computation of Special Functions,
!     Wiley, 1996,
!     ISBN: 0-471-11963-6,
!     LC: QA351.C45.
!
!   Modified on:
!     25 July 2017, by Yingyi Liu
!
!   Description:
!     This is a derived version from the subroutine JY01B.
!
!   Parameters:
!        Input :  x   --- Argument of Jn(x) & Yn(x) ( x � 0 )
!        Output:  BJ0 --- J0(x)
!                 BJ1 --- J1(x)
!  ============================================================================
!
        SUBROUTINE JY01BJ(X,BJ0,BJ1)

        IMPLICIT NONE
        REAL*8 PI,X,BJ0,BJ1,T,T2,A0,P0,Q0,TA0,P1,Q1,TA1

        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           BJ0=1.0D0
           BJ1=0.0D0
           RETURN
        ELSE IF (X.LE.4.0D0) THEN
           T=X/4.0D0
           T2=T*T
           BJ0=((((((-.5014415D-3*T2+.76771853D-2)*T2-.0709253492D0)*T2+.4443584263D0)    &
               *T2-1.7777560599D0)*T2+3.9999973021D0)*T2-3.9999998721D0)*T2+1.0D0
           BJ1=T*(((((((-.1289769D-3*T2+.22069155D-2)*T2-.0236616773D0)*T2+.1777582922D0)   &
               *T2-.8888839649D0)*T2+2.6666660544D0)*T2-3.9999999710D0)*T2+1.9999999998D0)
        ELSE
           T=4.0D0/X
           T2=T*T
           A0=DSQRT(2.0D0/(PI*X))
           P0=((((-.9285D-5*T2+.43506D-4)*T2-.122226D-3)*T2+.434725D-3)*T2-.4394275D-2)*T2+.999999997D0
           Q0=T*(((((.8099D-5*T2-.35614D-4)*T2+.85844D-4)*T2-.218024D-3)*T2+.1144106D-2)*T2-.031249995D0)
           TA0=X-.25D0*PI
           BJ0=A0*(P0*DCOS(TA0)-Q0*DSIN(TA0))
           P1=((((.10632D-4*T2-.50363D-4)*T2+.145575D-3)*T2-.559487D-3)*T2+.7323931D-2)*T2+1.000000004D0
           Q1=T*(((((-.9173D-5*T2+.40658D-4)*T2-.99941D-4)*T2+.266891D-3)*T2-.1601836D-2)*T2+.093749994D0)
           TA1=X-.75D0*PI
           BJ1=A0*(P1*DCOS(TA1)-Q1*DSIN(TA1))
        ENDIF

        RETURN
        END

!
!  ============================================================================
!        Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
!                 Y1(x), and their derivatives
!   License:
!
!     This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!     they give permission to incorporate this routine into a user program
!     provided that the copyright is acknowledged.
!
!   Author:
!
!     Shanjie Zhang, Jianming Jin
!
!   Reference:
!
!     Shanjie Zhang, Jianming Jin,
!     Computation of Special Functions,
!     Wiley, 1996,
!     ISBN: 0-471-11963-6,
!     LC: QA351.C45.
!
!   Modified on:
!     25 July 2017, by Yingyi Liu
!
!   Description:
!     This is a derived version from the subroutine JY01B.
!
!   Parameters:
!        Input :  x   --- Argument of Jn(x) & Yn(x) ( x � 0 )
!        Output:  BJ0 --- J0(x)
!                 BJ1 --- J1(x)
!                 BY0 --- Y0(x)
!                 BY1 --- Y1(x)
!  ============================================================================
!
        SUBROUTINE JY01BY0(X,BY0)

        IMPLICIT NONE
        REAL*8 PI,X,BJ0,BY0,T,T2,A0,P0,Q0,TA0

        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           BJ0=1.0D0
           BY0=-1.0D+300
           RETURN
        ELSE IF (X.LE.4.0D0) THEN
           T=X/4.0D0
           T2=T*T
           BJ0=((((((-.5014415D-3*T2+.76771853D-2)*T2-.0709253492D0)*T2+.4443584263D0)   &
               *T2-1.7777560599D0)*T2+3.9999973021D0)*T2-3.9999998721D0)*T2+1.0D0
           BY0=(((((((-.567433D-4*T2+.859977D-3)*T2-.94855882D-2)*T2+.0772975809D0)*T2    &
                -.4261737419D0)*T2+1.4216421221D0)*T2-2.3498519931D0)*T2+1.0766115157)*T2+.3674669052D0
           BY0=2.0D0/PI*DLOG(X/2.0D0)*BJ0+BY0
        ELSE
           T=4.0D0/X
           T2=T*T
           A0=DSQRT(2.0D0/(PI*X))
           P0=((((-.9285D-5*T2+.43506D-4)*T2-.122226D-3)*T2+.434725D-3)*T2-.4394275D-2)*T2+.999999997D0
           Q0=T*(((((.8099D-5*T2-.35614D-4)*T2+.85844D-4)*T2-.218024D-3)*T2+.1144106D-2)*T2-.031249995D0)
           TA0=X-.25D0*PI
           BY0=A0*(P0*DSIN(TA0)+Q0*DCOS(TA0))
        ENDIF

        RETURN
        END

!
!  ============================================================================
!        Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
!                 Y1(x), and their derivatives
!   License:
!
!     This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!     they give permission to incorporate this routine into a user program
!     provided that the copyright is acknowledged.
!
!   Author:
!
!     Shanjie Zhang, Jianming Jin
!
!   Reference:
!
!     Shanjie Zhang, Jianming Jin,
!     Computation of Special Functions,
!     Wiley, 1996,
!     ISBN: 0-471-11963-6,
!     LC: QA351.C45.
!
!   Modified on:
!     25 July 2017, by Yingyi Liu
!
!   Description:
!     This is a derived version from the subroutine JY01B.
!
!   Parameters:
!        Input :  x   --- Argument of Jn(x) & Yn(x) ( x � 0 )
!        Output:  BY1 --- Y1(x)
!  ============================================================================
!
        SUBROUTINE JY01BY1(X,BY1)

        IMPLICIT NONE
        REAL*8 PI,X,BJ1,BY1,T,T2,A0,P1,Q1,TA1

        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           BJ1=0.0D0
           BY1=-1.0D+300
           RETURN
        ELSE IF (X.LE.4.0D0) THEN
           T=X/4.0D0
           T2=T*T
           BJ1=T*(((((((-.1289769D-3*T2+.22069155D-2)*T2-.0236616773D0)*T2+.1777582922D0)   &
               *T2-.8888839649D0)*T2+2.6666660544D0)*T2-3.9999999710D0)*T2+1.9999999998D0)
           BY1=((((((((.6535773D-3*T2-.0108175626D0)*T2+.107657606D0)*T2-.7268945577D0)*T2+3.1261399273D0)    &
                *T2-7.3980241381D0)*T2+6.8529236342D0)*T2+.3932562018D0)*T2-.6366197726D0)/X
           BY1=2.0D0/PI*DLOG(X/2.0D0)*BJ1+BY1
        ELSE
           T=4.0D0/X
           T2=T*T
           A0=DSQRT(2.0D0/(PI*X))
           P1=((((.10632D-4*T2-.50363D-4)*T2+.145575D-3)*T2-.559487D-3)*T2+.7323931D-2)*T2+1.000000004D0
           Q1=T*(((((-.9173D-5*T2+.40658D-4)*T2-.99941D-4)*T2+.266891D-3)*T2-.1601836D-2)*T2+.093749994D0)
           TA1=X-.75D0*PI
           BY1=A0*(P1*DSIN(TA1)+Q1*DCOS(TA1))
        ENDIF

        RETURN
        END

end module FinGrnExtSubs_module
