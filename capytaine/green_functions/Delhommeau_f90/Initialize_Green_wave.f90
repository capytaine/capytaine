! Copyright (C) 2017-2019 Matthieu Ancellin
! See LICENSE file at <https://github.com/mancellin/capytaine>
MODULE INITIALIZE_GREEN_WAVE

  USE CONSTANTS

  IMPLICIT NONE

  PUBLIC :: GG
  PUBLIC :: INITIALIZE_TABULATED_INTEGRALS
  PUBLIC :: FF

CONTAINS

!------------------------------------------------------------------------------

  FUNCTION GG(ZZ)
    ! Estimation of exp(z)·E1(z) where E1(z) = ∫_z^∞ exp(-t)/t dt
    ! See p.367 of G. Delhommeau thesis (referenced as [Del]).

    ! This routine is only used in the following subroutine to compute the tabulated integrals

    ! The computation is done is single precision.

    COMPLEX(KIND=PRE),      INTENT(IN)  :: ZZ
    COMPLEX(KIND=PRE)                   :: GG
    COMPLEX(KIND=KIND(1e0))             :: G, Y, Z

    Z = CMPLX(ZZ, KIND=KIND(1e0))
    IF (REAL(Z) < -16.0) THEN                                      ! Case 1 p. 368 in [Del]
      Y = 1./Z
      G = Y*(1.+Y*(-1.+Y*(2.+Y*(-6.+Y*(24.+Y*(-120.))))))
    ELSE IF (ABS(AIMAG(Z)) > 10.0) THEN                            ! Case 3 p. 368 in [Del]
      G = 0.711093/(Z+0.415775)+0.278518/(Z+2.29428)+0.010389/(Z+6.2900)
    ELSE IF (REAL(Z) > -0.5) THEN                                  ! Case 2 p. 368 in [Del]
      G = -(LOG(Z)*(.1E+01+Z*(0.23721365E+00+Z*(0.206543E-01+Z*(0.763297E-03+     &
        Z*0.97087007E-05))))+0.5772156649E+00*(0.99999207E+00+                    &
        Z*(-0.149545886E+01+Z*(0.41806426E-01+Z*(-0.3000591E-01+                  &
        Z*(0.19387339E-02+Z*(-0.51801555E-03)))))))/(0.1E+01+                     &
        Z*(-0.76273617E+00+Z*(0.28388363E+00+Z*(-0.66786033E-01+Z*(0.12982719E-01 &
        +Z*(-0.8700861E-03+Z*0.2989204E-03))))))
    ELSE                                                           ! Case 4 p. 369 in [Del]
      IF (AIMAG(Z) < 0) THEN
        G = ((((((( (1.000000, 1.3935496E-06)*Z+ (15.82958, -20.14222))   &
          *Z+ (-70.52863, -227.9511))*Z+ (-985.4221, -226.6272))*Z        &
          + (-1202.318, 1580.907))*Z+ (953.2441, 1342.447))*Z             &
          + (417.3716, -196.6665))*Z+ (-9.881266, -24.24952))/            &
          (((((((( (1.000000, 0.0000000E+00)*Z+ (16.83184, -20.14481))*Z  &
          + (-55.66969, -248.1167))*Z+ (-1068.640, -434.4707))*Z          &
          + (-2082.250, 1522.471))*Z+ (383.3455, 2730.378))*Z             &
          + (1216.791, 351.7189))*Z+ (115.3926, -161.2647))*Z             &
          + (-3.777369, -4.510900))
      ELSE
        G = ((((((( (1.000000, -1.3935496E-06)*Z+ (15.82958, 20.14222))   &
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
    GG = CMPLX(G, KIND=PRE)

  END FUNCTION GG

!------------------------------------------------------------------------------

  SUBROUTINE INITIALIZE_TABULATED_INTEGRALS(IR, JZ, NPINTE, XR, XZ, APD)
    ! Compute the tabulated integral for the wave part of the Green function.
    ! These tables are independant of the mesh and of the frequency.
    ! They are only initialised once at the beginning of the computation.

    ! References:
    ! [1] Delhommeau, Amélioration des codes de calcul de diffraction-radiation, 2èmes journées de l'hydrodynamique, 1989
    ! [2] Babarit and Delhommeau, Theoretical and numerical aspects of the open source BEM solver NEMOH, EWTEC 2015
    ! [3] Xie et al., Comparison of existing methods for the calculation of the infinite water depth free-surface Green
    !                 function for the wave-structure interaction problem, Applied Ocean Research, 2018

    ! Inputs
    INTEGER, INTENT(IN) :: IR     != 328
    INTEGER, INTENT(IN) :: JZ     != 46
    INTEGER, INTENT(IN) :: NPINTE

    ! Outputs
    REAL(KIND=PRE), DIMENSION(IR),           INTENT(OUT) :: XR
    REAL(KIND=PRE), DIMENSION(JZ),           INTENT(OUT) :: XZ
    REAL(KIND=PRE), DIMENSION(IR, JZ, 2, 2), INTENT(OUT) :: APD

    ! Local variables
    INTEGER :: I, J, K
    REAL(KIND=PRE) :: THETA(NPINTE), CQT(NPINTE)
    REAL(KIND=PRE) :: CT
    COMPLEX(KIND=PRE) :: C1, C2, ZETA, CEX

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
      THETA(K) = -PI/2 + (K-1.0)/(NPINTE-1.0)*PI  ! The 1.0 are here on purpose to force the recasting of K and NPINTE as floats.
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
          ZETA = CMPLX(XZ(J), XR(I)*CT, KIND=PRE)
          IF (REAL(ZETA) <= -30.0) THEN
            CEX = (0.0, 0.0)
          ELSE
            CEX = EXP(ZETA)
          ENDIF
          IF (AIMAG(ZETA) < 0) THEN
            C1 = CQT(K)*(GG(ZETA) - II*PI*CEX - 1.0/ZETA)
          ELSE
            C1 = CQT(K)*(GG(ZETA) + II*PI*CEX - 1.0/ZETA)
          END IF
          C2 = CQT(K)*CEX
          APD(I, J, 1, 1) = APD(I, J, 1, 1) + CT*AIMAG(C1) ! named D_1(Z, X) in [1, 2]
          APD(I, J, 1, 2) = APD(I, J, 1, 2) + REAL(C1)     ! named D_2(Z, X) in [1, 2]
          APD(I, J, 2, 1) = APD(I, J, 2, 1) + CT*AIMAG(C2) ! named Z_1(Z, X) in [1, 2]
          APD(I, J, 2, 2) = APD(I, J, 2, 2) + REAL(C2)     ! named Z_2(Z, X) in [1, 2]
        END DO
      END DO
    END DO

    RETURN

  END SUBROUTINE INITIALIZE_TABULATED_INTEGRALS

!-------------------------------------------------------------------------------!

  FUNCTION FF(XTT, AK, AM)
    ! A function that will be Prony-decomposed for the finite-depth Green function.
    ! See the method "compute_exponential_decomposition" in bem/nemoh.py.

    ! Input
    REAL(KIND=PRE), INTENT(IN) :: XTT, AK, AM

    ! Local variables
    REAL(KIND=PRE) :: COEF, TOL, A, B, C, D, E, F, FF

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

!----------------------------------------------------------------

END MODULE INITIALIZE_GREEN_WAVE
