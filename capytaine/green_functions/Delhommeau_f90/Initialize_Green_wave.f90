! Copyright (C) 2017-2019 Matthieu Ancellin
! See LICENSE file at <https://github.com/mancellin/capytaine>
MODULE INITIALIZE_GREEN_WAVE

  USE CONSTANTS

  IMPLICIT NONE

  PUBLIC :: EXP_E1
  PUBLIC :: INITIALIZE_TABULATED_INTEGRALS
  PUBLIC :: FF

CONTAINS

!------------------------------------------------------------------------------

  FUNCTION EXP_E1(ZZ)
    ! Estimation of exp(z)·E1(z) where E1(z) = ∫_z^∞ exp(-t)/t dt
    ! See p.367 of G. Delhommeau thesis (referenced as [Del]).

    ! This routine is only used in the following subroutine to compute the tabulated integrals

    ! The computation is done is single precision.

    COMPLEX(KIND=PRE),      INTENT(IN)  :: ZZ
    COMPLEX(KIND=PRE)                   :: EXP_E1
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
    EXP_E1 = CMPLX(G, KIND=PRE)

  END FUNCTION EXP_E1

!------------------------------------------------------------------------------

  SUBROUTINE INITIALIZE_TABULATED_INTEGRALS(NB_POINTS_X, NB_POINTS_Z, NB_INTEGRATION_POINTS, X, Z, TABULATION)
    ! Compute the tabulated integrals for the wave part of the Green function.
    ! These tables are independent of the mesh and of the frequency.
    ! They are only initialised once at the beginning of the program.

    ! References:
    ! [1] Delhommeau, Amélioration des codes de calcul de diffraction-radiation, 2èmes journées de l'hydrodynamique, 1989
    ! [2] Babarit and Delhommeau, Theoretical and numerical aspects of the open source BEM solver NEMOH, EWTEC 2015
    ! [3] Xie et al., Comparison of existing methods for the calculation of the infinite water depth free-surface Green
    !                 function for the wave-structure interaction problem, Applied Ocean Research, 2018

    ! Inputs
    INTEGER, INTENT(IN) :: NB_POINTS_X     != 328
    INTEGER, INTENT(IN) :: NB_POINTS_Z     != 46
    INTEGER, INTENT(IN) :: NB_INTEGRATION_POINTS

    ! Outputs
    REAL(KIND=PRE), DIMENSION(NB_POINTS_X),                    INTENT(OUT) :: X
    REAL(KIND=PRE), DIMENSION(NB_POINTS_Z),                    INTENT(OUT) :: Z
    REAL(KIND=PRE), DIMENSION(NB_POINTS_X, NB_POINTS_Z, 2, 2), INTENT(OUT) :: TABULATION

    ! Local variables
    INTEGER :: I, J, K
    REAL(KIND=PRE) :: THETA(NB_INTEGRATION_POINTS), CQT(NB_INTEGRATION_POINTS)
    REAL(KIND=PRE) :: COSTHETA
    COMPLEX(KIND=PRE) :: JZETA, ZETA, EXPZETA

    ! Initialize Z axis
    DO J = 1, NB_POINTS_Z
      Z(J) = -AMIN1(10**(J/5.0-6), 10**(J/8.0-4.5), 16.)
    END DO

    ! Initialize X axis
    X(1) = 0.0
    DO I = 2, NB_POINTS_X
      IF (I < 40) THEN
        X(I) = AMIN1(10**((I-1.0)/5-6), 4.0/3.0 + ABS(I-32)/3.0)
      ELSE
        X(I) = 4.0/3.0 + ABS(I-32)/3.0
      ENDIF
    END DO

    ! Initialize THETA and CQT for the integration between -pi/2 and pi/2 with the Simpson rule.
    DO K = 1, NB_INTEGRATION_POINTS
      THETA(K) = -PI/2 + (K-1.0)/(NB_INTEGRATION_POINTS-1.0)*PI  ! The 1.0 are here on purpose to force the recasting of K and NB_INTEGRATION_POINTS as floats.
      IF ((K == 1) .OR. (K == NB_INTEGRATION_POINTS)) THEN
        CQT(K) = PI/(3*(NB_INTEGRATION_POINTS-1))
      ELSEIF (MOD(K,2)==0) THEN
        CQT(K) = 4.0/(3*(NB_INTEGRATION_POINTS-1))*PI
      ELSE
        CQT(K) = 2.0/(3*(NB_INTEGRATION_POINTS-1))*PI
      ENDIF
    ENDDO

    ! Compute tabulation
    TABULATION(:, :, :, :) = 0.0
    DO I = 1, NB_POINTS_X
      DO J = 1, NB_POINTS_Z
        DO K = 1, NB_INTEGRATION_POINTS
          COSTHETA = COS(THETA(K))
          ZETA = CMPLX(Z(J), X(I)*COSTHETA, KIND=PRE)
          IF (REAL(ZETA) <= -30.0) THEN
            EXPZETA = (0.0, 0.0)
          ELSE
            EXPZETA = EXP(ZETA)
          ENDIF
          IF (AIMAG(ZETA) < 0) THEN
            JZETA = EXP_E1(ZETA) - II*PI*EXPZETA
          ELSE
            JZETA = EXP_E1(ZETA) + II*PI*EXPZETA
          END IF
          TABULATION(I, J, 1, 1) = TABULATION(I, J, 1, 1) + CQT(K)*COSTHETA*AIMAG(JZETA - 1.0/ZETA) ! D_1 in [1, 2, 3]
          TABULATION(I, J, 2, 1) = TABULATION(I, J, 2, 1) + CQT(K)*COSTHETA*AIMAG(EXPZETA)          ! D_2 in [1, 2, 3]
#ifdef XIE_CORRECTION
          TABULATION(I, J, 1, 2) = TABULATION(I, J, 1, 2) + CQT(K)*REAL(JZETA)                      ! Z_1 in [3]
#else
          TABULATION(I, J, 1, 2) = TABULATION(I, J, 1, 2) + CQT(K)*REAL(JZETA - 1.0/ZETA)           ! Z_1 in [1, 2]
#endif
          TABULATION(I, J, 2, 2) = TABULATION(I, J, 2, 2) + CQT(K)*REAL(EXPZETA)                    ! Z_2 in [1, 2, 3]
        END DO
      END DO
    END DO

    RETURN

  END SUBROUTINE INITIALIZE_TABULATED_INTEGRALS

!-------------------------------------------------------------------------------!

  FUNCTION FF(XTT, AK, AM)
    ! A function that will be Prony-decomposed for the finite-depth Green function.
    ! See the method "find_best_exponential_decomposition" in green_functions/delhommeau.py.

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
