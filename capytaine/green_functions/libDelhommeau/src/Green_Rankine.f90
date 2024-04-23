! Copyright (C) 2017-2024 Matthieu Ancellin
! See LICENSE file at <https://github.com/capytaine/libDelhommeau>
MODULE GREEN_RANKINE

  USE FLOATING_POINT_PRECISION, ONLY: PRE
  USE CONSTANTS

  IMPLICIT NONE

  ! The index of the following node when going around a face.
  INTEGER, PRIVATE, DIMENSION(4), PARAMETER :: NEXT_NODE = (/ 2, 3, 4, 1 /)

CONTAINS

  PURE SUBROUTINE COMPUTE_INTEGRAL_OF_RANKINE_SOURCE                &
      (M,                                                           &
      Face_nodes, Face_center, Face_normal, Face_area, Face_radius, &
      S0, VS0)
    ! Estimate the integral S0 = ∫∫ 1/MM' dS(M') over a face
    ! and its derivative VS0 with respect to M.

    ! Based on formulas A6.1 and A6.3 (p. 381 to 383)
    ! in G. Delhommeau thesis (referenced below as [Del]).

    ! Inputs
    REAL(KIND=PRE), DIMENSION(3),    INTENT(IN) :: M
    REAL(KIND=PRE), DIMENSION(4, 3), INTENT(IN) :: Face_nodes
    REAL(KIND=PRE), DIMENSION(3),    INTENT(IN) :: Face_center, Face_normal
    REAL(KIND=PRE),                  INTENT(IN) :: Face_area, Face_radius

    ! Outputs
    REAL(KIND=PRE),               INTENT(OUT) :: S0
    REAL(KIND=PRE), DIMENSION(3), INTENT(OUT) :: VS0

    ! Local variables
    INTEGER                         :: L
    REAL(KIND=PRE)                  :: RO, GZ, DK, GY
    REAL(KIND=PRE), DIMENSION(4)    :: RR
    REAL(KIND=PRE), DIMENSION(3, 4) :: DRX
    REAL(KIND=PRE)                  :: ANT, DNT, ANL, DNL, ALDEN, AT
    REAL(KIND=PRE), DIMENSION(3)    :: PJ, GYX, ANTX, ANLX, DNTX

    RO = NORM2(M(1:3) - Face_center(1:3)) ! Distance from center of mass of the face to M.

    IF (RO > 7*Face_radius) THEN
      ! Asymptotic value if face far away from M
      S0       = Face_area/RO
      VS0(1:3) = (Face_center(1:3) - M)*S0/RO**2

    ELSE

      GZ = DOT_PRODUCT(M(1:3) - Face_center(1:3), Face_normal(1:3)) ! Called Z in [Del]

      DO CONCURRENT (L = 1:4)
        RR(L) = NORM2(M(1:3) - Face_nodes(L, 1:3))       ! Distance from vertices of Face to M.
        DRX(:, L) = (M(1:3) - Face_nodes(L, 1:3))/RR(L)  ! Normed vector from vertices of Face to M.
      END DO

      S0 = ZERO
      VS0(:) = ZERO

      DO L = 1, 4
        DK = NORM2(Face_nodes(NEXT_NODE(L), :) - Face_nodes(L, :))    ! Distance between two consecutive points, called d_k in [Del]
        IF (DK >= REAL(1e-3, PRE)*Face_radius) THEN
          PJ(:) = (Face_nodes(NEXT_NODE(L), :) - Face_nodes(L, :))/DK ! Normed vector from one corner to the next
          ! The following GYX(1:3) are called (a,b,c) in [Del]
          GYX(1) = Face_normal(2)*PJ(3) - Face_normal(3)*PJ(2)
          GYX(2) = Face_normal(3)*PJ(1) - Face_normal(1)*PJ(3)
          GYX(3) = Face_normal(1)*PJ(2) - Face_normal(2)*PJ(1)
          GY = DOT_PRODUCT(M - Face_nodes(L, :), GYX)                                    ! Called Y_k in  [Del]

          ANT = 2*GY*DK                                                                  ! Called N^t_k in [Del]
          DNT = (RR(NEXT_NODE(L))+RR(L))**2 - DK*DK + 2*ABS(GZ)*(RR(NEXT_NODE(L))+RR(L)) ! Called D^t_k in [Del]
          ANL = RR(NEXT_NODE(L)) + RR(L) + DK                                            ! Called N^l_k in [Del]
          DNL = RR(NEXT_NODE(L)) + RR(L) - DK                                            ! Called D^l_k in [Del]
          ALDEN = LOG(ANL/DNL)

          IF (ABS(GZ) >= REAL(1e-4, PRE)*Face_radius) THEN
            AT = ATAN(ANT/DNT)
          ELSE
            AT = 0.
          ENDIF

          ANLX(:) = DRX(:, NEXT_NODE(L)) + DRX(:, L)                    ! Called N^l_k_{x,y,z} in [Del]

          ANTX(:) = 2*DK*GYX(:)                                         ! Called N^t_k_{x,y,z} in [Del]
          DNTX(:) = 2*(RR(NEXT_NODE(L)) + RR(L) + ABS(GZ))*ANLX(:) &
            + 2*SIGN(ONE, GZ)*(RR(NEXT_NODE(L)) + RR(L))*Face_normal(:) ! Called D^t_k_{x,y,z} in [Del]

          IF (ABS(GY) < 1e-5) THEN
            ! Edge case where the singularity is on the boundary of the face (GY = 0, ALDEN = infty).
            ! This case seems to only occur when computating the free surface elevation,
            ! so no fix has been implemented for VS0, which is not needed then.
            S0 = S0 - 2*AT*ABS(GZ)
          ELSE
            ! General case
            S0 = S0 + GY*ALDEN - 2*AT*ABS(GZ)
          END IF

          VS0(:) = VS0(:) + ALDEN*GYX(:)     &
            - 2*SIGN(ONE, GZ)*AT*Face_normal(:)   &
            + GY*(DNL-ANL)/(ANL*DNL)*ANLX(:) &
            - 2*ABS(GZ)*(ANTX(:)*DNT - DNTX(:)*ANT)/(ANT*ANT+DNT*DNT)
        END IF
      END DO
    END IF

  END SUBROUTINE COMPUTE_INTEGRAL_OF_RANKINE_SOURCE

  ! =========================

  PURE SUBROUTINE COMPUTE_ASYMPTOTIC_RANKINE_SOURCE &
      (M,                                      &
      Face_center, Face_area,                  &
      S0, VS0)
    ! Same as above, but always use the approximate aymptotic value.

    ! Inputs
    REAL(KIND=PRE), DIMENSION(3), INTENT(IN) :: M
    REAL(KIND=PRE), DIMENSION(3), INTENT(IN) :: Face_center
    REAL(KIND=PRE),               INTENT(IN) :: Face_area

    ! Outputs
    REAL(KIND=PRE),               INTENT(OUT) :: S0
    REAL(KIND=PRE), DIMENSION(3), INTENT(OUT) :: VS0

    ! Local variables
    REAL(KIND=PRE) :: RO
    REAL(KIND=PRE), DIMENSION(3) :: VO

    VO = Face_center - M
    RO = NORM2(VO) ! Distance from center of mass of the face to M.

    IF (RO > REAL(1e-7, KIND=PRE)) THEN
      ! Asymptotic value if face far away from M
      S0       = Face_area/RO
      VS0(1:3) = VO*S0/RO**2
    ELSE
      ! Singularity...
      S0 = ZERO
      VS0(1:3) = ZERO
    END IF

  END SUBROUTINE COMPUTE_ASYMPTOTIC_RANKINE_SOURCE

  ! ====================================

END MODULE GREEN_RANKINE
