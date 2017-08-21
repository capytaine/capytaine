MODULE Green_1

  IMPLICIT NONE

  REAL, PARAMETER :: PI = 3.141592653588979 ! π

  ! The index of the following node when going around a face.
  INTEGER, PRIVATE, DIMENSION(4), PARAMETER :: NEXT_NODE = (/ 2, 3 ,4 ,1 /)

CONTAINS

  SUBROUTINE COMPUTE_S0 &
      (M, &
      Face_nodes, Face_Center, Face_Normal, Face_area, Face_radius, &
      S0, VS0)
    ! Estimate the integral over the face S0 = ∫∫ 1/MM' dS(M')
    ! and its derivative with respect to M.

    ! Based on formulas A6.1 and A6.3 (p. 381 to 383)
    ! in G. Delhommeau thesis (referenced below as [Del]).

    ! Inputs
    REAL, DIMENSION(3),    INTENT(IN) :: M
    REAL, DIMENSION(3, 4), INTENT(IN) :: Face_nodes
    REAL, DIMENSION(3),    INTENT(IN) :: Face_center, Face_Normal
    REAL,                  INTENT(IN) :: Face_area, Face_radius

    ! Outputs
    REAL,               INTENT(OUT) :: S0
    REAL, DIMENSION(3), INTENT(OUT) :: VS0

    ! Local variables
    INTEGER               :: L
    REAL                  :: RO, GZ, DK, GY
    REAL, DIMENSION(4)    :: RR
    REAL, DIMENSION(3, 4) :: DRX
    REAL                  :: ANT, DNT, ANL, DNL, ALDEN, AT
    REAL, DIMENSION(3)    :: PJ, GYX, ANTX, ANLX, DNTX

    RO = NORM2(M(1:3) - Face_center(1:3)) ! Distance from center of mass of the face to M.

    IF (RO > 7*Face_radius) THEN
      ! Asymptotic value if face far away from M
      S0       = Face_area/RO
      VS0(1:3) = (Face_center(1:3) - M)*S0/RO**2

    ELSE

      GZ = DOT_PRODUCT(M(1:3) - Face_center(1:3), Face_normal(1:3)) ! Called Z in [Del]

      DO L = 1, 4
        RR(L) = NORM2(M(1:3) - Face_nodes(1:3, L))       ! Distance from vertices of Face to M.
        DRX(:, L) = (M(1:3) - Face_nodes(1:3, L))/RR(L)  ! Normed vector from vertices of Face to M.
      END DO

      S0 = 0.0
      VS0(:) = 0.0

      DO L = 1, 4
        DK = NORM2(Face_nodes(:, NEXT_NODE(L)) - Face_nodes(:, L))    ! Distance between two consecutive points, called d_k in [Del]
        IF (DK >= 1E-3*Face_radius) THEN
          PJ(:) = (Face_nodes(:, NEXT_NODE(L)) - Face_nodes(:, L))/DK ! Normed vector from one corner to the next
          ! Called (a,b,c) in [Del]
          GYX(1) = Face_normal(2)*PJ(3) - Face_normal(3)*PJ(2)
          GYX(2) = Face_normal(3)*PJ(1) - Face_normal(1)*PJ(3)
          GYX(3) = Face_normal(1)*PJ(2) - Face_normal(2)*PJ(1)
          GY = DOT_PRODUCT(M - Face_nodes(:, L), GYX) ! Called Y_k in  [Del]

          ANT = 2*GY*DK                                                  ! Called N^t_k in [Del]
          DNT = (RR(NEXT_NODE(L))+RR(L))**2 - DK*DK + 2.0*ABS(GZ)*(RR(NEXT_NODE(L))+RR(L)) ! Called D^t_k in [Del]
          ANL = RR(NEXT_NODE(L)) + RR(L) + DK                                     ! Called N^l_k in [Del]
          DNL = RR(NEXT_NODE(L)) + RR(L) - DK                                     ! Called D^l_k in [Del]
          ALDEN = ALOG(ANL/DNL)

          IF (ABS(GZ) >= 1.E-4*Face_radius) THEN
            AT = ATAN(ANT/DNT)
          ELSE
            AT = 0.
          ENDIF

          ANLX(:) = DRX(:, NEXT_NODE(L)) + DRX(:, L)                    ! Called N^l_k_{x,y,z} in [Del]

          ANTX(:) = 2*DK*GYX(:)                                ! Called N^t_k_{x,y,z} in [Del]
          DNTX(:) = 2*(RR(NEXT_NODE(L)) + RR(L) + ABS(GZ))*ANLX(:) &
            + 2*SIGN(1.0, GZ)*(RR(NEXT_NODE(L)) + RR(L))*Face_normal(:) ! Called D^t_k_{x,y,z} in [Del]

          S0 = S0 + GY*ALDEN - 2*AT*ABS(GZ)

          VS0(:) = VS0(:) + ALDEN*GYX(:)     &
            - 2*SIGN(1.0, GZ)*AT*Face_normal(:)   &
            + GY*(DNL-ANL)/(ANL*DNL)*ANLX(:) &
            - 2*ABS(GZ)*(ANTX(:)*DNT - DNTX(:)*ANT)/(ANT*ANT+DNT*DNT)
        END IF
      END DO
    END IF

  END SUBROUTINE COMPUTE_S0

  ! =========================

  SUBROUTINE COMPUTE_ASYMPTOTIC_S0 &
      (M, &
      Face_Center, Face_area, &
      S0, VS0)
    ! Same as above, but always use the approximate aymptotic value.

    ! Inputs
    REAL, DIMENSION(3),    INTENT(IN) :: M
    REAL, DIMENSION(3),    INTENT(IN) :: Face_center
    REAL,                  INTENT(IN) :: Face_area

    ! Outputs
    REAL,               INTENT(OUT) :: S0
    REAL, DIMENSION(3), INTENT(OUT) :: VS0

    ! Local variables
    REAL                  :: RO

    RO = NORM2(M(1:3) - Face_center(1:3)) ! Distance from center of mass of the face to M.

    IF (RO > 1e-7) THEN
      ! Asymptotic value if face far away from M
      S0       = Face_area/RO
      VS0(1:3) = (Face_center(1:3) - M)*S0/RO**2
    ELSE
      S0 = 0.0
      VS0(1:3) = 0.0
    END IF

  END SUBROUTINE COMPUTE_ASYMPTOTIC_S0

  ! ===============================

  SUBROUTINE VAV                  &
      ( X0I, &
      Face_nodes, Face_Center, Face_Normal, Face_area, Face_radius, &
      depth, MK, &
      FSP, FSM, VSP, VSM)
    ! Main subroutine of the module, called in SOLVE_BEM.f90 and FREESURFACE.f90.

    ! Inputs
    REAL, DIMENSION(3),    INTENT(IN) :: X0I  ! Coordinates of the source point.
    REAL, DIMENSION(3, 4), INTENT(IN) :: Face_nodes
    REAL, DIMENSION(3),    INTENT(IN) :: Face_center, Face_Normal
    REAL,                  INTENT(IN) :: Face_area, Face_radius
    REAL,                  INTENT(IN) :: depth
    INTEGER,               INTENT(IN) :: MK

    ! Outputs
    REAL,               INTENT(OUT) :: FSP, FSM ! Integral of the Green function over the panel.
    REAL, DIMENSION(3), INTENT(OUT) :: VSP, VSM ! Gradient of the integral of the Green function with respect to X0I.

    ! Local variables
    REAL, DIMENSION(3)  :: XI
    REAL                :: S0, S1
    REAL, DIMENSION(3)  :: VS0, VS1

    XI(:) = X0I(:)
    CALL COMPUTE_S0(XI, Face_nodes, Face_Center, Face_Normal, Face_area, Face_radius, S0, VS0)

    ! Reflected problem across the free surface/sea bottom.
    XI(:) = X0I(:)
    XI(3) = -X0I(3) - 2*depth
    CALL COMPUTE_S0(XI, Face_nodes, Face_Center, Face_Normal, Face_area, Face_radius, S1, VS1)
    VS1(3) = -VS1(3)

    ! IF (Mesh%Isym == NO_Y_SYMMETRY) THEN
    ! Add up the contributions of the two problems.
    FSP    = -S0     - MK*S1
    VSP(:) = -VS0(:) - MK*VS1(:)
    FSM    = 0.0
    VSM(:) = 0.0

    ! ELSE IF (Mesh%Isym == Y_SYMMETRY) THEN
    !   ! Add up the contributions of the two problems...
    !   FSP    = -S0     - MK*S1
    !   VSP(:) = -VS0(:) - MK*VS1(:)
    !   FSM    = -S0     - MK*S1
    !   VSM(:) = -VS0(:) - MK*VS1(:)

    !   !... and do some more.

    !   ! Reflected problem across the symmetry plane (xOz)
    !   XI(1) = X0I(1)
    !   XI(2) = -X0I(2)
    !   XI(3) = X0I(3)
    !   CALL COMPUTE_S0(XI, Face, S0, VS0)
    !   VS0(2) = -VS0(2)

    !   ! Reflected problem across the symmetry plane (xOz) and across the free surface/sea bottom.
    !   XI(1) = X0I(1)
    !   XI(2) = -X0I(2)
    !   XI(3) = -X0I(3) - 2*depth
    !   CALL COMPUTE_S0(XI, Face, S1, VS1)
    !   VS1(2:3) = -VS1(2:3)

    !   ! Add up the new results.
    !   FSP    = FSP    - S0     - MK*S1
    !   VSP(:) = VSP(:) - VS0(:) - MK*VS1(:)
    !   FSM    = FSM    + S0     + MK*S1
    !   VSM(:) = VSM(:) + VS0(:) + MK*VS1(:)

    ! END IF

    FSP    = FSP/(4*PI)
    FSM    = FSM/(4*PI)
    VSP(:) = VSP(:)/(4*PI)
    VSM(:) = VSM(:)/(4*PI)

    RETURN
  END SUBROUTINE VAV
END MODULE Green_1
