! Copyright 2026 Capytaine developers
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
module Green_Rankine

  use floating_point_precision, only: pre
  use constants

  implicit none

  ! The index of the following node when going around a face.
  integer, private, dimension(4), parameter :: next_node = (/ 2, 3, 4, 1 /)

contains

  ! =================================================================

  pure subroutine integral_of_Rankine            &
      (M,                                        &
      face_nodes, face_center, face_normal,      &
      face_area, face_radius,                    &
      derivative_with_respect_to_first_variable, &
      int_G, int_nabla_G)
    ! When M lies on the face, the gradient is the one-sided limit (see exact_integral_of_Rankine):
    ! - derivative with respect to M (adjoint double layer K, full gradient): limit from the side of the normal of the
    !   face (usually the fluid), which gives the +½ n term of K;
    ! - derivative with respect to the source point (double layer D): limit from the other side, which gives the +½
    !   term of D.

    real(kind=pre), dimension(3),    intent(in) :: M
    real(kind=pre), dimension(4, 3), intent(in) :: face_nodes
    real(kind=pre), dimension(3),    intent(in) :: face_center, face_normal
    real(kind=pre),                  intent(in) :: face_area, face_radius
    logical,                         intent(in) :: derivative_with_respect_to_first_variable

    ! Outputs
    real(kind=pre),               intent(out) :: int_G
    real(kind=pre), dimension(3), intent(out) :: int_nabla_G

    ! Local
    real(kind=pre) :: r0

    r0 = norm2(M(1:3) - face_center(1:3))

    if (r0 > 7*face_radius) then
      call one_point_integral_of_Rankine           &
        (M,                                        &
        face_center, face_area,                    &
        derivative_with_respect_to_first_variable, &
        int_G, int_nabla_G)
    else
      call exact_integral_of_Rankine               &
        (M,                                        &
        face_nodes, face_center, face_normal,      &
        face_area, face_radius,                    &
        derivative_with_respect_to_first_variable, &
        merge(ONE, -ONE, derivative_with_respect_to_first_variable), &
        int_G, int_nabla_G)
    end if
    ! No jump in the far field (r0 > 7*face_radius): M cannot be on the face.
  end subroutine

  ! =================================================================

  pure subroutine exact_integral_of_Rankine      &
      (M,                                        &
      face_nodes, face_center, face_normal,      &
      face_area, face_radius,                    &
      derivative_with_respect_to_first_variable, &
      limit_side,                                &
      int_G, int_nabla_G)
    ! Estimate the integral S0 = ∫∫ 1/MM' dS(M') over a face
    ! and its derivative VS0 with respect to M.

    ! Based on formulas A6.1 and A6.3 (p. 381 to 383)
    ! in G. Delhommeau thesis (referenced below as [Del]).

    ! The gradient is discontinuous across the face: when M crosses the face, its normal component jumps by
    ! -4π (the jump of the single layer potential). When M is in the plane of the face (|GZ| < 1e-4*face_radius),
    ! limit_side selects the returned value:
    !   limit_side =  0: principal value (no jump term),
    !   limit_side = +1: limit when M approaches the face from the side of face_normal,
    !   limit_side = -1: limit when M approaches the face from the other side.
    ! The jump term is added wherever M lies on the face.

    ! Inputs
    real(kind=pre), dimension(3),    intent(in) :: M
    real(kind=pre), dimension(4, 3), intent(in) :: face_nodes
    real(kind=pre), dimension(3),    intent(in) :: face_center, face_normal
    real(kind=pre),                  intent(in) :: face_area, face_radius
    logical,                         intent(in) :: derivative_with_respect_to_first_variable
    real(kind=pre),                  intent(in) :: limit_side

    ! Outputs
    real(kind=pre),               intent(out) :: int_G
    real(kind=pre), dimension(3), intent(out) :: int_nabla_G

    ! Local variables
    integer                         :: L
    real(kind=pre)                  :: GZ, DK, GY
    real(kind=pre), dimension(4)    :: RR
    real(kind=pre), dimension(3, 4) :: DRX
    real(kind=pre)                  :: ANT, DNT, ANL, DNL, ALDEN, AT, AT_side
    real(kind=pre), dimension(3)    :: PJ, GYX, ANTX, ANLX, DNTX

    GZ = dot_product(M(1:3) - face_center(1:3), face_normal(1:3)) ! Called Z in [Del]

    do concurrent (L = 1:4)
      RR(L) = norm2(M(1:3) - face_nodes(L, 1:3))       ! Distance from vertices of face to M.
      DRX(:, L) = (M(1:3) - face_nodes(L, 1:3))/RR(L)  ! Normed vector from vertices of face to M.
    end do

    int_G = ZERO
    int_nabla_G(:) = ZERO

    do L = 1, 4
      DK = norm2(face_nodes(next_node(L), :) - face_nodes(L, :))    ! Distance between two consecutive points, called d_k in [Del]
      if (DK >= real(1e-3, PRE)*face_radius) then
        PJ(:) = (face_nodes(next_node(L), :) - face_nodes(L, :))/DK ! Normed vector from one corner to the next
        ! The following GYX(1:3) are called (a,b,c) in [Del]
        GYX(1) = face_normal(2)*PJ(3) - face_normal(3)*PJ(2)
        GYX(2) = face_normal(3)*PJ(1) - face_normal(1)*PJ(3)
        GYX(3) = face_normal(1)*PJ(2) - face_normal(2)*PJ(1)
        GY = dot_product(M - face_nodes(L, :), GYX)                                    ! Called Y_k in  [Del]

        ANT = 2*GY*DK                                                                  ! Called N^t_k in [Del]
        DNT = (RR(NEXT_NODE(L))+RR(L))**2 - DK*DK + 2*abs(GZ)*(RR(NEXT_NODE(L))+RR(L)) ! Called D^t_k in [Del]
        ANL = RR(NEXT_NODE(L)) + RR(L) + DK                                            ! Called N^l_k in [Del]
        DNL = RR(NEXT_NODE(L)) + RR(L) - DK                                            ! Called D^l_k in [Del]
        ALDEN = log(ANL/DNL)

        if (abs(GZ) >= real(1e-4, PRE)*face_radius) then
          AT = atan(ANT/DNT)
          AT_side = sign(ONE, GZ)
        else if ((limit_side /= ZERO) .and. (abs(GY) > real(1e-8, PRE)*face_radius)) then
          ! In the plane of the face, one-sided limit.
          ! The sum over the edges of atan(ANT/DNT) is half the solid angle under which the face is seen from M:
          ! π inside the face, 0 outside (e.g. on a coplanar neighboring face), π/2 on an edge.
          AT = atan(ANT/DNT)
          AT_side = limit_side
        else
          ! In the plane of the face, principal value.
          ! (Or |GY| small: M is on the line of the edge. Outside of the edge, the term is zero. On the edge, it is ±π/2
          ! depending on the side from which M approaches the edge: zero is their mean and gives π/2 in total.)
          AT = ZERO
          AT_side = ZERO
        endif

        ANLX(:) = DRX(:, NEXT_NODE(L)) + DRX(:, L)                    ! Called N^l_k_{x,y,z} in [Del]

        ANTX(:) = 2*DK*GYX(:)                                         ! Called N^t_k_{x,y,z} in [Del]
        DNTX(:) = 2*(RR(NEXT_NODE(L)) + RR(L) + abs(GZ))*ANLX(:)      &
          + 2*sign(ONE, GZ)*(RR(NEXT_NODE(L)) + RR(L))*face_normal(:) ! Called D^t_k_{x,y,z} in [Del]

        if (abs(GY) < 1e-5) then
          ! Edge case where the singularity is on the boundary of the face (GY = 0, ALDEN = infty).
          ! This case seems to only occur when computating the free surface elevation,
          ! so no fix has been implemented for nabla_G, which is not needed then.
          int_G = int_G - 2*AT*abs(GZ)
        else
          ! General case
          int_G = int_G + GY*ALDEN - 2*AT*abs(GZ)
        end if

        int_nabla_G(:) = int_nabla_G(:) + ALDEN*GYX(:) &
          - 2*AT_side*AT*face_normal(:)                &
          + GY*(DNL-ANL)/(ANL*DNL)*ANLX(:)             &
          - 2*abs(GZ)*(ANTX(:)*DNT - DNTX(:)*ANT)/(ANT*ANT+DNT*DNT)
      end if
    end do

    if (.not. derivative_with_respect_to_first_variable) then
      int_nabla_G(:) = - int_nabla_G(:)
    end if
  end subroutine exact_integral_of_rankine

  ! =================================================================

  pure subroutine one_point_integral_of_Rankine  &
      (M,                                        &
      face_center, face_area,                    &
      derivative_with_respect_to_first_variable, &
      int_G, int_nabla_G)
    ! Same as above, but always use the approximate asymptotic value.

    ! Inputs
    real(kind=pre), dimension(3), intent(in) :: M
    real(kind=pre), dimension(3), intent(in) :: face_center
    real(kind=pre),               intent(in) :: face_area
    logical,                      intent(in) :: derivative_with_respect_to_first_variable

    ! Outputs
    real(kind=pre),               intent(out) :: int_G
    real(kind=pre), dimension(3), intent(out) :: int_nabla_G

    ! Local variables
    real(kind=pre) :: ro
    real(kind=pre), dimension(3) :: vo

    vo = face_center - m
    ro = norm2(vo) ! distance from center of mass of the face to m.

    if (ro > real(1e-7, kind=pre)) then
      ! Asymptotic value if face far away from M
      int_G            = face_area/ro
      int_nabla_G(1:3) = vo*int_G/ro**2
    ELSE
      ! Singularity...
      int_G = ZERO
      int_nabla_G(1:3) = ZERO
    END IF

    if (.not. derivative_with_respect_to_first_variable) then
      int_nabla_G(:) = -int_nabla_G(:)
    endif
  end subroutine one_point_integral_of_Rankine

  ! =================================================================

  pure subroutine integral_of_reflected_Rankine  &
      (M,                                        &
      face_nodes, face_center, face_normal,      &
      face_area, face_radius,                    &
      derivative_with_respect_to_first_variable, &
      reflection_coefs,                          &
      reversed_limit_side,                       &
      int_G, int_nabla_G)

    real(kind=pre), dimension(3),    intent(in) :: M
    real(kind=pre), dimension(4, 3), intent(in) :: face_nodes
    real(kind=pre), dimension(3),    intent(in) :: face_center, face_normal
    real(kind=pre),                  intent(in) :: face_area, face_radius
    logical,                         intent(in) :: derivative_with_respect_to_first_variable
    real(kind=pre), dimension(2),    intent(in) :: reflection_coefs
    logical,                         intent(in) :: reversed_limit_side
    ! See exact_integral_of_reflected_Rankine

    ! Outputs
    real(kind=pre),               intent(out) :: int_G
    real(kind=pre), dimension(3), intent(out) :: int_nabla_G

    ! Local
    real(kind=pre) :: r1

    r1 = hypot(norm2(M(1:2) - face_center(1:2)),  &
               reflection_coefs(1) * M(3) + reflection_coefs(2) - face_center(3))

    if (r1 > 7*face_radius) then
      call one_point_integral_of_reflected_Rankine &
        (M,                                        &
        face_center, face_area,                    &
        derivative_with_respect_to_first_variable, &
        reflection_coefs,                          &
        int_G, int_nabla_G)
    else
      call exact_integral_of_reflected_Rankine     &
        (M,                                        &
        face_nodes, face_center, face_normal,      &
        face_area, face_radius,                    &
        derivative_with_respect_to_first_variable, &
        reflection_coefs,                          &
        reversed_limit_side,                       &
        int_G, int_nabla_G)
    end if
  end subroutine integral_of_reflected_Rankine

  ! =================================================================

  pure subroutine exact_integral_of_reflected_Rankine &
      (M,                                             &
      face_nodes, face_center, face_normal,           &
      face_area, face_radius,                         &
      derivative_with_respect_to_first_variable,      &
      reflection_coefs,                               &
      reversed_limit_side,                            &
      int_g, int_nabla_g)
    ! Reflecting with respect to the free surface

    ! When M lies on a face, it approaches it from the side limit_side, chosen as in integral_of_Rankine
    ! (+1 for the derivative with respect to M, -1 otherwise), or the opposite if reversed_limit_side.
    ! (reversed_limit_side is meant for the free surface image with high_freq, see build_matrices.)
    ! The side reflected_side from which the reflected point approaches the face is:
    !
    ! Derivative with respect to M (adjoint double layer K, or full gradient):
    !   The reflected point is R M + b, with R = diag(1, 1, reflection_coefs(1)). When M = M0 + η n with η -> 0
    !   with the sign of limit_side, the reflected point is (R M0 + b) + η R n, which approaches the plane of the
    !   face from the side sign(n·Rn)*limit_side. This is the genuine one-sided limit of the reflected term.
    !   For instance, for a face on the free surface (n = ±e_z, R n = -n) and the free surface image, the reflected
    !   point comes from the other side, the jump of the image term has the same direction as the one of the
    !   direct term, and the total coefficient is 1 instead of 1/2.
    !
    ! Derivative with respect to the source point (double layer D, direct method):
    !   The coefficient is not the limit of the double layer potential. It is the free term of the boundary integral
    !   equation, to which each singularity of the Green function that lies on the boundary contributes ½.
    !   This is obtained with the reflected side equal to limit_side: for a face on the free surface, the coefficient is 1.
    !   (The genuine one-sided limit would give 0, as the normal derivative of 1/r + 1/r_reflected vanishes on z = 0.)

    ! Inputs
    real(kind=pre), dimension(3),    intent(in) :: M
    real(kind=pre), dimension(4, 3), intent(in) :: face_nodes
    real(kind=pre), dimension(3),    intent(in) :: face_center, face_normal
    real(kind=pre),                  intent(in) :: face_area, face_radius
    logical,                         intent(in) :: derivative_with_respect_to_first_variable
    real(kind=pre), dimension(2),    intent(in) :: reflection_coefs
    logical,                         intent(in) :: reversed_limit_side

    ! Outputs
    real(kind=pre),               intent(out) :: int_g
    real(kind=pre), dimension(3), intent(out) :: int_nabla_g

    ! local
    real(kind=pre), dimension(3)              :: reflected_M
    real(kind=pre)                            :: limit_side, reflected_side, n_dot_Rn

    reflected_M(1:2) = M(1:2)
    reflected_M(3)   = reflection_coefs(1) * M(3) + reflection_coefs(2)

    ! For instance
    ! reflection_coefs = [-1.0, 0.0] for free surface symmetry
    ! reflection_coefs = [-1.0, -2*h] for sea bottom symmetry

    limit_side = merge(ONE, -ONE, derivative_with_respect_to_first_variable)
    if (reversed_limit_side) limit_side = -limit_side

    if (derivative_with_respect_to_first_variable) then
      n_dot_Rn = face_normal(1)**2 + face_normal(2)**2 + reflection_coefs(1)*face_normal(3)**2
      if (abs(n_dot_Rn) < real(1e-6, PRE)) then
        ! The reflected point moves parallel to the face (n·Rn = 0), it does not cross it.
        reflected_side = ZERO
      else
        reflected_side = sign(ONE, n_dot_Rn) * limit_side
      end if
    else
      reflected_side = limit_side
    end if

    call exact_integral_of_Rankine &
      (reflected_M,                &
      face_nodes, face_center,     &
      face_normal, face_area,      &
      face_radius,                 &
      .true.,                      &
      reflected_side,              &
      int_g, int_nabla_g)

    ! If we mirrored M, we mirror the gradient
    int_nabla_g(3) = reflection_coefs(1) * int_nabla_g(3)

    if (.not. derivative_with_respect_to_first_variable) then
      int_nabla_g(1:2) = -int_nabla_g(1:2)
      if (reflection_coefs(1) > 0) then
        int_nabla_G(3) = -int_nabla_g(3)
      endif
    end if

  end subroutine exact_integral_of_reflected_Rankine

  ! =================================================================

  pure subroutine one_point_integral_of_reflected_Rankine &
      (M,                                                 &
      face_center, face_area,                             &
      derivative_with_respect_to_first_variable,          &
      reflection_coefs,                                   &
      int_G, int_nabla_G)

    ! Inputs
    real(kind=pre), dimension(3), intent(in) :: M
    real(kind=pre), dimension(3), intent(in) :: face_center
    real(kind=pre),               intent(in) :: face_area
    logical,                      intent(in) :: derivative_with_respect_to_first_variable
    real(kind=pre), dimension(2), intent(in) :: reflection_coefs

    ! outputs
    real(kind=pre),               intent(out) :: int_G
    real(kind=pre), dimension(3), intent(out) :: int_nabla_G

    ! local
    real(kind=pre), dimension(3)              :: reflected_M

    reflected_M(1:2) = M(1:2)
    reflected_M(3)   = reflection_coefs(1) * M(3) + reflection_coefs(2)

    call one_point_integral_of_Rankine(    &
      reflected_M, face_center, face_area, &
      .true.,                              &
      int_G, int_nabla_G)

    int_nabla_G(3) = reflection_coefs(1) * int_nabla_G(3)

    if (.not. derivative_with_respect_to_first_variable) then
      int_nabla_G(1:2) = -int_nabla_G(1:2)
      if (reflection_coefs(1) > 0) then
        int_nabla_G(3) = -int_nabla_G(3)
      endif
    end if

  end subroutine one_point_integral_of_reflected_Rankine

  ! =================================================================

end module Green_Rankine
