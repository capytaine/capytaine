! Copyright (C) 2017-2025 Matthieu Ancellin
! See LICENSE file at <https://github.com/capytaine/libDelhommeau>
module Green_Rankine

  use floating_point_precision, only: pre
  use constants
  use mesh_types

  implicit none

  ! The index of the following node when going around a face.
  integer, private, dimension(4), parameter :: next_node = (/ 2, 3, 4, 1 /)

contains

  ! =================================================================

  pure subroutine integral_of_Rankine(field_point, source_face, &
                                derivative_with_respect_to_first_variable, &
                                int_G, int_nabla_G)
    real(kind=pre), dimension(3), intent(in) :: field_point
    type(face_type), intent(in) :: source_face
    logical, intent(in) :: derivative_with_respect_to_first_variable
    real(kind=pre), intent(out) :: int_G
    real(kind=pre), dimension(3), intent(out) :: int_nabla_G

    ! Local
    real(kind=pre) :: r0

    r0 = norm2(field_point - source_face%center)

    if (r0 > 7*source_face%radius) then
      call one_point_integral_of_Rankine           &
        (field_point,                              &
        source_face,                               &
        derivative_with_respect_to_first_variable, &
        int_G, int_nabla_G)
    else
      call exact_integral_of_Rankine               &
        (field_point,                              &
        source_face,                               &
        derivative_with_respect_to_first_variable, &
        int_G, int_nabla_G)
    end if
  end subroutine integral_of_Rankine

  ! =================================================================

  pure subroutine exact_integral_of_Rankine      &
      (field_point, face,                       &
      derivative_with_respect_to_first_variable, &
      int_G, int_nabla_G)
    ! Estimate the integral S0 = ∫∫ 1/MM' dS(M') over a face
    ! and its derivative VS0 with respect to M.

    ! Based on formulas A6.1 and A6.3 (p. 381 to 383)
    ! in G. Delhommeau thesis (referenced below as [Del]).

    ! Inputs
    real(kind=pre), dimension(3), intent(in) :: field_point
    type(face_type), intent(in) :: face
    logical,                         intent(in) :: derivative_with_respect_to_first_variable

    ! Outputs
    real(kind=pre),               intent(out) :: int_G
    real(kind=pre), dimension(3), intent(out) :: int_nabla_G

    ! Local variables
    integer                         :: L
    real(kind=pre)                  :: GZ, DK, GY
    real(kind=pre), dimension(4)    :: RR
    real(kind=pre), dimension(3, 4) :: DRX
    real(kind=pre)                  :: ANT, DNT, ANL, DNL, ALDEN, AT
    real(kind=pre), dimension(3)    :: PJ, GYX, ANTX, ANLX, DNTX

    GZ = dot_product(field_point - face%center, face%normal) ! Called Z in [Del]

    do concurrent (L = 1:4)
      RR(L) = norm2(field_point - face%vertices(L, :))       ! Distance from vertices of face to field point
      DRX(:, L) = (field_point - face%vertices(L, :))/RR(L)  ! Normed vector from vertices of face to field point
    end do

    int_G = ZERO
    int_nabla_G(:) = ZERO

    do L = 1, 4
      DK = norm2(face%vertices(next_node(L), :) - face%vertices(L, :))    ! Distance between two consecutive points, called d_k in [Del]
      if (DK >= real(1e-3, PRE)*face%radius) then
        PJ(:) = (face%vertices(next_node(L), :) - face%vertices(L, :))/DK ! Normed vector from one corner to the next
        ! The following GYX(1:3) are called (a,b,c) in [Del]
        GYX(1) = face%normal(2)*PJ(3) - face%normal(3)*PJ(2)
        GYX(2) = face%normal(3)*PJ(1) - face%normal(1)*PJ(3)
        GYX(3) = face%normal(1)*PJ(2) - face%normal(2)*PJ(1)
        GY = dot_product(field_point - face%vertices(L, :), GYX)                       ! Called Y_k in  [Del]

        ANT = 2*GY*DK                                                                  ! Called N^t_k in [Del]
        DNT = (RR(NEXT_NODE(L))+RR(L))**2 - DK*DK + 2*abs(GZ)*(RR(NEXT_NODE(L))+RR(L)) ! Called D^t_k in [Del]
        ANL = RR(NEXT_NODE(L)) + RR(L) + DK                                            ! Called N^l_k in [Del]
        DNL = RR(NEXT_NODE(L)) + RR(L) - DK                                            ! Called D^l_k in [Del]
        ALDEN = log(ANL/DNL)

        if (abs(GZ) >= real(1e-4, PRE)*face%radius) then
          AT = atan(ANT/DNT)
        else
          AT = 0.
        endif

        ANLX(:) = DRX(:, NEXT_NODE(L)) + DRX(:, L)                    ! Called N^l_k_{x,y,z} in [Del]

        ANTX(:) = 2*DK*GYX(:)                                         ! Called N^t_k_{x,y,z} in [Del]
        DNTX(:) = 2*(RR(NEXT_NODE(L)) + RR(L) + abs(GZ))*ANLX(:)      &
          + 2*sign(ONE, GZ)*(RR(NEXT_NODE(L)) + RR(L))*face%normal(:) ! Called D^t_k_{x,y,z} in [Del]

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
          - 2*sign(ONE, GZ)*AT*face%normal(:)          &
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
      (field_point, face,                       &
      derivative_with_respect_to_first_variable, &
      int_G, int_nabla_G)
    ! Same as above, but always use the approximate asymptotic value.

    ! Inputs
    real(kind=pre), dimension(3), intent(in) :: field_point
    type(face_type), intent(in) :: face
    logical,    intent(in) :: derivative_with_respect_to_first_variable

    ! Outputs
    real(kind=pre),               intent(out) :: int_G
    real(kind=pre), dimension(3), intent(out) :: int_nabla_G

    ! Local variables
    real(kind=pre) :: ro
    real(kind=pre), dimension(3) :: vo

    vo = face%center - field_point
    ro = norm2(vo) ! distance from center of mass of the face to field point

    if (ro > real(1e-7, kind=pre)) then
      ! Asymptotic value if face far away from M
      int_G            = face%area/ro
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
      (field_point, face,                       &
      derivative_with_respect_to_first_variable, &
      reflection_coefs,                          &
      int_G, int_nabla_G)

    real(kind=pre), dimension(3), intent(in) :: field_point
    type(face_type), intent(in) :: face
    logical,                         intent(in) :: derivative_with_respect_to_first_variable
    real(kind=pre), dimension(2),    intent(in) :: reflection_coefs

    ! Outputs
    real(kind=pre),               intent(out) :: int_G
    real(kind=pre), dimension(3), intent(out) :: int_nabla_G

    ! Local
    real(kind=pre) :: r1

    r1 = hypot(norm2(field_point(1:2) - face%center(1:2)),  &
               reflection_coefs(1) * field_point(3) + reflection_coefs(2) - face%center(3))

    if (r1 > 7*face%radius) then
      call one_point_integral_of_reflected_Rankine &
        (field_point, face,                        &
        derivative_with_respect_to_first_variable, &
        reflection_coefs,                          &
        int_G, int_nabla_G)
    else
      call exact_integral_of_reflected_Rankine     &
        (field_point, face,                        &
        derivative_with_respect_to_first_variable, &
        reflection_coefs,                          &
        int_G, int_nabla_G)
    end if
  end subroutine integral_of_reflected_Rankine

  ! =================================================================

  pure subroutine exact_integral_of_reflected_Rankine &
      (field_point, face,                            &
      derivative_with_respect_to_first_variable,     &
      reflection_coefs,                              &
      int_g, int_nabla_g)
    ! Reflecting with respect to the free surface

    ! Inputs
    real(kind=pre), dimension(3), intent(in) :: field_point
    type(face_type), intent(in) :: face
    logical,                      intent(in) :: derivative_with_respect_to_first_variable
    real(kind=pre), dimension(2), intent(in) :: reflection_coefs

    ! Outputs
    real(kind=pre),               intent(out) :: int_g
    real(kind=pre), dimension(3), intent(out) :: int_nabla_g

    ! local
    real(kind=pre), dimension(3) :: reflected_point

    ! Create reflected point
    reflected_point(1:2) = field_point(1:2)
    reflected_point(3) = reflection_coefs(1) * field_point(3) + reflection_coefs(2)

    ! For instance
    ! reflection_coefs = [-1.0, 0.0] for free surface symmetry
    ! reflection_coefs = [-1.0, -2*h] for sea bottom symmetry

    call exact_integral_of_Rankine &
      (reflected_point, face,      &
      .true.,                      &
      int_g, int_nabla_g)

    ! Adjust the gradient components for the reflected point
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
      (field_point, face,                                &
      derivative_with_respect_to_first_variable,         &
      reflection_coefs,                                  &
      int_G, int_nabla_G)

    ! Inputs
    real(kind=pre), dimension(3), intent(in) :: field_point
    type(face_type), intent(in) :: face
    logical,                      intent(in) :: derivative_with_respect_to_first_variable
    real(kind=pre), dimension(2), intent(in) :: reflection_coefs

    ! outputs
    real(kind=pre),               intent(out) :: int_G
    real(kind=pre), dimension(3), intent(out) :: int_nabla_G

    ! local
    real(kind=pre), dimension(3) :: reflected_point

    ! Create reflected point
    reflected_point(1:2) = field_point(1:2)
    reflected_point(3) = reflection_coefs(1) * field_point(3) + reflection_coefs(2)

    call one_point_integral_of_Rankine(    &
      reflected_point, face,               &
      .true.,                             &
      int_G, int_nabla_G)

    ! Adjust the gradient components based on reflection
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
