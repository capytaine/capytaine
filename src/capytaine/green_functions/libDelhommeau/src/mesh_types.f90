! Copyright (C) 2017-2025 Matthieu Ancellin
! See LICENSE file at <https://github.com/capytaine/libDelhommeau>

module mesh_types

  use floating_point_precision, only: pre
  use constants

  implicit none

  ! Custom type to represent a face/panel in the mesh
  type face_type
    real(kind=pre), dimension(3)   :: center      ! Center point of the face
    real(kind=pre), dimension(3)   :: normal      ! Normal vector of the face
    real(kind=pre), dimension(4,3) :: vertices    ! Coordinates of the 4 vertices
    real(kind=pre)                 :: area        ! Area of the face
    real(kind=pre)                 :: radius      ! Radius of the circumscribed circle
    real(kind=pre), dimension(:,:), allocatable :: quad_points  ! Integration points
    real(kind=pre), dimension(:), allocatable   :: quad_weights ! Integration weights
  end type face_type

  contains

  function create_face(face_vertices, face_center, face_normal, face_area, face_radius, &
                      face_quad_points, face_quad_weights) result(face)
    real(kind=pre), dimension(4,3), intent(in) :: face_vertices
    real(kind=pre), dimension(3), intent(in)   :: face_center
    real(kind=pre), dimension(3), intent(in)   :: face_normal
    real(kind=pre), intent(in)                 :: face_area
    real(kind=pre), intent(in)                 :: face_radius
    real(kind=pre), dimension(:,:), intent(in) :: face_quad_points
    real(kind=pre), dimension(:), intent(in)   :: face_quad_weights
    type(face_type) :: face

    face%vertices = face_vertices
    face%center = face_center
    face%normal = face_normal
    face%area = face_area
    face%radius = face_radius

    ! Allocate and copy quadrature data
    allocate(face%quad_points(size(face_quad_points,1), size(face_quad_points,2)))
    allocate(face%quad_weights(size(face_quad_weights)))
    face%quad_points = face_quad_points
    face%quad_weights = face_quad_weights
  end function create_face

  subroutine destroy_face(face)
    type(face_type), intent(inout) :: face
    if (allocated(face%quad_points)) deallocate(face%quad_points)
    if (allocated(face%quad_weights)) deallocate(face%quad_weights)
  end subroutine destroy_face

end module mesh_types
