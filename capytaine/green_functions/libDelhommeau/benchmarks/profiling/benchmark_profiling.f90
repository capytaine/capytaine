program benchmark_profiling

use omp_lib
use ieee_arithmetic

use matrices, only: build_matrices
use delhommeau_integrals, only: default_r_spacing, default_z_spacing, construct_tabulation
use constants, only: pre  ! Floating point precision
use constants, only: nb_tabulated_values
use old_prony_decomposition, only: lisc

implicit none

integer(kind=8) :: starting_time, final_time, clock_rate

integer, parameter :: nb_faces = 4096
integer, parameter :: nb_vertices = 4*nb_faces
integer, parameter :: nb_quadrature_points = 1

real(kind=pre) :: wavenumber, depth

! Geometry of the mesh
real(kind=pre), dimension(nb_vertices, 3) :: vertices
integer, dimension(nb_faces, 4) :: faces
real(kind=pre), dimension(nb_faces, 3) :: face_center
real(kind=pre), dimension(nb_faces, 3) :: face_normal
real(kind=pre), dimension(nb_faces) :: face_area
real(kind=pre), dimension(nb_faces) :: face_radius
real(kind=pre), dimension(nb_faces, nb_quadrature_points, 3) :: quadrature_points
real(kind=pre), dimension(nb_faces, nb_quadrature_points) :: quadrature_weights

! Tabulation of the integrals used in the Green function
integer, parameter :: tabulation_grid_shape = 1
integer, parameter :: tabulation_nb_integration_points = 1000
integer, parameter :: tabulation_nr = 676
integer, parameter :: tabulation_nz = 372
real(kind=pre), dimension(tabulation_nr)                       :: tabulated_r
real(kind=pre), dimension(tabulation_nz)                       :: tabulated_z
real(kind=pre), dimension(tabulation_nr, tabulation_nz, nb_tabulated_values)  :: tabulated_integrals

integer, parameter :: gf_singularities = 0

! Prony decomposition for the finite depth Green function
integer, parameter :: nexp_max = 31
integer :: nexp
real, dimension(nexp_max) :: ambda_f32, ar_f32
real(kind=pre), dimension(nexp_max) :: ambda, ar

integer i
real(kind=pre), dimension(3) :: coeffs

! The interaction matrices to be computed
complex(kind=pre), dimension(:, :), allocatable :: S
complex(kind=pre), dimension(:, :, :), allocatable :: K

print*, "-- Run libdelhommeau/benchmark/profiling/profiling.f90"

call RANDOM_INIT(.true.,.true.)

allocate(S(nb_faces, nb_faces))
allocate(K(nb_faces, nb_faces, 1))

tabulated_r(:) = default_r_spacing(tabulation_nr, 100d0, tabulation_grid_shape)
tabulated_z(:) = default_z_spacing(tabulation_nz, -251d0, tabulation_grid_shape)
tabulated_integrals(:, :, :) = construct_tabulation(tabulated_r, tabulated_z, tabulation_nb_integration_points)

wavenumber = 1.0

call random_panels(nb_faces, vertices, faces, face_center, face_normal, face_area, face_radius)
quadrature_points = reshape(face_center, shape(quadrature_points))
quadrature_weights = reshape(face_area, shape(quadrature_weights))

print'(a)', "depth kind           elapsed time (s)"

do i=1, 1

   if (i .eq. 1) then
      depth = ieee_value(depth, ieee_positive_inf)
   else
      depth = 50.
      call lisc(real(wavenumber*depth*tanh(wavenumber*depth)), real(wavenumber*depth), ambda_f32, ar_f32, nexp)
      ambda(:) = real(ambda_f32(:), kind=pre)
      ar(:) = real(ar_f32(:), kind=pre)
      nexp = nexp + 1
      ambda(nexp) = 0.0
      ar(nexp) = 2.0
   end if

   call system_clock(count_rate=clock_rate)

   coeffs = [1d0, 1d0, 1d0]
   call system_clock(starting_time)
   call build_matrices(                                              &
        nb_faces, face_center, face_normal,                          &
        nb_vertices, nb_faces, vertices, faces,                      &
        face_center, face_normal, face_area, face_radius,            &
        nb_quadrature_points, quadrature_points, quadrature_weights, &
        wavenumber, depth,                                           &
        coeffs,                                                      &
        tabulation_nb_integration_points, tabulation_grid_shape,     &
        tabulated_r, tabulated_z, tabulated_integrals,               &
        nexp, ambda, ar,                                             &
        .false., gf_singularities, .true.,                           &
        S, K)
   call system_clock(final_time)

   print'(1F5.0,a16,1ES16.6)', depth, " full ", real(final_time - starting_time)/clock_rate

   ! coeffs = [0d0, 0d0, 1d0]
   ! call system_clock(starting_time)
   ! call build_matrices(                                                   &
   !      nb_faces, face_center, face_normal,                               &
   !      nb_vertices, nb_faces, vertices, faces,                           &
   !      face_center, face_normal, face_area, face_radius,                 &
   !      nb_quadrature_points, quadrature_points, quadrature_weights,      &
   !      wavenumber, depth,                                                &
   !      coeffs,                                                           &
   !      tabulation_nb_integration_points, tabulation_grid_shape,          &
   !      tabulated_r, tabulated_z, tabulated_integrals,                    &
   !      nexp, ambda, ar,                                                  &
   !      .false., .true.,                                                  &
   !      S, K)
   ! call system_clock(final_time)
   !
   ! ! print'(1F5.0,a16,1ES16.6)', depth, " wave_only ", real(final_time - starting_time)/clock_rate
   !
   ! coeffs = [0d0, 0d0, 1d0]
   ! call system_clock(starting_time)
   ! call build_matrices(                                                       &
   !      nb_faces, face_center, face_normal,                                   &
   !      nb_vertices, nb_faces, vertices, faces,                               &
   !      face_center, face_normal, face_area, face_radius,                     &
   !      nb_quadrature_points, quadrature_points, quadrature_weights,          &
   !      wavenumber, depth,                                                    &
   !      coeffs,                                                               &
   !      tabulation_nb_integration_points, tabulation_grid_shape,              &
   !      tabulated_r, tabulated_z, tabulated_integrals,                        &
   !      nexp, ambda, ar,                                                      &
   !      .true., gf_singularities, .true.,                                     &
   !      S, K)
   ! call system_clock(final_time)
   !
   ! print'(1F5.0,a16,1ES16.6)', depth, " half_wave_only ", real(final_time - starting_time)/clock_rate

end do

contains

  subroutine random_panels(nb_faces, vertices, faces, face_center, face_normal, face_area, face_radius)
    integer, intent(in) :: nb_faces

    integer, dimension(nb_faces, 4), intent(out) :: faces
    real(kind=pre), dimension(4*nb_faces, 3), intent(out) :: vertices
    real(kind=pre), dimension(nb_faces, 3), intent(out) :: face_center
    real(kind=pre), dimension(nb_faces, 3), intent(out) :: face_normal
    real(kind=pre), dimension(nb_faces), intent(out) :: face_area
    real(kind=pre), dimension(nb_faces), intent(out) :: face_radius

    real(kind=pre), dimension(3, 2) :: vertex_shifts
    integer :: i

    faces = transpose(reshape([(i, i=1,(nb_faces*4),1)], [4, nb_faces]))

    call random_number(face_area)
    face_area = 1.0 + face_area

    face_radius = face_area * sqrt(2.0)/2.0

    call random_number(face_center)
    face_center(:, 1) = 20*(face_center(:, 1) - 0.5)
    face_center(:, 2) = 20*(face_center(:, 2) - 0.5)
    face_center(:, 3) = -10*face_center(:, 3) - face_radius(:)

    call random_number(face_normal)
    face_normal(:, :) = face_normal(:, :) - 0.5
    do i = 1, nb_faces
      face_normal(i, :) = face_normal(i, :)/norm2(face_normal(i, :))
    enddo

    do i = 1, nb_faces
      vertex_shifts = face_radius(i) * two_orthogonal_vector(face_normal(i, :))
      vertices(4*(i-1)+1, :) = face_center(i, :) + vertex_shifts(:, 1)
      vertices(4*(i-1)+2, :) = face_center(i, :) + vertex_shifts(:, 2)
      vertices(4*(i-1)+3, :) = face_center(i, :) - vertex_shifts(:, 1)
      vertices(4*(i-1)+4, :) = face_center(i, :) - vertex_shifts(:, 2)
    enddo
  end subroutine

  pure function two_orthogonal_vector(n) result(vecs)
    ! Given a normal vector `n`, returns two other vectors
    ! such that the three of them is an orthonormal basis
    real(kind=pre), dimension(3), intent(in) :: n
    real(kind=pre), dimension(3, 2) :: vecs

    vecs(:, 1) = [n(2), -n(1), 0d0]
    if (norm2(vecs(:, 1)) < abs(1e-5)) then
      vecs(:, 1) = [n(3), 0d0, -n(1)]
    endif

    ! Cross product
    vecs(1, 2) = n(2) * vecs(3, 1) - n(3) * vecs(2, 1)
    vecs(2, 2) = n(3) * vecs(1, 1) - n(1) * vecs(3, 1)
    vecs(3, 2) = n(1) * vecs(2, 1) - n(2) * vecs(1, 1)

    vecs(:, 1) = vecs(:, 1)/norm2(vecs(:, 1))
    vecs(:, 2) = vecs(:, 2)/norm2(vecs(:, 2))
  end function

end program benchmark_profiling
