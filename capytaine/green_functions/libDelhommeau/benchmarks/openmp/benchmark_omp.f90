program benchmark_omp

use omp_lib
use ieee_arithmetic

use matrices, only: build_matrices
use delhommeau_integrals, only: default_r_spacing, default_z_spacing, construct_tabulation
use constants, only: pre  ! Floating point precision
use constants, only: nb_tabulated_values
use old_prony_decomposition, only: lisc

implicit none

integer(kind=8) :: starting_time, final_time, clock_rate

integer, parameter :: nb_faces = 8192
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
integer, parameter :: tabulation_nr = 328
integer, parameter :: tabulation_nz = 46
real(kind=pre), dimension(tabulation_nr)                       :: tabulated_r
real(kind=pre), dimension(tabulation_nz)                       :: tabulated_z
real(kind=pre), dimension(tabulation_nr, tabulation_nz, nb_tabulated_values)  :: tabulated_integrals

integer, parameter :: gf_singularities = 1  ! low_freq

integer, parameter :: finite_depth_method = 1
! Prony decomposition for the finite depth Green function
integer :: nexp
real(kind=pre), dimension(2, 31) :: prony_decomposition
real(kind=pre), dimension(1) :: dispersion_roots  ! Dummy, only for FinGreen3D

integer n_threads

! The interaction matrices to be computed
complex(kind=pre), dimension(:, :), allocatable :: S
complex(kind=pre), dimension(:, :, :), allocatable :: K

print*, "-- Run libdelhommeau/benchmark/openmp/benchmark_omp.f90"

call RANDOM_INIT(.true.,.true.)

allocate(S(nb_faces, nb_faces))
allocate(K(nb_faces, nb_faces, 1))

tabulated_r(:) = default_r_spacing(tabulation_nr, 100d0, tabulation_grid_shape)
tabulated_z(:) = default_z_spacing(tabulation_nz, -251d0, tabulation_grid_shape)
tabulated_integrals(:, :, :) = construct_tabulation(tabulated_r, tabulated_z, tabulation_nb_integration_points)

wavenumber = 1.0

if (.true.) then
   depth = ieee_value(depth, ieee_positive_inf)
else
   depth = 50.
   call lisc(real(wavenumber*depth*tanh(wavenumber*depth)), real(wavenumber*depth), nexp, prony_decomposition)
end if

call random_panels(nb_faces, vertices, faces, face_center, face_normal, face_area, face_radius)
quadrature_points = reshape(face_center, shape(quadrature_points))
quadrature_weights = reshape(face_area, shape(quadrature_weights))

! ! For debugging
! open (unit=4, file='vertices.dat', form='formatted')
! do i = 1, 4*nb_faces
!    write (4, *) vertices(i, :)
! end do

call system_clock(count_rate=clock_rate)

open(unit=210, file="benchmark_omp.csv")
print'(a)', "           kind num. thread(s) elapsed time (s)"
write(210, '(a)') "n_threads,elapsed_time,kind"

do n_threads = 1, OMP_GET_MAX_THREADS()
  call omp_set_num_threads(n_threads)

  call system_clock(starting_time)
  call build_matrices(                                           &
    nb_faces, face_center, face_normal,                          &
    nb_vertices, nb_faces, vertices, faces,                      &
    face_center, face_normal, face_area, face_radius,            &
    nb_quadrature_points, quadrature_points, quadrature_weights, &
    wavenumber, depth,                                           &
    tabulation_nb_integration_points, tabulation_grid_shape,     &
    tabulated_r, tabulated_z, tabulated_integrals,               &
    finite_depth_method, prony_decomposition, dispersion_roots,  &
    gf_singularities, .true.,                                    &
    S, K)
  call system_clock(final_time)

  print'(a15,1i15,a,1ES16.6)', " full", n_threads, ' ', real(final_time - starting_time)/clock_rate
  write(210,'(1i3,a,1ES12.6,a)') n_threads, ",", real(final_time - starting_time)/clock_rate, ",full"
enddo

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

end program benchmark_omp
