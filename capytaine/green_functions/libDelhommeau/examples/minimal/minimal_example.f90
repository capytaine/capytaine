program test

  use ieee_arithmetic

  use matrices, only: build_matrices, is_infinity
  use delhommeau_integrals, only: default_r_spacing, default_z_spacing, construct_tabulation
  use floating_point_precision, only: pre
  use constants, only: zero, nb_tabulated_values
  use old_prony_decomposition, only: lisc

  implicit none

  integer, parameter :: nb_faces = 2
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
  integer, parameter :: tabulation_grid_shape = 1   ! scaled_nemoh3 method
  integer, parameter :: tabulation_nb_integration_points = 251
  integer, parameter :: tabulation_nr = 676
  integer, parameter :: tabulation_nz = 372
  real(kind=pre), dimension(tabulation_nr)         :: tabulated_r
  real(kind=pre), dimension(tabulation_nz)         :: tabulated_z
  real(kind=pre), allocatable, dimension(:, :, :)  :: tabulated_integrals

  integer, parameter :: gf_singularities = 0  ! high_freq

  ! Prony decomposition for the finite depth Green function
  integer, parameter :: nexp_max = 31
  integer :: nexp
  real, dimension(nexp_max) :: ambda_f32, ar_f32
  real(kind=pre), dimension(nexp_max) :: ambda, ar

  ! The interaction matrices to be computed
  complex(kind=pre), dimension(nb_faces, nb_faces) :: S
  complex(kind=pre), dimension(nb_faces, nb_faces, 1) :: K

  integer :: i
  real(kind=pre), dimension(3) :: coeffs

  tabulated_r(:) = default_r_spacing(tabulation_nr, 100d0, tabulation_grid_shape)
  tabulated_z(:) = default_z_spacing(tabulation_nz, -251d0, tabulation_grid_shape)
  allocate(tabulated_integrals(tabulation_nr, tabulation_nz, nb_tabulated_values))
  tabulated_integrals = construct_tabulation(tabulated_r, tabulated_z, tabulation_nb_integration_points)

  depth = ieee_value(depth, ieee_positive_inf)

  vertices = reshape([  &
    0.0,  1.0,  1.0,  0.0,  1.0,  2.0,  2.0,  1.0,  &
    0.0,  0.0,  1.0,  1.0,  0.0,  0.0,  1.0,  1.0,  &
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0  &
    ], shape(vertices))

  faces = reshape([  &
    1, 5, &
    2, 6, &
    3, 7, &
    4, 8  &
    ], shape(faces))

  face_center = reshape([  &
    0.5,  1.5,  &
    0.5,  0.5,  &
    -1.0, -1.0  &
    ], shape(face_center))

  face_normal = reshape([  &
    0.0, 0.0,  &
    0.0, 0.0,  &
    1.0, 1.0   &
    ], shape(face_normal))

  face_area = [1.0, 1.0]
  face_radius = [0.71, 0.71]

  quadrature_points = reshape(face_center, shape(quadrature_points))
  quadrature_weights = reshape(face_area, shape(quadrature_weights))

  print*, "-- Run libdelhommeau/examples/minimal/minimal_example.f90"

  coeffs = [1d0, 0d0, 0d0]
  call build_matrices(                                           &
    nb_faces, face_center, face_normal,                          &
    nb_vertices, nb_faces, vertices, faces,                      &
    face_center, face_normal, face_area, face_radius,            &
    nb_quadrature_points, quadrature_points, quadrature_weights, &
    ZERO, depth, coeffs,                                         &
    tabulation_nb_integration_points, tabulation_grid_shape,     &
    tabulated_r, tabulated_z, tabulated_integrals,               &
    nexp, ambda, ar,                                             &
    .true., gf_singularities, .true.,                            &
    S, K)
  print*, "Rankine part: S"
  do i = 1, nb_faces
    print"(4ES20.12)", S(i, :)
  enddo
  print*, "Rankine part: K"
  do i = 1, nb_faces
    print"(4ES20.12)", K(i, :, 1)
  enddo

  wavenumber = 1d0
  coeffs = [1d0, 1d0, 1d0]
  call build_matrices(                                           &
    nb_faces, face_center, face_normal,                          &
    nb_vertices, nb_faces, vertices, faces,                      &
    face_center, face_normal, face_area, face_radius,            &
    nb_quadrature_points, quadrature_points, quadrature_weights, &
    wavenumber, depth, coeffs,                                   &
    tabulation_nb_integration_points, tabulation_grid_shape,     &
    tabulated_r, tabulated_z, tabulated_integrals,               &
    nexp, ambda, ar,                                             &
    .true., gf_singularities, .true.,                            &
    S, K)
  print*, "k=1.0, h=infty: S"
  do i = 1, nb_faces
    print"(4ES20.12)", S(i, :)
  enddo
  print*, "k=1.0, h=infty: K"
  do i = 1, nb_faces
    print"(4ES20.12)", K(i, :, 1)
  enddo

  wavenumber = 2d0
  coeffs = [1d0, 1d0, 1d0]
  call build_matrices(                                           &
    nb_faces, face_center, face_normal,                          &
    nb_vertices, nb_faces, vertices, faces,                      &
    face_center, face_normal, face_area, face_radius,            &
    nb_quadrature_points, quadrature_points, quadrature_weights, &
    wavenumber, depth, coeffs,                                   &
    tabulation_nb_integration_points, tabulation_grid_shape,     &
    tabulated_r, tabulated_z, tabulated_integrals,               &
    nexp, ambda, ar,                                             &
    .true., gf_singularities, .true.,                            &
    S, K)
  print*, "k=2.0, h=infty: S"
  do i = 1, nb_faces
    print"(4ES20.12)", S(i, :)
  enddo
  print*, "k=2.0, h=infty: K"
  do i = 1, nb_faces
    print"(4ES20.12)", K(i, :, 1)
  enddo

  ! finite depth

  depth = 2.0
  wavenumber = 1d0
  coeffs = [1d0, 1d0, 1d0]

  call lisc(real(wavenumber*depth*tanh(wavenumber*depth)), real(wavenumber*depth), ambda_f32, ar_f32, nexp)
  ambda(:) = real(ambda_f32(:), kind=pre)
  ar(:) = real(ar_f32(:), kind=pre)
  nexp = nexp + 1
  ambda(nexp) = 0.0
  ar(nexp) = 2.0

  call build_matrices(                                           &
    nb_faces, face_center, face_normal,                          &
    nb_vertices, nb_faces, vertices, faces,                      &
    face_center, face_normal, face_area, face_radius,            &
    nb_quadrature_points, quadrature_points, quadrature_weights, &
    wavenumber, depth, coeffs,                                   &
    tabulation_nb_integration_points, tabulation_grid_shape,     &
    tabulated_r, tabulated_z, tabulated_integrals,               &
    nexp, ambda, ar,                                             &
    .true., gf_singularities, .true.,                            &
    S, K)
  print*, "k=1.0, h=2.0: S"
  do i = 1, nb_faces
    print"(4ES20.12)", S(i, :)
  enddo
  print*, "k=1.0, h=2.0: K"
  do i = 1, nb_faces
    print"(4ES20.12)", K(i, :, 1)
  enddo


  depth = 2.0
  wavenumber = 2d0
  coeffs = [1d0, 1d0, 1d0]

  call lisc(real(wavenumber*depth*tanh(wavenumber*depth)), real(wavenumber*depth), ambda_f32, ar_f32, nexp)
  ambda(:) = real(ambda_f32(:), kind=pre)
  ar(:) = real(ar_f32(:), kind=pre)
  nexp = nexp + 1
  ambda(nexp) = 0.0
  ar(nexp) = 2.0

  call build_matrices(                                           &
    nb_faces, face_center, face_normal,                          &
    nb_vertices, nb_faces, vertices, faces,                      &
    face_center, face_normal, face_area, face_radius,            &
    nb_quadrature_points, quadrature_points, quadrature_weights, &
    wavenumber, depth, coeffs,                                   &
    tabulation_nb_integration_points, tabulation_grid_shape,     &
    tabulated_r, tabulated_z, tabulated_integrals,               &
    nexp, ambda, ar,                                             &
    .true., gf_singularities, .true.,                            &
    S, K)
  print*, "k=2.0, h=2.0: S"
  do i = 1, nb_faces
    print"(4ES20.12)", S(i, :)
  enddo
  print*, "k=2.0, h=2.0: K"
  do i = 1, nb_faces
    print"(4ES20.12)", K(i, :, 1)
  enddo

end program test
