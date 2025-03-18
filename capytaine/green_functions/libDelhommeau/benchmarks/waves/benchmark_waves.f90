program benchmark_waves

use constants
use ieee_arithmetic
use old_prony_decomposition, only: lisc
use delhommeau_integrals
use green_wave

implicit none

integer, parameter :: n_samples = 1000
integer :: i_sample

integer(kind=8) :: starting_time, final_time, clock_rate, clock_rate_in_ns

real(kind=pre) :: wavenumber, depth

integer, parameter :: tabulation_nb_integration_points = 1001
integer, parameter :: tabulation_grid_shape = SCALED_NEMOH3_GRID
integer, parameter :: tabulation_nr = 676
integer, parameter :: tabulation_nz = 372
real(kind=pre), dimension(tabulation_nr) :: tabulated_r_range
real(kind=pre), dimension(tabulation_nz) :: tabulated_z_range
real(kind=pre), dimension(tabulation_nr, tabulation_nz, nb_tabulated_values) :: tabulated_integrals

integer, parameter :: gf_singularities = LOW_FREQ

integer, parameter :: finite_depth_method = NEWER_FINITE_DEPTH
integer :: nexp
real(kind=pre), dimension(2, 31) :: prony_decomposition
real(kind=pre), dimension(1) :: dispersion_roots  ! Dummy

real(kind=pre), dimension(:, :), allocatable  :: points
real(kind=pre), dimension(4, 3) :: face_nodes
real(kind=pre), dimension(1, 3) :: face_quadrature_points
real(kind=pre), dimension(1) :: face_quadrature_weights
real(kind=pre), dimension(3) :: face_center
real(kind=pre), dimension(3) :: face_normal
real(kind=pre) :: face_area
real(kind=pre) :: face_radius
logical, parameter :: derivative_with_respect_to_first_variable = .true.

complex(kind=pre), dimension(:), allocatable :: S
complex(kind=pre), dimension(:, :), allocatable :: VS

print*, "-- Run libDelhommeau/benchmarks/waves/benchmark_waves.f90"

dispersion_roots(:) = 0.0

call system_clock(count_rate=clock_rate)
clock_rate_in_ns = clock_rate / 1000000000

allocate(points(3, n_samples))
call random_number(points)
points(3, :) = - points(3, :)

face_nodes(1, :) = [0d0, 0d0, 0d0]
face_nodes(2, :) = [1d0, 0d0, 0d0]
face_nodes(3, :) = [1d0, 1d0, 0d0]
face_nodes(4, :) = [0d0, 1d0, 0d0]
face_center(:) = [0.5d0, 0.5d0, 0d0]
face_quadrature_points(1, :) = [0.5d0, 0.5d0, 0d0]
face_quadrature_weights(1) = 1d0
face_normal(:) = [0d0, 0d0, 1d0]
face_area = 1d0
face_radius = sqrt(2d0)/2

tabulated_r_range(:) = default_r_spacing(tabulation_nr, 100d0, tabulation_grid_shape)
tabulated_z_range(:) = default_z_spacing(tabulation_nz, -251d0, tabulation_grid_shape)
!tabulated_integrals = construct_tabulation(tabulated_r_range, tabulated_z_range, tabulation_nb_integration_points)

allocate(S(n_samples))
allocate(VS(3, n_samples))

depth = 10d0
wavenumber = 1d0
call lisc(real(wavenumber*depth*tanh(wavenumber*depth)), real(wavenumber*depth), nexp, prony_decomposition)
call system_clock(starting_time)
do i_sample = 1, n_samples
!  call integral_of_wave_part                                     &
!      (points(:, i_sample),                                      &
!      face_nodes,                                                &
!      face_center, face_normal, face_area, face_radius,          &
!      face_quadrature_points, face_quadrature_weights,           &
!      wavenumber, depth,                                         &
!      tabulation_nb_integration_points, tabulation_grid_shape,   &
!      tabulated_r_range, tabulated_z_range, tabulated_integrals, &
!      gf_singularities,                                          &
!      finite_depth_method, prony_decomposition, dispersion_roots,&
!      derivative_with_respect_to_first_variable,                 &
!      S(i_sample), VS(:, i_sample)                               &
!      )
enddo
call system_clock(final_time)
print*, "Infinite depth :", (final_time - starting_time)/clock_rate_in_ns/n_samples, " ns"

end program
