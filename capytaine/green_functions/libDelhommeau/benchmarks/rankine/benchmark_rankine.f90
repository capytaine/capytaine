program benchmark_rankine

use constants, only: pre  ! Floating point precision
use constants, only: nb_tabulated_values
use green_rankine

implicit none

integer, parameter :: n_samples = 1000
integer :: i_sample

integer(kind=8) :: starting_time, final_time, clock_rate, clock_rate_in_ns
real(kind=pre), dimension(:, :), allocatable  :: points
real(kind=pre), dimension(4, 3) :: face_nodes
real(kind=pre), dimension(3) :: face_center
real(kind=pre), dimension(3) :: face_normal
real(kind=pre) :: face_area
real(kind=pre) :: face_radius

real(kind=pre), dimension(:), allocatable :: S
real(kind=pre), dimension(:, :), allocatable :: VS

print*, "-- Run libDelhommeau/benchmarks/rankine/benchmark_rankine.f90"

call system_clock(count_rate=clock_rate)
clock_rate_in_ns = clock_rate / 1000000000

allocate(points(3, n_samples))
call random_number(points)

face_nodes(1, :) = [0d0, 0d0, 0d0]
face_nodes(2, :) = [1d0, 0d0, 0d0]
face_nodes(3, :) = [1d0, 1d0, 0d0]
face_nodes(4, :) = [0d0, 1d0, 0d0]
face_center(:) = [0.5d0, 0.5d0, 0d0]
face_normal(:) = [0d0, 0d0, 1d0]
face_area = 1d0
face_radius = sqrt(2d0)/2

allocate(S(n_samples))
allocate(VS(3, n_samples))


call system_clock(starting_time)
do i_sample = 1, n_samples
    call COMPUTE_INTEGRAL_OF_RANKINE_SOURCE(                          &
        points(:, i_sample),                                          &
        face_nodes, face_center, face_normal, face_area, face_radius, &
        S(i_sample), VS(:, i_sample)                                  &
    )
enddo
call system_clock(final_time)
print*, "Rankine             :", (final_time - starting_time)/clock_rate_in_ns/n_samples, " ns"


call system_clock(starting_time)
do i_sample = 1, n_samples
    call COMPUTE_ASYMPTOTIC_RANKINE_SOURCE( &
        points(:, i_sample),                &
        face_nodes, face_area,              &
        S(i_sample), VS(:, i_sample)        &
    )
enddo
call system_clock(final_time)
print*, "Asymptotic Rankine  :", (final_time - starting_time)/clock_rate_in_ns/n_samples, " ns"

end program
