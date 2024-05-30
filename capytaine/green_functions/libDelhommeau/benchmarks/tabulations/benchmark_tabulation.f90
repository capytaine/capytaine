program benchmark_tabulation

use constants, only: pre  ! Floating point precision
use constants, only: nb_tabulated_values
use delhommeau_integrals

implicit none

integer, parameter :: tabulation_grid_shape = 1
integer, parameter :: n_tabulation = 4
integer, parameter :: n_samples = 1000
integer :: i_sample
integer :: i_tabulation
integer, dimension(n_tabulation) :: tabulation_nr
integer, dimension(n_tabulation) :: tabulation_nz


integer(kind=8) :: starting_time, final_time, clock_rate, clock_rate_in_ns

real(kind=pre), dimension(:), allocatable          :: r, z, tabulated_r, tabulated_z
real(kind=pre), dimension(:, :, :), allocatable    :: tabulated_integrals
real(kind=pre), dimension(:, :), allocatable       :: integrals

print*, "-- Run libdelhommeau/benchmark/tabulations/benchmark_tabulation.f90"

call system_clock(count_rate=clock_rate)
clock_rate_in_ns = clock_rate / 1000000000

allocate(r(n_samples))
allocate(z(n_samples))
allocate(integrals(n_samples, nb_tabulated_values))

call random_number(r)
r = 100*r

call random_number(z)
z = -10*z


call system_clock(starting_time)
do i_sample = 1, n_samples
    integrals(i_sample, :) = numerical_integration(r(i_sample), z(i_sample), 2501)
enddo
call system_clock(final_time)
print*, "Integration (fine)  :", (final_time - starting_time)/clock_rate_in_ns/n_samples, " ns"


call system_clock(starting_time)
do i_sample = 1, n_samples
    integrals(i_sample, :) = numerical_integration(r(i_sample), z(i_sample), 251)
enddo
call system_clock(final_time)
print*, "Integration (coarse):", (final_time - starting_time)/clock_rate_in_ns/n_samples, " ns"


tabulation_nr = [32, 328, 656, 656]
tabulation_nz = [4, 48, 92, 184]

do i_tabulation = 1, n_tabulation

  allocate(tabulated_r(tabulation_nr(i_tabulation)))
  allocate(tabulated_z(tabulation_nz(i_tabulation)))
  allocate(tabulated_integrals(tabulation_nr(i_tabulation), tabulation_nz(i_tabulation), nb_tabulated_values))

  tabulated_r(:) = default_r_spacing(tabulation_nr(i_tabulation), 100d0, tabulation_grid_shape)
  tabulated_z(:) = default_z_spacing(tabulation_nz(i_tabulation), -251d0, tabulation_grid_shape)
  tabulated_integrals(:, :, :) = construct_tabulation(tabulated_r, tabulated_z, 1001)

  call system_clock(starting_time)
  do i_sample = 1, n_samples
      integrals(i_sample, :) = pick_in_default_tabulation(r(i_sample), z(i_sample), &
        tabulation_grid_shape, tabulated_r, tabulated_z, tabulated_integrals)
  enddo
  call system_clock(final_time)
  print*, "Interpolate         :", (final_time - starting_time)/clock_rate_in_ns/n_samples, " ns"

  deallocate(tabulated_r, tabulated_z, tabulated_integrals)

enddo


call system_clock(starting_time)
do i_sample = 1, n_samples
    integrals(i_sample, :) = asymptotic_approximations(r(i_sample), z(i_sample))
enddo
call system_clock(final_time)
print*, "Asymptotic          :", (final_time - starting_time)/clock_rate_in_ns/n_samples, " ns"


end program
