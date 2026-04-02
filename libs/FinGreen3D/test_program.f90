! Just test that the library compiles and run

program test
  use fingreen3d_module, only: fingreen3d_routine, dispersion

  implicit none

  integer, parameter :: n = 10
  integer, parameter :: nk = 200
  real(8), parameter :: g = 9.81

  real(8) :: omega, water_depth

  real(8), dimension(nk) :: roots_of_dispersion_relationship
  real(8), dimension(n) :: r, x3, xi3
  complex(8), dimension(n, 3) :: gf
  integer :: i

  omega = 1d0
  water_depth = 10d0

  call dispersion(roots_of_dispersion_relationship, nk, omega, water_depth)

  do i = 1, n
    r(i) = real(i, kind=8)/n
    x3(i) = -real(i, kind=8)/n
    xi3(i) = -2*real(i, kind=8)/n
  end do

  !$OMP PARALLEL DO PRIVATE(i)
  do i = 1, n
    call fingreen3d_routine(r(i), x3(i), xi3(i), sqrt(g*omega), roots_of_dispersion_relationship, nk, water_depth, gf(i, :))
  end do

  print*, real(gf(:, 1))
  print*, real(gf(:, 2))
  print*, real(gf(:, 3))

end program test
