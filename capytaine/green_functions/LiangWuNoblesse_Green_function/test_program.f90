program test
  use GreenFuncMod, only: HavelockGF
  implicit none

  integer, parameter :: n = 10

  real(8), dimension(n) :: r
  complex(8), dimension(0:3, n) :: gf
  integer :: i

  do i = 1, n
    r(i) = real(i, kind=8)/n
  end do

  !$OMP PARALLEL DO PRIVATE(i)
  do i = 1, n
    call HavelockGF(r(i), 0d0, 1d0, gf(:, i))
  end do

  print*, real(gf(0, :))

end program test
