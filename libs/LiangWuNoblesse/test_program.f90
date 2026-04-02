! Just test that the library compiles and run

program test
  use LiangWuNoblesseWaveTerm, only: HavelockGF
  implicit none

  integer, parameter :: n = 10

  real(8), dimension(n) :: r, z
  complex(8), dimension(n) :: gf, gf_r
  integer :: i

  do i = 1, n
    r(i) = real(i, kind=8)/n
  end do

  do i = 1, n
    z(i) = -real(i, kind=8)/n
  end do

  !$OMP PARALLEL DO PRIVATE(i)
  do i = 1, n
    call HavelockGF(r(i), z(i), gf(i), gf_r(i))
  end do

  print*, real(gf(:))

end program test
