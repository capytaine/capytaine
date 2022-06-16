! Copyright (C) 2022 Matthieu Ancellin
! See LICENSE file at <https://github.com/mancellin/capytaine>
!
! This module contains functions to evaluate the following integrals
! D1 = Re[ ∫(-i cosθ)(J(ζ) - 1/ζ)dθ ]
! D2 = Re[ ∫(-i cosθ)(e^ζ)dθ ]
! #ifdef XIE_CORRECTION
! Z1 = Re[ ∫(J(ζ))dθ ]
! #else
! Z1 = Re[ ∫(J(ζ) - 1/ζ)dθ ]
! #endif
! Z2 = Re[ ∫(e^ζ)dθ ]
! where ζ depends on θ, as well as two additional parameters `r ∈ [0, +∞)` and `z ∈ (-∞, 0]`.
!
! They are required for the evaluation of the Green function and its gradient.
!
module delhommeau_integrals

  implicit none

  integer, parameter :: pre = kind(1d0) ! floating point numbers precision
  real(kind=pre), parameter :: pi = 3.141592653588979 ! π
  complex(kind=pre), parameter :: ii = (0, 1)         ! imaginary unit

  public :: numerical_integration
  public :: asymptotic_approximations
  public :: construct_tabulation
  public :: default_r_spacing, default_z_spacing
  public :: pick_in_default_tabulation

contains

  pure function numerical_integration(r, z, nb_integration_points) result(integrals)
    ! Compute the integrals by numerical integration, with `nb_integration_points` points.

    ! input
    real(kind=pre), intent(in) :: r
    real(kind=pre), intent(in) :: z
    integer,        intent(in), optional :: nb_integration_points

    ! output
    real(kind=pre), dimension(2, 2) :: integrals

    ! local variables
    integer :: n, k
    real(kind=pre) :: theta, delta_theta, cos_theta
    complex(kind=pre) :: zeta, exp_zeta, jzeta

    if (present(nb_integration_points)) then
      n = nb_integration_points
    else
      n = 251
    endif

    delta_theta = pi/(n-1.0)  ! step of the integration

    ! initial values
    integrals(:, :) = 0.0
    theta = -pi/2

    do k = 1, n
      theta = theta + delta_theta
      cos_theta = cos(theta)
      zeta = cmplx(z, r*cos_theta, kind=pre)
      exp_zeta = exp(zeta)
      if (aimag(zeta) < 0) then
        jzeta = exp_e1(zeta) - ii*pi*exp_zeta
      else
        jzeta = exp_e1(zeta) + ii*pi*exp_zeta
      end if
      integrals(1, 1) = integrals(1, 1) + delta_theta * cos(theta) * aimag(jzeta - 1.0/zeta)
      integrals(2, 1) = integrals(2, 1) + delta_theta * cos(theta) * aimag(exp_zeta)
      ! #ifdef xie_correction
      ! integrals(1, 2) = integrals(1, 2) + delta_theta * cos(theta) * real(jzeta)
      ! #else
      integrals(1, 2) = integrals(1, 2) + delta_theta * cos(theta) * real(jzeta - 1.0/zeta)
      ! #endif
      integrals(2, 2) = integrals(2, 2) + delta_theta * cos(theta) * real(exp_zeta)
    enddo

  contains

    pure function exp_e1(zz)
      ! estimation of exp(z)·E1(z) where E1(z) = ∫_z^∞ exp(-t)/t dt
      ! see p.367 of G. Delhommeau thesis (referenced as [del]).
      ! the computation is done is single precision.

      complex(kind=pre), intent(in)  :: zz
      complex(kind=pre)              :: exp_e1
      complex(kind=kind(1e0))        :: g, y, z

      z = cmplx(zz, kind=kind(1e0))
      if (real(z) < -16.0) then                                      ! case 1 p. 368 in [del]
        y = 1./z
        g = y*(1.+y*(-1.+y*(2.+y*(-6.+y*(24.+y*(-120.))))))
      else if (abs(aimag(z)) > 10.0) then                            ! case 3 p. 368 in [del]
        g = 0.711093/(z+0.415775)+0.278518/(z+2.29428)+0.010389/(z+6.2900)
      else if (real(z) > -0.5) then                                  ! case 2 p. 368 in [del]
        g = -(log(z)*(.1e+01+z*(0.23721365e+00+z*(0.206543e-01+z*(0.763297e-03+     &
          z*0.97087007e-05))))+0.5772156649e+00*(0.99999207e+00+                    &
          z*(-0.149545886e+01+z*(0.41806426e-01+z*(-0.3000591e-01+                  &
          z*(0.19387339e-02+z*(-0.51801555e-03)))))))/(0.1e+01+                     &
          z*(-0.76273617e+00+z*(0.28388363e+00+z*(-0.66786033e-01+z*(0.12982719e-01 &
          +z*(-0.8700861e-03+z*0.2989204e-03))))))
      else                                                           ! case 4 p. 369 in [del]
        if (aimag(z) < 0) then
          g = ((((((( (1.000000, 1.3935496e-06)*z+ (15.82958, -20.14222))   &
            *z+ (-70.52863, -227.9511))*z+ (-985.4221, -226.6272))*z        &
            + (-1202.318, 1580.907))*z+ (953.2441, 1342.447))*z             &
            + (417.3716, -196.6665))*z+ (-9.881266, -24.24952))/            &
            (((((((( (1.000000, 0.0000000e+00)*z+ (16.83184, -20.14481))*z  &
            + (-55.66969, -248.1167))*z+ (-1068.640, -434.4707))*z          &
            + (-2082.250, 1522.471))*z+ (383.3455, 2730.378))*z             &
            + (1216.791, 351.7189))*z+ (115.3926, -161.2647))*z             &
            + (-3.777369, -4.510900))
        else
          g = ((((((( (1.000000, -1.3935496e-06)*z+ (15.82958, 20.14222))   &
            *z+ (-70.52863, 227.9511))*z+ (-985.4221, 226.6272))*z          &
            + (-1202.318, -1580.907))*z+ (953.2441, -1342.447))*z           &
            + (417.3716, 196.6665))*z+ (-9.881266, 24.24952))/              &
            (((((((( (1.000000, 0.0000000e+00)*z+ (16.83184, 20.14481))*z   &
            + (-55.66969, 248.1167))*z+ (-1068.640, 434.4707))*z            &
            + (-2082.250, -1522.471))*z+ (383.3455, -2730.378))*z           &
            + (1216.791, -351.7189))*z+ (115.3926, 161.2647))*z             &
            + (-3.777369, 4.510900))
        end if
      end if
      exp_e1 = cmplx(g, kind=pre)

    end function exp_e1

  end function numerical_integration

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function asymptotic_approximations(r, z) result(integrals)
    ! Compute the integrals using an approximate expression for large r and |z|
    real(kind=pre), intent(in) :: r
    real(kind=pre), intent(in) :: z

    real(kind=pre), dimension(2, 2) :: integrals

    real(kind=pre) :: r1, expz_sqr, sin_kr, cos_kr

    r1 = hypot(r, z)

    expz_sqr = exp(z) * sqrt(2*pi/r)
    cos_kr  = cos(r - pi/4)
    sin_kr  = sin(r - pi/4)

    integrals(1, 1) = pi*(expz_sqr*(cos_kr - sin_kr/(2*r)) - r/r1**3)
    integrals(2, 1) =     expz_sqr*(sin_kr + cos_kr/(2*r))
    ! #ifdef xie_correction
    ! integrals(1, 2) = pi*(-expz_sqr*sin_kr + z/r1**3 - one/r1)
    ! #else
    integrals(1, 2) = pi*(-expz_sqr*sin_kr + z/r1**3)
    ! #endif
    integrals(2, 2) =     expz_sqr*cos_kr

  end function asymptotic_approximations


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function construct_tabulation(r_range, z_range, nb_integration_points) result(tabulation)
    integer,        intent(in) :: nb_integration_points
    real(kind=pre), dimension(:), intent(in) :: r_range
    real(kind=pre), dimension(:), intent(in) :: z_range
    real(kind=pre), dimension(size(r_range), size(z_range), 2, 2) :: tabulation

    integer :: i, j

    do concurrent (j = 1:size(z_range))
      do concurrent (i = 1:size(r_range))
        tabulation(i, j, :, :) = numerical_integration(r_range(i), z_range(j), nb_integration_points)
      enddo
    enddo

  end function construct_tabulation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function default_z_spacing(nz)
    integer, intent(in) :: nz
    real(kind=pre), dimension(nz) :: default_z_spacing
    integer :: j
    do concurrent (j = 1:nz)
      default_z_spacing(j) = -min(10**(j/5.0-6.0), 10**(j/8.0-4.5))
      ! change of slope at z = -1e-2
    enddo
  end function default_z_spacing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function default_r_spacing(nr)
    integer, intent(in) :: nr
    real(kind=pre), dimension(nr) :: default_r_spacing
    integer :: i
    do concurrent (i = 1:nr)
      default_r_spacing(i) = min(                         &
                                 10**((i-1.0)/5.0 - 6.0), &
                                 abs(i-32)/3.0 + 4.0/3.0  &
                               )
      ! change of slope at r = 1.0
    enddo
  end function default_r_spacing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function pick_in_default_tabulation(r, z, r_range, z_range, tabulation) result(integrals)
    ! inputs
    real(kind=pre), intent(in) :: r
    real(kind=pre), intent(in) :: z
    real(kind=pre), dimension(:), intent(in) :: r_range
    real(kind=pre), dimension(:), intent(in) :: z_range
    real(kind=pre), dimension(size(r_range), size(z_range), 2, 2), intent(in) :: tabulation

    ! output
    real(kind=pre), dimension(2, 2) :: integrals

    ! local variables
    integer :: i, j

    i = max(2, min(size(r_range)-1, nearest_r_index(r)))
    j = max(2, min(size(z_range)-1, nearest_z_index(z)))

    integrals(:, :) = lagrange_polynomial_interpolation( &
      r, z,                                              &
      r_range(i-1:i+1), z_range(j-1:j+1),                &
      tabulation(i-1:i+1, j-1:j+1, :, :)                 &
      )

  contains

    pure function nearest_r_index(r)
      real(kind=pre), intent(in) :: r
      integer :: nearest_r_index

      if (r < 1e-6) then
        nearest_r_index = 2
      else if (r < 1.0) then
        nearest_r_index = int(5*(log10(r) + 6) + 1)
      else
        nearest_r_index = int(3*r + 28)
      endif
    end function

    pure function nearest_z_index(z)
      real(kind=pre), intent(in) :: z
      integer :: nearest_z_index
      real(kind=pre) :: absz

      absz = abs(z)

      if (absz > 1e-2) then
        nearest_z_index = int(8*(log10(absz) + 4.5))
      else
        nearest_z_index = int(5*(log10(absz) + 6))
      endif
    end function

    pure function lagrange_polynomial_interpolation(r, z, local_r_range, local_z_range, local_tabulation) result(integrals)
      ! inputs
      real(kind=pre),                        intent(in) :: r, z
      real(kind=pre), dimension(3),          intent(in) :: local_r_range
      real(kind=pre), dimension(3),          intent(in) :: local_z_range
      real(kind=pre), dimension(3, 3, 2, 2), intent(in) :: local_tabulation

      ! output
      real(kind=pre), dimension(2, 2) :: integrals

      ! local variable
      real(kind=pre), dimension(3) :: xl, zl

      xl(1) = pl2(local_r_range(2), local_r_range(3), local_r_range(1), r)
      xl(2) = pl2(local_r_range(3), local_r_range(1), local_r_range(2), r)
      xl(3) = pl2(local_r_range(1), local_r_range(2), local_r_range(3), r)
      zl(1) = pl2(local_z_range(2), local_z_range(3), local_z_range(1), z)
      zl(2) = pl2(local_z_range(3), local_z_range(1), local_z_range(2), z)
      zl(3) = pl2(local_z_range(1), local_z_range(2), local_z_range(3), z)

      integrals(1, 1) = dot_product(xl, matmul(local_tabulation(:, :, 1, 1), zl))
      integrals(2, 1) = dot_product(xl, matmul(local_tabulation(:, :, 2, 1), zl))
      integrals(1, 2) = dot_product(xl, matmul(local_tabulation(:, :, 1, 2), zl))
      integrals(2, 2) = dot_product(xl, matmul(local_tabulation(:, :, 2, 2), zl))
    end function lagrange_polynomial_interpolation

    pure function pl2(u1, u2, u3, xu)
      real(kind=pre), intent(in) :: u1, u2, u3, xu
      real(kind=pre) :: pl2
      pl2 = ((xu-u1)*(xu-u2))/((u3-u1)*(u3-u2))
    end function pl2

  end function pick_in_default_tabulation

end module delhommeau_integrals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
    use delhommeau_integrals

    integer :: i, j
    real(kind(1d0)) :: r, z
    real(kind(1d0)), dimension(328)           :: r_range
    real(kind(1d0)), dimension(46)            :: z_range
    real(kind(1d0)), dimension(328, 46, 2, 2) :: tabulation

    r_range = default_r_spacing(328)
    z_range = default_z_spacing(46)
    tabulation = construct_tabulation(r_range, z_range, 250)

    i = 300
    j = 40
    r = r_range(i)
    z = z_range(j)
    print*, tabulation(i, j, :, :)
    print*, numerical_integration(r, z)
    print*, asymptotic_approximations(r, z)
    print*, pick_in_default_tabulation(r, z, r_range, z_range, tabulation)

end program
