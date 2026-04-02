! Copyright (C) 2022-2024 Matthieu Ancellin
! See LICENSE file at <https://github.com/capytaine/libDelhommeau>
!
! This module contains functions to evaluate the following integrals
! for a given range of values of `r ∈ [0, +∞)` and `z ∈ (-∞, 0]`.
! I(1) = Re[ 𝒢^+ ] = 2/π Re[ ∫ e^ζ (E(ζ) + iπ) dθ ]
! I(2) = Re[ 𝒢^- ] = 2/π Re[ ∫ (e^ζ (E(ζ) + iπ) - 1/ζ) dθ ]
! I(3) = Im[ 𝒢^+ ] = Im[ 𝒢^- ] =  2/π Re[ ∫(e^ζ) dθ ]
! I(4) = Re[ ∂𝒢^+/∂r ] = 2 Re[ ∫ (i cosθ) (e^ζ (E(ζ) + iπ) - 1/ζ) dθ ]
! I(5) = Im[ ∂𝒢^+/∂r ] = 2 Re[ ∫ (i cosθ) (e^ζ) dθ ]
! where ζ = z + i r cos θ.
!
! They are required for the evaluation of the Green function and its gradient.
!
module delhommeau_integrals

  use floating_point_precision, only: pre
  use constants

  implicit none

  public :: numerical_integration
  public :: asymptotic_approximations
  public :: construct_tabulation
  public :: default_r_spacing, default_z_spacing
  public :: pick_in_default_tabulation

  !private  ! Other functions are private by default

contains

  pure function numerical_integration(r, z, nb_integration_points) result(integrals)
    ! Compute the integrals by numerical integration, with `nb_integration_points` points.

    ! input
    real(kind=pre), intent(in) :: r
    real(kind=pre), intent(in) :: z
    integer,        intent(in) :: nb_integration_points

    ! output
    real(kind=pre), dimension(nb_tabulated_values) :: integrals

    ! local variables
    integer :: k, n
    real(kind=pre) :: theta, delta_theta, cos_theta
    complex(kind=pre) :: zeta, exp_zeta, jzeta

    ! initial values
    integrals(:) = 0.0

    ! Should have an odd number of points for Simpson rule
    if (mod(nb_integration_points, 2) == 0) then
      n = nb_integration_points + 1
    else
      n = nb_integration_points
    endif

    do k = 1, n
      theta = -pi/2 + pi*(k-1.0)/(n-1.0)
      ! The 1.0 are here on purpose to force the recasting of k and n as floats.

      ! Step of the integration using Simpson rule
      if ((k == 1) .or. (k == n)) then
        delta_theta = pi/(3*(n-1.0))
      elseif (mod(k, 2) == 0) then
        delta_theta = 4*pi/(3*(n-1.0))
      else
        delta_theta = 2*pi/(3*(n-1.0))
      endif

      cos_theta = cos(theta)
      zeta = cmplx(z, r*cos_theta, kind=pre)
      IF (REAL(ZETA) <= -30.0) THEN
        exp_zeta = (0.0, 0.0)
      else
        exp_zeta = exp(zeta)
      endif
      jzeta = exp_e1(zeta) + ii*pi*exp_zeta
      integrals(1) = integrals(1) + delta_theta * real(jzeta)
      integrals(2) = integrals(2) + delta_theta * real(jzeta - 1.0/zeta)
      integrals(3) = integrals(3) + delta_theta * real(exp_zeta)
      integrals(4) = integrals(4) - delta_theta * cos_theta * aimag(jzeta - 1.0/zeta)
      integrals(5) = integrals(5) - delta_theta * cos_theta * aimag(exp_zeta)
    enddo

    integrals(1) = integrals(1)/PI
    integrals(2) = integrals(2)/PI
    integrals(4) = integrals(4)/PI
    integrals(:) = 2*integrals(:)

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
    ! Evaluate the wave part of legacy's Delhommeau Green function
    ! using an approximate expression for large r and |z|
    real(kind=pre), intent(in) :: r
    real(kind=pre), intent(in) :: z

    real(kind=pre), dimension(nb_tabulated_values) :: integrals

    real(kind=pre) :: r1, expz_sqr, sin_kr, cos_kr

    r1 = hypot(r, z)

    expz_sqr = exp(z) * sqrt(2*pi/r)
    cos_kr  = cos(r - pi/4)
    sin_kr  = sin(r - pi/4)

    integrals(1) = -expz_sqr*sin_kr + z/r1**3 - 1/r1
    integrals(2) = -expz_sqr*sin_kr + z/r1**3
    integrals(3) =  expz_sqr*cos_kr
    integrals(4) = -expz_sqr*(cos_kr - sin_kr/(2*r)) + r/r1**3
    integrals(5) = -expz_sqr*(sin_kr + cos_kr/(2*r))
    integrals(:) = 2*integrals(:)

  end function asymptotic_approximations

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function construct_tabulation(r_range, z_range, nb_integration_points) result(tabulation)
    real(kind=pre), dimension(:), intent(in) :: r_range
    real(kind=pre), dimension(:), intent(in) :: z_range
    integer,        intent(in) :: nb_integration_points
    real(kind=pre), dimension(size(r_range), size(z_range), nb_tabulated_values) :: tabulation

    integer :: i, j

    do concurrent (j = 1:size(z_range))
      do concurrent (i = 1:size(r_range))
        tabulation(i, j, :) = numerical_integration(r_range(i), z_range(j), nb_integration_points)
      enddo
    enddo

  end function construct_tabulation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function default_r_spacing(nr, rmax, method)
    integer, intent(in) :: nr
    real(kind=pre), intent(in) :: rmax
    integer, intent(in) :: method
    real(kind=pre), dimension(nr) :: default_r_spacing

    ! Reference parameters from Nemoh 3 model
    integer, parameter :: nr_ref = 676
    integer, parameter :: index_of_1_ref = 81  ! index of the change of slope

    ! local variables
    integer :: i, index_of_1

    default_r_spacing(1) = 0.0

    if (method == LEGACY_GRID) then
      do concurrent (i = 2:nr)
        default_r_spacing(i) = min(                    &
                                 10**((i-1.0)/5.0 - 6.0), &
                                 abs(i-32)/3.0 + 4.0/3.0  &
                               )
      enddo
    else
      ! change of slope at r = 1.0 that is i=index_of_1
      index_of_1 = nint(nr*1.0/nr_ref*index_of_1_ref)
      do concurrent (i = 1:nr)
        if (i < index_of_1) then
          ! Exponential spacing
          default_r_spacing(i) = 10.0**(-10*(1-real(i-1)/(index_of_1-1)))
        else
          ! Linear spacing
          default_r_spacing(i) = (rmax-1.0)/(nr-index_of_1)*(i-index_of_1)+1.0
        endif
      enddo
    endif
  end function default_r_spacing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function default_z_spacing(nz, zmin, method)
    integer, intent(in) :: nz
    real(kind=pre), intent(in) :: zmin
    integer, intent(in) :: method
    real(kind=pre), dimension(nz) :: default_z_spacing

    integer :: j
    real(kind=pre) :: dz

    if (method == LEGACY_GRID) then
      do concurrent (j = 1:nz)
        default_z_spacing(j) = -min(10**(j/5.0-6.0), 10**(j/8.0-4.5))
        ! change of slope at z = -1e-2
      enddo
      if (nz == 46) then  ! For consistency with Nemoh 2...
        default_z_spacing(46) = -16.0
      endif
    else
      dz = (log10(abs(zmin))+10.0)/nz
      do concurrent (j = 1:nz)
        default_z_spacing(j) = -min(10**(dz*j-10.0), abs(zmin))
      enddo
    endif
  end function default_z_spacing

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function pick_in_default_tabulation(r, z, method, r_range, z_range, tabulation) result(interpolated_values)
    ! inputs
    real(kind=pre), intent(in) :: r
    real(kind=pre), intent(in) :: z
    integer, intent(in) :: method
    real(kind=pre), dimension(:), intent(in) :: r_range
    real(kind=pre), dimension(:), intent(in) :: z_range
    real(kind=pre), dimension(size(r_range), size(z_range), nb_tabulated_values), intent(in) :: tabulation

    ! output
    real(kind=pre), dimension(nb_tabulated_values) :: interpolated_values

    ! local variables
    integer :: i, j

    i = max(2, min(size(r_range)-1, nearest_r_index(r, r_range, method)))
    j = max(2, min(size(z_range)-1, nearest_z_index(z, z_range, method)))

    interpolated_values(:) = lagrange_polynomial_interpolation( &
      r, z,                                                     &
      r_range(i-1:i+1), z_range(j-1:j+1),                       &
      tabulation(i-1:i+1, j-1:j+1, :)                           &
      )
  end function pick_in_default_tabulation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function nearest_r_index(r, r_range, method)
    real(kind=pre), intent(in) :: r
    integer, intent(in) :: method
    real(kind=pre), dimension(:), intent(in) :: r_range
    integer :: nearest_r_index

    ! Reference parameters from Nemoh 3 model
    integer, parameter :: nr_ref = 676
    integer, parameter :: index_of_1_ref = 81  ! index of the change of slope

    ! local variables
    integer :: index_of_1
    real(kind=pre) :: rmax

    if (method == LEGACY_GRID) then
      if (r < 1e-6) then
        nearest_r_index = 1
      else if (r < 1.0) then
        nearest_r_index = int(5*(log10(r) + 6) + 1)
      else
        nearest_r_index = int(3*r + 28)
      endif
    else
      index_of_1 = nint(real(size(r_range)*index_of_1_ref)/nr_ref)
      rmax = r_range(size(r_range))

      if (r < 1e-10) then
        nearest_r_index = 1
      else if (r < 1.0) then
        nearest_r_index = nint((log10(r)/10.0 + 1.0)*(index_of_1-1) + 1)
      else
        nearest_r_index = nint((r - 1)*(size(r_range) - index_of_1)/(rmax - 1) + index_of_1)
      endif
    endif
  end function

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function nearest_z_index(z, z_range, method)
    real(kind=pre), intent(in) :: z
    integer, intent(in) :: method
    real(kind=pre), dimension(:), intent(in) :: z_range
    integer :: nearest_z_index

    ! local parameters
    real(kind=pre) :: absz
    real(kind=pre) :: dz
    integer :: nz

    absz = abs(z)
    nz = size(z_range)

    if (method == LEGACY_GRID) then
      if (absz > 1e-2) then
        nearest_z_index = int(8*(log10(absz) + 4.5))
      else
        nearest_z_index = int(5*(log10(absz) + 6))
      endif
    else
      dz = (log10(abs(z_range(nz)))+10.0)/nz
      nearest_z_index = nint((log10(absz)+10)/dz)
    endif
  end function

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function lagrange_polynomial_interpolation(         &
      r, z, local_r_range, local_z_range, local_tabulation &
      ) result(interpolated_values)
    ! inputs
    real(kind=pre),                        intent(in) :: r, z
    real(kind=pre), dimension(3),          intent(in) :: local_r_range
    real(kind=pre), dimension(3),          intent(in) :: local_z_range
    real(kind=pre), dimension(3, 3, nb_tabulated_values), intent(in) :: local_tabulation

    ! output
    real(kind=pre), dimension(nb_tabulated_values) :: interpolated_values

    ! local variable
    integer :: k
    real(kind=pre), dimension(3) :: xl, zl

    xl(1) = pl2(local_r_range(2), local_r_range(3), local_r_range(1), r)
    xl(2) = pl2(local_r_range(3), local_r_range(1), local_r_range(2), r)
    xl(3) = pl2(local_r_range(1), local_r_range(2), local_r_range(3), r)
    zl(1) = pl2(local_z_range(2), local_z_range(3), local_z_range(1), z)
    zl(2) = pl2(local_z_range(3), local_z_range(1), local_z_range(2), z)
    zl(3) = pl2(local_z_range(1), local_z_range(2), local_z_range(3), z)

    do concurrent (k=1:nb_tabulated_values)
      interpolated_values(k) = dot_product(xl, matmul(local_tabulation(:, :, k), zl))
    enddo

    contains

    pure function pl2(u1, u2, u3, xu)
      real(kind=pre), intent(in) :: u1, u2, u3, xu
      real(kind=pre) :: pl2
      pl2 = ((xu-u1)*(xu-u2))/((u3-u1)*(u3-u2))
    end function pl2

  end function lagrange_polynomial_interpolation

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module delhommeau_integrals
