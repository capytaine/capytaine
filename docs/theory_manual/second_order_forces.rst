
===================
Second order forces
===================

Generalities
============

Equation of motion in perturbation theory
-----------------------------------------

Time domain and frequency domain
--------------------------------

Time domain interpretation of first order excitation force resulting from incoming plane wave with direction :math:`\beta` and frequency :math:`\omega`

.. math::
   F^{(1)}(t, \beta) = \Re \left( \hat F^{(1)} (\omega, \beta) e^{- i \omega t} \right)

In a mixture of plane waves with amplitudes :math:`a(\omega, \beta)`, the total excitation force is computed linearly as:

.. math::
   F^{(1)}(t) = \sum_i \Re \left( a(\omega_i, \beta_i) \hat F^{(1)}(\omega_i, \beta_i) e^{- i \omega_i t} \right)


Time domain interpretation of second order excitation force resulting from the interaction of two incoming plane waves of respectives directions :math:`\beta_1` and :math:`\beta_2` and respective frequencies :math:`\omega_1` and :math:`\omega_2`

.. math::
   F^{(2+)}(t, \beta_1, \beta_2) = \Re \left( \hat F^{(2+)} (\omega_1, \omega_2, \beta_1, \beta_2) e^{- i (\omega_1 + \omega_2) t} \right)

.. math::
   F^{(2-)}(t, \beta_1, \beta_2) = \Re \left( \hat F^{(2-)} (\omega_1, \omega_2, \beta_1, \beta_2) e^{- i |\omega_1 - \omega_2| t} \right)

In a mixture of plane waves, the total excitation force is computed bi-linearly as:

.. math::
   :nowrap:

   \begin{align*}
      F^{(2+)}(t) & = \sum_{i, j} \Re \left( a(\omega_i, \beta_i) \, a(\omega_j, \beta_j) \, \hat F^{(2+)}(\omega_i, \omega_j, \beta_i, \beta_j) \, e^{- i (\omega_i + \omega_j) t} \right) \\
      F^{(2-)}(t) & = \sum_{i, j} \Re \left( a(\omega_i, \beta_i) \, a(\omega_j, \beta_j)^* \, \hat F^{(2-)}(\omega_i, \omega_j, \beta_i, \beta_j) \, e^{- i |\omega_i - \omega_j| t} \right)
   \end{align*}

where the star is the complex conjugate.

The mean drift force is the second-order force associated with the difference frequency of a frequency with itself, leading to a time-independant term.

.. math::
   \hat{F}^{\text{mean drift}} (\omega, \beta_1, \beta_2) = \hat F^{(2-)} (\omega, \omega, \beta_1, \beta_2)


Output format
-------------

   TODO

Evaluation of second order forces
=================================

Far-field mean drift force
--------------------------

The expression of the mean drift force for the degrees of freedom Surge and Sway is the following:

.. math::
   \left\langle{\begin{matrix} F_x \\ F_y \end{matrix}}\right\rangle =
   -2 \pi \rho \omega\binom{\cos \beta}{\sin \beta} \Im( H(\beta))
   -2 \pi \rho \frac{k\left(k_0 h\right)^2}{h\left[\left(k h\right)^2-\left(k_0 h\right)^2+k_0 h\right]} \int_0^{2 \pi}|H(\theta)|^2\binom{\cos \theta}{\sin \theta} d \theta

where :math:`\beta` is the wave direction, :math:`H` the Kochin function, :math:`k` the wavenumber, :math:`h` the water depth and :math:`k_0` the deep water wavenumber.

.. note::
   The coefficient in front of the integral above can be very large at high frequency and can make the result very sensitive to small numerical inaccuracies in the Kochin function.
   Unfortunately, numerical incurracies can be common at high frequency, that is when the wavelength is small with respect to the mesh resolution.
   In other words, mesh convergence should be checked carefully for the mean drift force at high frequency.

Here is the expression for the Yaw moment:

.. math::
   \left\langle{M_z}\right\rangle = 2 \pi \frac{\rho \omega}{k}\Re (\dot H(\beta)) -
   \frac{2 \pi \rho (k_0h)^2}{h[(kh)^2 - (k_0h)^2 + k_0h]}\Im (\int_0^{2 \pi} H(\theta)^* \dot H(\theta) \mathrm{d} \theta )

The Kochin function has to be rebuild from the contributions of all the radiation problems and the diffraction problem:

.. math::
   H(\theta) = e^{i\frac{\pi}{2}} ( H_D(\theta) + \sum_{k=1}^6 X_k H_{R_k} (\theta))

where :math:`X_k` is the motion RAO of the body corresponding to the degree of freedom :math:`k`,
:math:`H_{R_k}` is the Kochin function associated with the radiated potential of degree of freedom :math:`k`
and :math:`H_{D}` is the Kochin function associated with the diffracted potential.
