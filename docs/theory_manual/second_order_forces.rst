
===================
Second order forces
===================

Generalities
============

Equation of motion in perturbation theory
-----------------------------------------

The full equations of motion of a body interacting with waves are linearized using perturbation theory.
Denoting by :math:`M` the inertia matrix, :math:`A` the added mass matrix,
:math:`B` the radiation damping matrix and :math:`S` the hydrostatic stiffness
matrix, we have at first order in frequency domain:

.. math::
   \left(- \omega^2 (M + A(\omega)) - i \omega B(\omega) + S\right) \hat q^{(1)}(\omega) = \hat F^{(1)}_\text{excitation}(\omega)

(or the Cummins equation equivalent in time domain) and at second order

.. math::
   \left(- \omega^2 (M + A(\omega)) - i \omega B(\omega) + S\right) \hat q^{(2)}(\omega) = \hat F^{(2)}_\text{excitation}(\hat q^{(1)}, \omega) + ...

where :math:`...` encompasses all possible other second-order terms, such as gyroscopic effect or non-linear internal stiffness for elastic bodies.

Both equations hold for all frequencies :math:`\omega`, with the same operator on the left-hand side.
They differ by the frequency content of their right-hand side.
For a single incoming wave of frequency :math:`\omega_i`, the first order force :math:`\hat F^{(1)}` is non-zero at :math:`\omega_i`, while the second-order force :math:`\hat F^{(2)}` is non-zero at frequencies :math:`0` and :math:`2\omega_i`.

This document is about the computation of :math:`F^{(2)}_\text{excitation}(q^{(1)})`.

Time domain and frequency domain
--------------------------------

Time domain interpretation of first order excitation force resulting from incoming plane wave with direction :math:`\beta` and frequency :math:`\omega`

.. math::
   F^{(1)}(t, \beta) = \Re \left( \hat F^{(1)} (\omega, \beta) e^{- i \omega t} \right)

In a mixture of plane waves with amplitudes :math:`a(\omega, \beta)`, the total excitation force is computed linearly as:

.. math::
   F^{(1)}(t) = \sum_i \Re \left( a(\omega_i, \beta_i) \hat F^{(1)}(\omega_i, \beta_i) e^{- i \omega_i t} \right)


In a mixture of plane waves, the total second order excitation force is computed bi-linearly as:

.. math::
   :nowrap:

   \begin{align}
      F^{(2)}(t) & = \sum_{i, j} \Re \left( a(\omega_i, \beta_i) \, a(\omega_j, \beta_j) \, \hat F^{(2+)}(\omega_i, \omega_j, \beta_i, \beta_j) \, e^{- i (\omega_i + \omega_j) t} \right) \\
                 & + \sum_{i, j} a(\omega_i, \beta_i) \, a(\omega_j, \beta_j)^* \, \hat F^{(2-)}(\omega_i, \omega_j, \beta_i, \beta_j) \, e^{- i (\omega_i - \omega_j) t}
   \end{align}

where the star is the complex conjugate and the sums are over all the pairs :math:`(i, j)`, including :math:`i = j`.
The transfer functions :math:`\hat F^{(2+)}` and :math:`\hat F^{(2-)}` are called QTFs (Quadratic Transfer Functions).
The sum frequency transfer function is symmetric and the difference frequency transfer function is Hermitian:

.. math::
   \hat F^{(2+)}(\omega_j, \omega_i, \beta_j, \beta_i) = \hat F^{(2+)}(\omega_i, \omega_j, \beta_i, \beta_j)

.. math::
   \hat F^{(2-)}(\omega_j, \omega_i, \beta_j, \beta_i) = \hat F^{(2-)}(\omega_i, \omega_j, \beta_i, \beta_j)^*.

Hence the terms :math:`(i, j)` and :math:`(j, i)` of the difference frequency sum are complex conjugate, and their sum is real without taking the real part.
The exponent :math:`\omega_i - \omega_j` of a single term can be negative.
For :math:`\omega_i \neq \omega_j`, only the sum of the two :math:`i,j` and :math:`j,i` terms has a physical meaning: it is a term oscillating at the positive difference frequency :math:`|\omega_i - \omega_j|`.

The mean drift force is the second-order force associated with the difference frequency of a frequency with itself, leading to a time-independent term, that is the mean value of the second-order force:

.. math::
   \langle F^{(2)} \rangle = \sum_k \sum_{i, j} a(\omega_k, \beta_i) \, a(\omega_k, \beta_j)^* \, \hat F^{\text{mean drift}}(\omega_k, \beta_i, \beta_j),

with

.. math::
   \hat{F}^{\text{mean drift}} (\omega_k, \beta_i, \beta_j) = \hat F^{(2-)} (\omega_k, \omega_k, \beta_i, \beta_j)


Cross-direction components of the difference-frequency terms
------------------------------------------------------------

The quadratic second-order forces are built from products :math:`u v` of two first-order quantities
(pressure and displacement, two components of the velocity, ...).
Let :math:`\hat u_i` and :math:`\hat v_i` be the complex amplitudes of :math:`u` and :math:`v` in a single wave of direction :math:`\beta_i` and unit amplitude.
In a mixture of waves of the same frequency :math:`\omega_k` and amplitudes :math:`a_i = a(\omega_k, \beta_i)`,

.. math::
   u(t) = \Re \left( \sum_i a_i \, \hat u_i \, e^{- i \omega_k t} \right), \qquad
   v(t) = \Re \left( \sum_j a_j \, \hat v_j \, e^{- i \omega_k t} \right),

and the mean value of their product is

.. math::
   \langle u v \rangle = \frac{1}{2} \Re \left( \sum_{i, j} a_i \, a_j^* \, \hat u_i \, \hat v_j^* \right)
                       = \sum_{i, j} a_i \, a_j^* \, Q_{ij}(\hat u, \hat v)

with

.. math::
   :label: diff_freq_product_rule

   Q_{ij}(\hat u, \hat v) = \frac{1}{4} \left( \hat u_i \, \hat v_j^* + \hat u_j^* \, \hat v_i \right).

The second equality is obtained by writing :math:`\Re(z) = (z + z^*)/2` and swapping :math:`i` and :math:`j` in the term with :math:`z^*`.
This choice of :math:`Q_{ij}` among all the matrices giving the same mean value is the one that is Hermitian, :math:`Q_{ji} = Q_{ij}^*`,
such that the sum above is real term by term for each pair :math:`\{i, j\}`.
Its diagonal is the mean value in a single wave of unit amplitude:

.. math::
   Q_{ii}(\hat u, \hat v) = \frac{1}{2} \Re \left( \hat u_i \, \hat v_i^* \right).

This product rule is applied to each quadratic term of the second-order force,
and the linear operations (sums, integration over the hull or along the waterline, multiplication by real coefficients) are applied to :math:`Q_{ij}` directly.
The resulting :math:`\hat F^{\text{mean drift}}(\omega, \beta_i, \beta_j)` is then also Hermitian in :math:`(i, j)`.

The product rule :eq:`diff_freq_product_rule` derived above for the mean-drift force is also valid for the other quadratic parts of the difference-frequency forces.

.. note::
   Complex factors, such as :math:`i \omega \rho` to get the pressure from the potential,
   have to be applied to the first-order amplitudes :math:`\hat u_i` before the product.
   Multiplying :math:`Q_{ij}` by a complex number does not correspond to any physical quantity, and would break the Hermitian symmetry.


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
   Unfortunately, numerical inacuracies can be common at high frequency, that is when the wavelength is small with respect to the mesh resolution.
   In other words, mesh convergence should be checked carefully for the mean drift force at high frequency.

Here is the expression for the Yaw moment:

.. math::
   \left\langle{M_z}\right\rangle = 2 \pi \frac{\rho \omega}{k}\Re (\dot H(\beta)) -
   \frac{2 \pi \rho (k_0h)^2}{h[(kh)^2 - (k_0h)^2 + k_0h]}\Im (\int_0^{2 \pi} H(\theta)^* \dot H(\theta) \mathrm{d} \theta )

The Kochin function has to be rebuilt from the contributions of all the radiation problems and the diffraction problem:

.. math::
   H(\theta) = e^{i\frac{\pi}{2}} ( H_D(\theta) + \sum_{k=1}^6 X_k H_{R_k} (\theta))

where :math:`X_k` is the motion RAO of the body corresponding to the degree of freedom :math:`k`,
:math:`H_{R_k}` is the Kochin function associated with the radiated potential of degree of freedom :math:`k`
and :math:`H_{D}` is the Kochin function associated with the diffracted potential.
