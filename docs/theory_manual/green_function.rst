================================
Expression of the Green function
================================

In infinite depth
=================

The boundary integral equations relate the potential :math:`\Phi` to the normal velocity
:math:`u \cdot n` via the Green function :math:`G`.
Let us now discuss the evaluation of this function for an infinite water depth.
See also [X18]_ for a review of several expression and evaluation methods for this expression.

Green function's Hankel form
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The infinite depth Green function takes the following form

.. math::
   G(\xi, x) = - \frac{1}{4 \pi} \left( \frac{1}{\|x - \xi\|} + k \mathcal{G}\left(k \sqrt{(x_1 - \xi_1)^2 + (x_2 - \xi_2)^2}, k (x_3 + \xi_3) \right) \right)
   :label: green_function_inf_depth

The first term of :math:`G` is the usual Green function for the 3D Laplace equation without our specific boundary conditions.
The :math:`\mathcal{G}` term is complex-valued and it is introduced to satisfy the boundary conditions :eq:`bc_fs`.

.. prf:property::
   The infinite depth Green function :math:`G` is scaling-invariant in the sense of

   .. math::

        \forall x, \xi, \quad G(x, \xi, k) = G(k \xi, k x, 1).

   as well as symmetric in the sense of

   .. math::

        \forall x, \xi, \quad G(x, \xi) = G(\xi, x).

   The first term of :eq:`green_function_inf_depth` is invariant under all rotations and translations, whereas the other terms are invariant under isometric transformations that don't change the vertical coordinate (reflection across a vertical plane, rotation around a vertical axis, translation following an horizontal vector).

Introducing the dimensionless variables :math:`r = k \sqrt{(\xi_1 - x_1)^2 + (\xi_2 - x_2)^2}` and :math:`z = k (x_3 + \xi_3)`, the wave term reads in its usual `Hankel transform` form:

.. math::
    \mathcal{G}(r, z) = \frac{1}{k} \int_0^\infty \frac{\kappa + k}{\kappa - k} \exp \left(\kappa z \right) J_0 \left(\kappa r \right) d \kappa
    :label: green_function_inf_depth_hankel_form

The integrand above is singular for :math:`\kappa = k`.
Avoiding the singularity with a proper integral path avoiding the singularity in the complex plane leads to

.. math::
    \mathcal{G}(r, z) & = \frac{1}{k} \mathop{PV} \int_0^\infty \frac{\kappa + k}{\kappa - k} \exp \left(\kappa z \right) J_0 \left(\kappa r \right) d \kappa
    & \qquad \qquad \qquad \qquad + 2 \pi i e^{z} J_0(r)

where :math:`\mathop{PV}` stands for the principal value of the integral.

Singularities extraction
~~~~~~~~~~~~~~~~~~~~~~~~

.. prf:lemma::
   :label: lipschitz-integral

    Lipschitz integral:

    .. math::
       \frac{1}{\sqrt{r^2 + z^2}} = \int_0^\infty \exp \left(\kappa z \right) J_0 \left(\kappa r \right) d \kappa

Using the Lipschitz integral, the wave part of the infinite depth Green function can be rewritten as:

.. math::
    \mathcal{G}(r, z) & = \frac{1}{\sqrt{\tilde{r}^2 + \tilde{z}^2}} + \frac{2}{k} \int_0^\infty \frac{k}{\kappa - k} \exp \left(\kappa z \right) J_0 \left(\kappa r \right) d \kappa
                      & = \frac{1}{\sqrt{\tilde{r}^2 + \tilde{z}^2}} + \mathcal{G}^+(r, z)
    :label: green_function_inf_depth_hankel_form_low_freq

that is

.. math::
   -4 \pi G(x, \xi) = \frac{1}{|x - \xi|} + \frac{1}{|x - S_0(\xi)|} + k \, \mathcal{G}^+(kr, kz)

where :math:`S_0(\xi)` is the reflection of :math:`\xi` across the free surface :math:`S_0(\xi) = (\xi_1, \xi_2, -\xi_3)`.

The newly introduced term is a similar to the Rankine term :math:`\frac{1}{|x - \xi|}` except one point has been mirrored across the free surface, and is thus referred to as the `reflected Rankine term`.
It is singular when `x = S_0(\xi)`, which only appends when computing the interaction of a panel horizontal on the free surface with itself.
Note that for horizontal panels on the free surface, there is another singularity in :math:`\mathcal{G}^+(kr, kz)` which needs to be handled.

.. note::
   For convenience, the reflected Rankine term :math:`\frac{1}{|x - S_0(\xi)|}` can also be written as :math:`\frac{1}{|S_0(x) - \xi|}`.
   This is especially useful when computing the integral of this term :math:`\int_{\Gamma} \frac{1}{|x - S_0(\xi)|} d\xi`.

Alternatively, the Lipschitz integral can also be used to rewrite :eq:`green_function_inf_depth_hankel_form` as

.. math::
    \mathcal{G}(r, z) & = - \frac{1}{\sqrt{\tilde{r}^2 + \tilde{z}^2}} + \frac{2}{k} \int_0^\infty \frac{\kappa}{\kappa - k} \exp \left(\kappa z \right) J_0 \left(\kappa r \right) d \kappa
                      & = - \frac{1}{\sqrt{\tilde{r}^2 + \tilde{z}^2}} + \mathcal{G}^-(r, z)
    :label: green_function_inf_depth_hankel_form_high_freq

that is

.. math::
   -4 \pi G(x, \xi) = \frac{1}{|x - \xi|} - \frac{1}{|x - s(\xi)|} + k \, \mathcal{G}^-(kr, kz)


The notation :math:`\mathcal{G}^-` and :math:`\mathcal{G}^+` are based on the one from [X18]_.
In Capytaine, the variants above are referred to as the `low-frequency variant` for :eq:`green_function_inf_depth_hankel_form_low_freq` and the `high-frequency variant` for :eq:`green_function_inf_depth_hankel_form_high_freq`, due to the asymptotic behavior of the Green function:

.. math::
    -4 \pi G(x, \xi) \rightarrow_{k \rightarrow 0} \frac{1}{|x - \xi|} + \frac{1}{|x - s(\xi)|}

that is

.. math::
    k \mathcal{G}^+(r, z) \rightarrow_{k \rightarrow 0} 0

and

.. math::
    -4 \pi G(x, \xi) \rightarrow_{k \rightarrow \infty} \frac{1}{|x - \xi|} - \frac{1}{|x - s(\xi)|}

that is

.. math::
    k \mathcal{G}^-(r, z) \rightarrow_{k \rightarrow \infty} 0

As discussed in [A24]_, working with the low-frequency variant is usually more accurate and is thus the default in recent version of Capytaine.

.. prf:property::
   :label: green_function_inf_depth_z_derivative

   From the definitions of :math:`\mathcal G^+` and :math:`\mathcal G^-`, we have

   .. math::
     \frac{\partial \mathcal G^+}{ \partial \tilde z} = \mathcal G ^- = \mathcal{G}^+(r, z) + \frac{2}{\sqrt{ r^2 + z^2}}



Guével-Delhommeau formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. prf:lemma::
   :label: integral-form-of-Bessel

   Integral form of the Bessel function:

   .. math::
      J_0(x) = \frac{2}{\pi} \Re \int_0^{\pi/2} \exp\left(i x \cos \theta \right) d\theta


Integrating the integral form of the Bessel function into :eq:`green_function_inf_depth_hankel_form_low_freq`

.. math::
    \mathcal{G}^+(r, z) & = \frac{4}{\pi} \mathop{PV} \int_0^\infty \Re \int_0^{\pi/2} \frac{\exp \left(\kappa \left(z + i r \cos\theta \right) \right)}{\kappa - k} d\theta d\kappa
    & \qquad \qquad \qquad \qquad + 4 i \Re \int_0^{\pi/2} \exp( z + i  r \cos \theta) d\theta

which can be rewritten as

.. math::
    \mathcal{G}^+(r, z) & = \frac{4}{\pi} \Re \left( \int^0_{-\pi/2} e^\zeta(r, z, \theta) \left[ E_1(\zeta(r, z, \theta)) + i \pi \right] ) \, \mathrm{d} \theta \right) \\
    & \qquad \qquad \qquad \qquad + 4 i \Re \left( \int^{0}_{-\pi/2} e^{\zeta (r, z, \theta)} \, \mathrm{d} \theta \right)
    :label: green_function_inf_depth_guevel_low_freq

where

.. math::
    \zeta (r, z, \theta) = z + i r \cos \theta.
    :label: def_zeta

and :math:`E_1` is the first exponential integral, defined as

.. math::
   E_1(\zeta) = \int_\zeta^\infty \frac{e^{-t}}{t} \mathrm{d} t.

where we used the following property

.. math::
   \int_0^\infty  \frac{e^{-k z}}{k - a} dk =
   \begin{cases}
     e^{-a z} (E_1(az) + i \pi) & \Im(z) \ge 0 \\
     e^{-a z} (E_1(az) - i \pi) & \Im(z) < 0
   \end{cases}

.. note::
   In [Del87]_ integrals of the form :math:`\int_{-\pi/2}^{\pi/2} \ldots d \theta` are used.
   Given the parity of the integrand with respect to :math:`\theta`, we prefer to simplify them as :math:` 2\int_{0}^{\pi/2} \ldots d \theta`.

.. prf:lemma::
   :label: integrate-one-over-zeta

    Using the integral form of the Bessel function :prf:ref:`integral-form-of-Bessel` into :prf:ref:`lipschitz-integral` gives

    .. math::
       \frac{2}{\pi} \Re \int_{0}^{\pi/2} \frac{1}{\zeta(\theta)} \, \mathrm{d} \theta = - \frac{1}{\sqrt{r^2 + z^2}}.
       :label: int_1_over_zeta

    which can also be found in [Del89]_

The above lemma allows to retrieve the expression of the Green function found e.g. in [BD15]_:

.. math::
    \mathcal{G}^-(r, z) & =  \frac{4}{\pi} \Re \left( \int^{0}_{-\pi/2} \left( e^\zeta(r, z, \theta) \left[ E_1(\zeta(r, z, \theta)) + i \pi \right] - \frac{1}{\zeta(r, z, \theta)} \right) \, \mathrm{d} \theta \right) \\
    & \qquad \qquad \qquad \qquad + 2 i \Re \left( \int^{\pi/2}_{-\pi/2} e^{\zeta (r, z, \theta)} \, \mathrm{d} \theta \right)
    :label: green_function_inf_depth_high_freq

Alternative formulations
~~~~~~~~~~~~~~~~~~~~~~~~

.. prf:lemma::

    The `zeroth order Bessel function of the first kind <https://personal.math.ubc.ca/~cbm/aands/page_360.htm>`_ :math:`J_0` and `the Struve function <https://personal.math.ubc.ca/~cbm/aands/page_496.htm>`_ :math:`H_0` are such that

    .. math::
        J_0(r) & = \frac{2}{\pi} \int_{0}^{\pi/2} \cos(r\cos(\theta)) \, \mathrm{d} \theta \\
        H_0(r) & = \frac{2}{\pi} \int_{0}^{\pi/2} \sin(r\cos(\theta)) \, \mathrm{d} \theta \\

    hence

    .. math::
        2 \int_{0}^{\pi/2} i e^{\zeta} \, \mathrm{d} \theta = \pi e^z \left(- H_0(r) + i J_0(r) \right)


The function :math:`\mathcal{G}` can also be rewritten as

.. math::
    \mathcal{G}(r, z) & = \frac{1}{\sqrt{r^2 + z^2}} + \frac{2}{\pi} \int^{\pi/2}_{-\pi/2} \Re \left( e^\zeta E_1(\zeta) \right) \, \mathrm{d} \theta + 2 \int^{\pi/2}_{-\pi/2} i e^{\zeta (r, z, \theta)} \, \mathrm{d} \theta \\
    & = \frac{1}{\sqrt{r^2 + z^2}} + \frac{2}{\pi} \int^{\pi/2}_{-\pi/2} \Re \left( e^\zeta E_1(\zeta) \right) \, \mathrm{d} \theta + 2 \pi e^z \left( - H_0(r) + i J_0(r) \right)

Noblesse [N82]_ splits the function :math:`\mathcal{G}` into a near field term :math:`N` and a wave field :math:`W` such that

.. math::
   N(r, z) & = \frac{1}{\sqrt{r^2 + z^2}} + \frac{2}{\pi} \int^{\pi/2}_{-\pi/2} \Re \left( e^\zeta E_1(\zeta)  \right) \, \mathrm{d} \theta  \\
   W(r, z) & = 2 \pi e^z \left( - H_0(r) + i J_0(r) \right)

Note that :math:`E_1`, :math:`J_0` and :math:`H_0` are available for instance in the `Scipy library <https://docs.scipy.org/doc/scipy/reference/special.html>`_.

The literature also contains alternative formulation with another variant of the :math:`\zeta` function as seen in the lemma below.

.. prf:lemma::

    For any function :math:`f`, the following two formulations of the integral are equivalent [Del89]_:

    .. math::
        \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} f \left(\zeta(\theta) \right) \mathrm{d} \theta =
        \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} f \left(\tilde{\zeta}(\theta) \right) \mathrm{d} \theta

    where :math:`\zeta` is defined in :eq:`def_zeta` and :math:`\tilde{\zeta}` is defined as

    .. math::
       \tilde{\zeta} (\theta) = k \left( x_3 + \xi_3 + i \left( (x_1 - \xi_1) \cos\theta + (x_2 - \xi_2) \sin\theta \right) \right).


.. prf:proof::
   Let us note that:

   .. math::
      :nowrap:

      \begin{align*}
      (x_1 - \xi_1) \cos(\theta) + (x_2 - \xi_2) \sin(\theta) & = \Re \left( \left( x_1 - \xi_1  + i (x_2 - \xi_2) \right) e^{-i \theta} \right) \\
                   & = \Re \left( r e^{i (\alpha - \theta)} \right) \\
                   & = r \cos \left( \alpha - \theta \right) \\
      \end{align*}

   where :math:`r` and :math:`\alpha` are defined by

   .. math::
      :nowrap:

      \[
          r e^{i \alpha} = (x_1 - \xi_1)  + i (x_2 - \xi_2).
      \]

   Finally note that:

    .. math::
        :nowrap:

        \[
            \int_{-\frac{\pi}{2}-\alpha}^{\frac{\pi}{2}-\alpha} f \left(\zeta(\theta) \right) \mathrm{d} \theta =
            \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} f \left(\zeta(\theta) \right) \mathrm{d} \theta
        \]


Gradient of the Green function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The gradient of the Green function with respect to its first variable (that is :math:`x`) can be written as

.. math::
   \nabla_1 G(x, \xi) = - \frac{1}{4 \pi} \left( - \frac{x - \xi}{\|x - \xi\|^3} + k
      \begin{pmatrix}
        \frac{\partial r}{\partial x_1} \frac{\partial \mathcal{G}}{\partial r} \\
        \frac{\partial r}{\partial x_2} \frac{\partial \mathcal{G}}{\partial r} \\
        \frac{\partial z}{\partial x_3} \frac{\partial \mathcal{G}}{\partial z}
      \end{pmatrix}
   \right)

with

.. math::
   \frac{\partial r}{\partial x_1} & = k^2 \frac{x_1 - \xi_1}{r} \\
   \frac{\partial r}{\partial x_2} & = k^2 \frac{x_2 - \xi_2}{r} \\
   \frac{\partial z}{\partial x_3} & = k.

or equivalently using :math:`\mathcal G^+` from :eq:`green_function_inf_depth_low_freq`:

.. math::
   \nabla_1 G(x, \xi) = - \frac{1}{4 \pi} \left( - \frac{x - \xi}{\|x - \xi\|^3} + k
   - \frac{1}{4 \pi} \left( - \frac{x - S_0(\xi)}{\|x - S_0(\xi)\|^3} + k
      \begin{pmatrix}
        \frac{\partial r}{\partial x_1} \frac{\partial \mathcal{G}}{\partial r} \\
        \frac{\partial r}{\partial x_2} \frac{\partial \mathcal{G}}{\partial r} \\
        \frac{\partial z}{\partial x_3} \frac{\partial \mathcal{G}}{\partial z}
      \end{pmatrix}
   \right)

The derivative of :math:`\mathcal G^+` with respect to :math:`z` can be handled with :prf:ref:`green_function_inf_depth_z_derivative`.
Let us now focus on the derivative with respect to :math:`r`:

From the definition of the exponential integral, we have

.. math::
    \frac{d}{d\zeta}\left(e^\zeta \left( E_1(\zeta) + i \pi \right)\right) = e^\zeta \left( E_1(\zeta) + i \pi \right) - 1/\zeta

hence

.. math::
   \frac{\partial \mathcal{G}^+}{\partial r} = & \frac{4}{\pi} \Re \left( \int_{0}^{\pi/2} \frac{\partial \zeta}{\partial r} \left( e^\zeta \left( E_1(\zeta) + i \pi \right) - \frac{1}{\zeta} \right) \, \mathrm{d}\theta \right) \\
   & \qquad \qquad \qquad \qquad + 4 i \Re \left( \int^{\pi/2}_{0} i \cos(\theta) e^{\zeta} \, \mathrm{d} \theta \right)
   :label: green_function_inf_depth_r_derivative

.. note::
    The derivative of :math:`G` with respect to :math:`x_1` and :math:`x_2` are antisymmetric in the sense of

    .. math::
       :nowrap:

        \[
        \frac{\partial G}{\partial x_1} (\xi, x) = - \frac{\partial G}{\partial x_1}(x, \xi).
        \]

    Its derivative with respect to :math:`x_3` is symmetric in infinite depth.

    In finite depth, some terms of the derivative with respect to :math:`x_3` are symmetric and some are antisymmetric.


Higher order derivatives
~~~~~~~~~~~~~~~~~~~~~~~~

From :eq:`green_function_inf_depth_z_derivative`, one has

.. math::
   \frac{\partial \mathcal{G}}{\partial z} &= \mathcal{G}(r, z) + \left( 1 + \frac{\partial}{\partial z} \right) \frac{1}{\sqrt{r^2 + z^2}} \\
   \frac{\partial^2 \mathcal{G}}{\partial z \partial r} &= \frac{\partial \mathcal{G}}{\partial r} + \left( \frac{\partial}{\partial r} + \frac{\partial^2}{\partial z \partial r} \right) \frac{1}{\sqrt{r^2 + z^2}}

and

.. math::
   \frac{\partial^2 \mathcal{G}}{\partial z^2} &= \mathcal{G}(r, z) + \left( 1 + 2 \frac{\partial}{\partial z} + \frac{\partial^2}{\partial z^2} \right)\frac{1}{\sqrt{r^2 + z^2}} \\
                                               &= \mathcal{G}(r, z) + \frac{1}{\sqrt{r^2 + z^2}} - 2 \frac{z}{(r^2 + z^2)^{3/2}} - \frac{r^2 - 2 z^2}{(r^2 + z^2)^{5/2}}

Since the Green function is solution of the Laplace equation, it follows that

.. math::
   \frac{\partial^2 \mathcal{G}}{\partial r^2} + \frac{1}{r} \frac{\partial \mathcal{G}}{\partial r} + \frac{\partial^2 \mathcal{G}}{\partial z^2} = 0

then

.. math::
   \frac{\partial^2 \mathcal{G}}{\partial r^2} = - \frac{1}{r} \frac{\partial \mathcal{G}}{\partial r} - \mathcal{G} - \left( 1 + 2 \frac{\partial}{\partial z} + \frac{\partial^2}{\partial z^2} \right)\frac{1}{\sqrt{r^2 + z^2}} \\

All higher order derivative can be expressed with the help of :math:`\mathcal{G}` and :math:`\frac{\partial \mathcal{G}}{\partial r}`.

.. note::
   The same derivation is done in e.g. [N20]_ using instead the function :math:`F = \mathcal{G} - \frac{1}{\sqrt{r^2 + z^2}}` for which the expressions are slightly simpler.


Delhommeau's method for computation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current version of Capytaine can use either the low-frequency variant :eq:`green_function_inf_depth_low_freq` or high-frequency variant :eq:`green_function_inf_depth_high_freq` when evaluating the Green function and its integral over a panel.
For this purpose, the following values needs to be computed:

.. math::
    I_1(r, z) & = \frac{4}{\pi} \Re \left( \int^{\pi/2}_{0} J(\zeta) \, \mathrm{d} \theta \right) \\
    I_2(r, z) & = \frac{4}{\pi} \Re \left( \int^{\pi/2}_{0} \left( J(\zeta) - \frac{1}{\zeta} \right) \, \mathrm{d} \theta \right) \\
    I_3(r, z) & = 4 \Re \left( \int^{\pi/2}_{0} e^{\zeta} \, \mathrm{d} \theta \right) \\
    I_4(r, z) & = \frac{4}{\pi} \Re \left( \int^{\pi/2}_{0} i \cos(\theta) \left( J(\zeta) - \frac{1}{\zeta} \right) \, \mathrm{d} \theta \right) \\
    I_5(r, z) & = 4 \Re \left( \int^{\pi/2}_{0} i \cos(\theta) e^{\zeta} \, \mathrm{d} \theta \right)

then :eq:`green_function_inf_depth_low_freq` and :eq:`green_function_inf_depth_r_derivative` reads

.. math::
   \mathcal{G}^+(r, z) & = I_1(r, z) + i I_3(r, z) \\
   \frac{\partial \mathcal{G}^+}{\partial r} & = I_4(r, z) + i I_5(r, z).

and :eq:`green_function_inf_depth_high_freq` reads (still using :eq:`green_function_inf_depth_r_derivative` for the derivative):

.. math::
   \mathcal{G}^-(r, z) = I_2(r, z) + i I_3(r, z).
   \frac{\partial \mathcal{G}^+}{\partial r} & = I_4(r, z) + i I_5(r, z).

To limit the computational cost of the evaluation of these integrals, they are precomputed for selected values of :math:`r` and :math:`z` and stored in a table.
When evaluating the Green function, the values of the integrals are retrieved by interpolating the values in the tables.

For large values of :math:`r` and :math:`z`, these integrals are asymptotically approximated by the following expressions:

.. math::
      I_1(r, z) & \simeq - 2 e^z \sqrt{\frac{2\pi}{r}} \sin(r - \pi/4) + \frac{2 z}{(r^2 + z^2)^{3/2}} - \frac{2}{\sqrt{r^2 + z^2}} \\
      I_2(r, z) & \simeq - 2 e^z \sqrt{\frac{2\pi}{r}} \sin(r - \pi/4) + \frac{2 z}{(r^2 + z^2)^{3/2}} \\
      I_3(r, z) & \simeq 2 e^z \sqrt{\frac{2\pi}{r}} \cos(r - \pi/4) \\
      I_4(r, z) & \simeq - 2 e^z \sqrt{\frac{2\pi}{r}} \left( \cos(r - \pi/4) - \frac{1}{2r} \sin(r-\pi/4) \right) + \frac{2 r}{(r^2 + z^2)^{3/2}} \\
      I_5(r, z) & \simeq - 2 e^z \sqrt{\frac{2\pi}{r}} \left( \sin(r - \pi/4) + \frac{1}{2r} \cos(r - \pi/4) \right)


Incorporating these asymptotic approximation in the expression of the Green function, one gets:

.. math::
    \mathcal{G}(r, z) \simeq & -\frac{1}{\sqrt{r^2 + z^2}} - 2 e^z \sqrt{\frac{2\pi}{r}} \left(\sin(r - \pi/4) - i\cos(r - \pi/4)\right) \\
   & \qquad\qquad\qquad\qquad + 2 \frac{z}{(r^2 + z^2)^{3/2}}
   :label: green_function_asymptotical_approx


Integration
~~~~~~~~~~~

TODO

As seen in :eq:`green_function_inf_depth_z_derivative`, new reflected-Rankine-type
terms might appear in the derivative of the Green wave term.
By default, they are integrated with the same method used for the same
numerical quadrature method as the rest of the wave term.
The setting ``gf_singularities="low_freq_with_rankine_term"`` is an attempt to
integrate them exactly using the same code as the main reflected Rankine term.

In finite depth
===============

Green function
~~~~~~~~~~~~~~

TODO

Gradient of the Green function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO
