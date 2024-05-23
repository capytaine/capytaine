=============
Theory manual
=============

This document presents several aspects of the theoretical background for Nemoh and Capytaine.
We refer also to [Del87]_ [Del89]_ [Del93]_ [BD15]_ and [AD18]_.

.. warning:: This document is a work in progress. It is incomplete and might
   contain errors.

.. contents:: Contents

Linear boundary value problem
=============================

Hypotheses
----------

1. The fluid is inviscid.
2. The fluid is incompressible: :math:`\nabla \cdot u = 0` with :math:`u` the flow velocity.
3. The flow is irrotational: :math:`\nabla \times u = 0`.
4. The wave amplitude is small with respect to the wavelength.
5. The amplitude of the body motion is small with respect to its dimension.
6. The sea bottom is flat. The water depth is denoted :math:`h`.

Mathematical problem
--------------------

The mass conservation equation can be rewritten as the Laplace equation

.. math::
   \nabla^2 \phi = 0,
   :label: laplace

where :math:`\phi` is the velocity potential defined as :math:`u = \nabla \phi`.

Since the problem is linear, we look for a solution in the frequency domain:

.. math::
   \phi = \mathrm{Re} \left( \Phi e^{-i \omega t} \right).


The partial differential equation :eq:`laplace` is completed with the following boundary conditions:

* linearized free surface boundary condition:

.. math::
   g \frac{\partial \Phi}{\partial z} - \omega^2 \Phi = 0, \qquad \text{on } z = 0.
   :label: bc_fs

* no velocity boundary condition at the (flat) sea bottom:

.. math::
   \frac{\partial \Phi}{\partial z} = 0, \qquad \text{on } z = -h.
   :label: bc_bottom

* a given velocity :math:`u` on the floating body surface :math:`\Gamma`:

.. math::
   \nabla \Phi \cdot n = u \cdot n, \qquad \text{on } \Gamma,
   :label: bc_body

where :math:`n` denotes the normal vector at the surface of the floating body.

.. * in the far field,
   .. math::
      \sqrt{R} \left( \frac{\partial \Phi}{\partial R} - i k \right) \left( \Phi - Phi_0 \right)
      \rightarrow 0, \qquad \text{when } R \rightarrow \infty,

The normal velocity on the floating body surface is the input of the problem.
It depends on the type of problem:

**Radiation problem**:
    For the radiation problem, the normal velocity on the body surface corresponds to the motion of the body along one of its degrees of freedom.
    The resolution of the Laplace problem allows to derive the added mass and the radiation damping associated with this degree of freedom (see also Post-processing_).

**Diffraction problem**:
    For the diffraction problem, the velocity on the floating body is given by the velocity of Airy's wave field.
    Once the problem has been solved, the linear Froude-Krylov force is computed by the integration of the pressure (:math:`p = i \rho \omega \Phi`) on the floating body (see also Post-processing_).

    The incoming Airy's wave fields is given by

    .. math::
       \Phi_0 = - i \frac{g}{\omega} \frac{\cosh (k (z+h))}{\cosh (k h)} e^{i k (x \cos \beta + y \sin \beta)}

    in finite depth, where the wave number :math:`k` is defined by the dispersion relation :math:`\omega^2 = k g \tanh (k h)`, and by

    .. math::
       \Phi_0 = - i \frac{g}{\omega} e^{k z} e^{i k (x \cos \beta + y \sin \beta)}

    in infinite depth, where the wave number :math:`k` is defined by :math:`\omega^2 = k g`.

    In the above equations, :math:`\beta` is the angle of the incoming wave.
    The angle :math:`\beta = 0` corresponds to waves propagating in the :math:`x` direction from :math:`x=-\infty` to :math:`x=+\infty`.
    The angle :math:`\beta = \pi/2` corresponds to waves propagating in the :math:`y` direction from :math:`y=-\infty` to :math:`y=+\infty`.


Integral problem
----------------

The partial differential equation can be rewritten as a boundary integral problem.
Let us introduce the Green function :math:`G(\xi, \cdot)`, which is solution of the partial differential equation:

.. math::
   \nabla^2 G(x; \xi) = \delta(\xi - x), \qquad \forall x,

where the :math:`\nabla` is meant as the derivative with respect to :math:`x`.

The above equation is associated with the boundary condition :eq:`bc_fs` and :eq:`bc_bottom`, where :math:`\xi` is a given point in the domain and :math:`\delta` is the Dirac distribution.

With the help of this Green function :math:`G`, the potential of the surface of the floating body :math:`\Gamma` can be rewritten as a function of a source distribution :math:`\sigma`:

.. math::
   \Phi(x) = \iint_\Gamma \sigma(\xi) G(x; \xi) \, \mathrm{dS}(\xi)
   :label: potential_representation

for all point :math:`x` in the fluid or on the hull of the floating body :math:`\Gamma`.

.. note:: There is a typo in equation :eq:`potential_representation` in [BD15]_.

The integral on the other boundaries of the domain is zero due to the properties of the Green function.

The differentiation of :eq:`potential_representation` differs depending whether :math:`x` is in the bulk of the fluid or on the hull.

On the hull, one has [Del87]_:

.. math::
   \frac{\partial \Phi}{\partial n}(x) = (u \cdot n)(x) = \frac{\sigma(x)}{2} + \iint_\Gamma \sigma(\xi) \, (\nabla G(x; \xi) \cdot n) \, \mathrm{dS}(\xi).
   :label: normal_velocity_on_hull_representation

where :math:`x` is a point on :math:`\Gamma` and :math:`n` is the vector normal to :math:`\Gamma` in :math:`x`.
For any vector :math:`t` tangential to :math:`\Gamma` at :math:`x`, one has

.. math::
   \frac{\partial \Phi}{\partial t}(x) = (u \cdot t)(x) = \iint_\Gamma \sigma(\xi) \, (\nabla G(x; \xi) \cdot t) \, \mathrm{dS}(\xi).
   :label: tangential_velocity_on_hull_representation

Finally, for :math:`x` in the bulk of the fluid, one has

.. math::
   \nabla \Phi(x) = u(x) = \iint_\Gamma \sigma(\xi) \, \nabla G(x; \xi) \, \mathrm{dS}(\xi).
   :label: velocity_in_bulk_representation

.. note:: Dimensional analysis:

    :math:`\Phi` is in m²·s¯¹.

    :math:`\sigma` is in m·s¯¹.

    :math:`G` is in m¯¹.

Expression of the Green function
================================

In infinite depth
-----------------

The integral problem above relates the potential :math:`\Phi` to the normal velocity
:math:`u \cdot n` via the Green function :math:`G`. Let us know discuss the evaluation of this
function for an infinite water depth.
See also [X18]_.

Green function
~~~~~~~~~~~~~~

The infinite depth Green function takes the following form

.. math::
   G(\xi, x) = - \frac{1}{4 \pi} \left( \frac{1}{\|x - \xi\|} + k \mathcal{G}\left(k \sqrt{(x_1 - \xi_1)^2 + (x_2 - \xi_2)^2}, k (x_3 + \xi_3) \right) \right)
   :label: green_function_inf_depth

.. proof:property::

   The function :math:`G` is symmetric in the sense of

   .. math::

        \forall x, \xi, \quad G(x, \xi) = G(\xi, x).


The first term of :math:`G` is the usual Green function for the 3D Laplace equation without our specific boundary conditions.
The :math:`\mathcal{G}` term is complex-valued and it is introduced to satisfy the boundary conditions :eq:`bc_fs`.

Introducing the dimensionless variables :math:`r = k \sqrt{(\xi_1 - x_1)^2 + (\xi_2 - x_2)^2}` and :math:`z = k (x_3 + \xi_3)`, this term reads

.. math::
    \mathcal{G}(r, z) & = \frac{1}{\sqrt{r^2 + z^2}} + \frac{2}{\pi} \Re \left( \int^{\pi/2}_{-\pi/2}  J(\zeta(r, z, \theta)) \, \mathrm{d} \theta \right) \\
    & \qquad \qquad \qquad \qquad + 2 i \Re \left( \int^{\pi/2}_{-\pi/2} e^{\zeta (r, z, \theta)} \, \mathrm{d} \theta \right)
    :label: green_function_inf_depth_low_freq

where

.. math::
    J(\zeta) = e^\zeta \left[ E_1(\zeta) + i \pi \right]

where :math:`E_1` is the first exponential integral, defined as

.. math::
    E_1(\zeta) = \int_\zeta^\infty \frac{e^{-t}}{t} dt,

and

.. math::
    \zeta (r, z, \theta) = z + i r \cos \theta.
    :label: def_zeta

The first term of :eq:`green_function_inf_depth_low_freq` is actually a Rankine-type singularity similar to the first term of :eq:`green_function_inf_depth`, except that one of the point has been reflected through the free surface.

Variants of the formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. _integrate_one_over_zeta:

.. proof:lemma::

    The following identity holds [Del89]_:

    .. math::
       \Re \int^{\pi/2}_{-\pi/2} \frac{1}{\zeta(\theta)} \, \mathrm{d} \theta = - \frac{\pi}{\sqrt{r^2 + z^2}}.
       :label: int_1_over_zeta

The above lemma allows to retrieve the expression of the Green function found e.g. in [BD15]_:

.. math::
    \mathcal{G}(r, z) & = - \frac{1}{\sqrt{r^2 + z^2}} + \frac{2}{\pi} \Re \left( \int^{\pi/2}_{-\pi/2} \left( J(\zeta(r, z, \theta)) - \frac{1}{\zeta(r, z, \theta)} \right) \, \mathrm{d} \theta \right) \\
    & \qquad \qquad \qquad \qquad + 2 i \Re \left( \int^{\pi/2}_{-\pi/2} e^{\zeta (r, z, \theta)} \, \mathrm{d} \theta \right)
    :label: green_function_inf_depth_high_freq

(Note the minus sign in front of the first term.)

.. proof:lemma::

    The `zeroth order Bessel function of the first kind <https://personal.math.ubc.ca/~cbm/aands/page_360.htm>`_ :math:`J_0` and `the Struve function <https://personal.math.ubc.ca/~cbm/aands/page_496.htm>`_ :math:`H_0` are such that

    .. math::
        J_0(r) & = \frac{1}{\pi} \int_{-\pi/2}^{\pi/2} \cos(r\cos(\theta)) \, \mathrm{d} \theta \\
        H_0(r) & = \frac{1}{\pi} \int_{-\pi/2}^{\pi/2} \sin(r\cos(\theta)) \, \mathrm{d} \theta \\

    hence

    .. math::
        \int_{-\pi/2}^{\pi/2} i e^{\zeta} \, \mathrm{d} \theta = \pi e^z \left(- H_0(r) + i J_0(r) \right)


The function :math:`\mathcal{G}` can also be rewritten as

.. math::
    \mathcal{G}(r, z) & = \frac{1}{\sqrt{r^2 + z^2}} + \frac{2}{\pi} \int^{\pi/2}_{-\pi/2} \Re \left( e^\zeta E_1(\zeta) \right) \, \mathrm{d} \theta + 2 \int^{\pi/2}_{-\pi/2} i e^{\zeta (r, z, \theta)} \, \mathrm{d} \theta \\
    & = \frac{1}{\sqrt{r^2 + z^2}} + \frac{2}{\pi} \int^{\pi/2}_{-\pi/2} \Re \left( e^\zeta E_1(\zeta) \right) \, \mathrm{d} \theta + 2 \pi e^z \left( - H_0(r) + i J_0(r) \right)

Noblesse [N82]_ splits the function :math:`\mathcal{G}` into a near field term :math:`N` and a wave field :math:`W` such that

.. math::
   N(r, z) & = \frac{1}{\sqrt{r^2 + z^2}} + \frac{2}{\pi} \int^{\pi/2}_{-\pi/2} \Re \left( e^\zeta E_1(\zeta)  \right) \, \mathrm{d} \theta  \\
   W(r, z) & = 2 \pi e^z \left( - H_0(r) + i J_0(r) \right)


Note that :math:`E_1`, :math:`J_0` and :math:`H_0` are available for instance in the `Scipy library <https://docs.scipy.org/doc/scipy/reference/special.html>`_.


.. proof:lemma::

    For any function :math:`f`, the following two formulations of the integral are equivalent [Del89]_:

    .. math::
        \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} f \left(\zeta(\theta) \right) \mathrm{d} \theta =
        \int_{-\frac{\pi}{2}}^{\frac{\pi}{2}} f \left(\tilde{\zeta}(\theta) \right) \mathrm{d} \theta

    where :math:`\zeta` is defined in :eq:`def_zeta` and :math:`\tilde{\zeta}` is defined as

    .. math::
       \tilde{\zeta} (\theta) = k \left( x_3 + \xi_3 + i \left( (x_1 - \xi_1) \cos\theta + (x_2 - \xi_2) \sin\theta \right) \right).

.. proof:proof::

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

The gradient of the Green function can be written as

.. math::
   \nabla_x G(\xi, x) = - \frac{1}{4 \pi} \left( - \frac{x - \xi}{\|x - \xi\|^3} + k
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
   \frac{\partial z}{\partial x_3} & = k

and, using the identity :math:`J'(\zeta) = J(\zeta) - 1/\zeta`,

.. math::
   \frac{\partial \mathcal{G}}{\partial r} = & - \frac{r}{(r^2 + z^2)^{3/2}} + \frac{2}{\pi} \Re \left( \int_{-\pi/2}^{\pi/2} i \cos(\theta) \left( J(\zeta) - \frac{1}{\zeta} \right) \, \mathrm{d}\theta \right) \\
   & \qquad \qquad \qquad \qquad + 2 i \Re \left( \int^{\pi/2}_{-\pi/2} i \cos(\theta) e^{\zeta} \, \mathrm{d} \theta \right)

and

.. math::
   \frac{\partial \mathcal{G}}{\partial z} = & - \frac{z}{(r^2 + z^2)^{3/2}} + \frac{2}{\pi} \Re \left( \int_{-\pi/2}^{\pi/2} \left( J(\zeta) - \frac{1}{\zeta} \right) \, \mathrm{d}\theta \right) \\
    & \qquad \qquad \qquad \qquad + 2 i \Re \left( \int^{\pi/2}_{-\pi/2} e^{\zeta } \, \mathrm{d} \theta \right) \\

that is, using :numref:`Lemma {number} <integrate_one_over_zeta>`

.. math::
   \frac{\partial \mathcal{G}}{\partial z} = \mathcal{G}(r, z) + \frac{2}{\sqrt{r^2 + z^2}} - \frac{z}{(r^2 + z^2)^{3/2}}
   :label: green_function_inf_depth_dGdz


.. ..note:: There seems to be a typo in the term of :eq:`green_function_inf_depth_deriv_2` in [Del89]_ and [BD15]_.

.. note::
    The derivative of :math:`G` with respect to :math:`x_1` and :math:`x_2` are antisymmetric in the sense of

    .. math::
       :nowrap:

        \[
        \frac{\partial G}{\partial x_1} (\xi, x) = - \frac{\partial G}{\partial x_1}(x, \xi).
        \]

    Its derivative with respect to :math:`x_3` is symmetric in infinite depth.

    In finite depth, some terms of the derivative with respect to :math:`x_3` are symmetric and some are antisymmetric.



Higher order derivative
~~~~~~~~~~~~~~~~~~~~~~~

From :eq:`green_function_inf_depth_dGdz`, one has

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
    I_1(r, z) & = \frac{2}{\pi} \Re \left( \int^{\pi/2}_{-\pi/2} J(\zeta) \, \mathrm{d} \theta \right) \\
    I_2(r, z) & = \frac{2}{\pi} \Re \left( \int^{\pi/2}_{-\pi/2} \left( J(\zeta) - \frac{1}{\zeta} \right) \, \mathrm{d} \theta \right) \\
    I_3(r, z) & = 2 \Re \left( \int^{\pi/2}_{-\pi/2} e^{\zeta} \, \mathrm{d} \theta \right) \\
    I_4(r, z) & = \frac{2}{\pi} \Re \left( \int^{\pi/2}_{-\pi/2} i \cos(\theta) \left( J(\zeta) - \frac{1}{\zeta} \right) \, \mathrm{d} \theta \right) \\
    I_5(r, z) & = 2 \Re \left( \int^{\pi/2}_{-\pi/2} i \cos(\theta) e^{\zeta} \, \mathrm{d} \theta \right)

then :eq:`green_function_inf_depth_low_freq` reads

.. math::
   \mathcal{G}(r, z) = \frac{1}{\sqrt{r^2 + z^2}} + I_1(r, z) + i I_3(r, z).

and :eq:`green_function_inf_depth_high_freq` reads

.. math::
   \mathcal{G}(r, z) = \frac{-1}{\sqrt{r^2 + z^2}} + I_2(r, z) + i I_3(r, z).

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

As seen in :eq:`green_function_inf_depth_dGdz`, new reflected-Rankine-type
terms might appear in the derivative of the Green wave term.
By default, they are integrated with the same method used for the same
numerical quadrature method as the rest of the wave term.
The setting ``gf_singularities="low_freq_with_rankine_term"`` is an attempt to
integrate them exactly using the same code as the main reflected Rankine term.

In finite depth
---------------

Green function
~~~~~~~~~~~~~~

TODO

Gradient of the Green function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO


Symmetries
----------

The first term of :eq:`green_function_inf_depth` is invariant under all rotations and translations, whereas the other terms are invariant under isometric transformations that don't change the vertical coordinate (reflection across a vertical plane, rotation around a vertical axis, translation following an horizontal vector).


Discretization
==============

The equations :eq:`potential_representation` and :eq:`normal_velocity_on_hull_representation` can be discretized using a collocation method.
Considering a mesh of the surface of the floating body :math:`\Gamma = \cup_i \Gamma_i`:

.. math::
   \Phi_i   & = \Phi(x_i), \\
   \sigma_i & = \sigma(x_i), \\
   u_i      & = (u \cdot n)(x_i) \\
   S_{ij}   & = \iint_{\Gamma_j} G(x_i, \xi) \mathrm{dS}(\xi), \\
   K_{ij}   & = \frac{\delta_{ij}}{2} + \iint_{\Gamma_j} \nabla_x G(x_i; \xi) \cdot n_i \, \mathrm{dS}(\xi),

where for all :math:`i`, :math:`x_i` is the center of the face :math:`\Gamma_i` and :math:`n_i` is its normal vector.
Each element of the matrices :math:`S` and :math:`K` can be seen as the interaction between two faces of the mesh.

.. note::
   :math:`K` should not be confused with the similar matrix :math:`D` defined as:

   .. math::
      D_{ij} = \frac{\delta_{ij}}{2} + \iint_{\Gamma_j} \nabla_\xi G(x_i; \xi) \cdot n_j \, \mathrm{dS}(\xi).

   Note that the derivation of :math:`G` is done with respect to a different variable.

   The matrix :math:`D` is used in the `direct` boundary integral equation, as e.g. in HAMS [Liu19]_.
   In the mathematical literature, :math:`D` is also referred to as the `double layer operator` and :math:`K` as the `adjoint double layer operator`.


The matrices :math:`S` and :math:`K` relates the vectors :math:`\Phi`, :math:`u` and :math:`\sigma` through the following approximations of :eq:`potential_representation` and :eq:`normal_velocity_on_hull_representation`:

.. math::
   \Phi = S \sigma, \qquad u = K \sigma.
   :label: discrete_BEM_problem

The resolution of the discrete problem with Nemoh consists of two main steps:

1. The evaluation of the coefficients of the complex-valued matrices :math:`S` and :math:`K`
2. The resolution of the complex-valued linear problem :math:`K \sigma = u`.

Once :math:`\sigma` has been computed, :math:`\Phi` can be easily deduced.
Then other magnitudes such as the Froude-Krylov forces or the added mass can be derived.

.. mermaid::
    :caption: A simplified flowchart of the internals of Capytaine solver

    flowchart TD;
        h[Water depth] --> gf(Assembling matrices);
        ω[Wave frequency ω] --> gf(Assembling matrices);
        m[Mesh] --> gf;
        gf -- K matrix --> ls(Linear solver);
        un[Normal velocity on hull] --> ls;
        gf -- S matrix --> mvp(Matrix vector product);
        ls -- sources distribution σ --> mvp;
        mvp -- potential distribution Φ --> int("Integrate on mesh");
        m --> int;
        int --> f["Hydrodynamic forces\n(aka added mass and radiation damping)"]

        classDef input fill:#DDDDDD,color:#333333,stroke:#444444
        classDef step fill:#88BBBB,color:#003333,stroke:#226666
        classDef output fill:#FFE3AA,color:#553900,stroke:#AA8439
        class ω,m,un,h input
        class gf,ls,mvp,int step
        class f output


Problem with forward speed
==========================

We refer to [D22]_ for a detailed description of the theory behind the approximate forward speed model used in Capytaine.

It relies on the following hypotheses:

1. The magnitude :math:`U` of the forward speed is small.
2. The body is thin enough, such that the flow around the body assuming a rigid free surface (also called *double-body flow*) can be approximated by :math:`\overrightarrow{u} = (-U, 0, 0)` in the reference frame of the body.

Then, the following modification are done to the solver to take forward speed into account:

1. **Doppler shift:** The frequency used in the computation is replaced by the *encounter frequency*

.. math::
   \omega_e = \omega - k U \cos (\beta)

where :math:`k` is the wavenumber and :math:`\beta` is the wave direction.
For this purpose, the ``wave_direction`` parameter can be passed to radiation problem.

2. **Normal velocity on hull:** The boundary condition on the body radiating with a dof defined by the displacement :math:`\delta\!r(x, y, z)` reads

.. math::
   \frac{\partial \phi}{\partial n} = - i \omega_e \delta\!r \cdot n + U \frac{\partial \delta\! r}{\partial x} \cdot n

The above relationship has currently only been implemented for the six dofs of single rigid bodies, as follows

+-------+---------------------+------------------------------------------------+
| Dof   | :math:`\delta \! r` | :math:`\frac{\partial \delta\! r}{\partial x}` |
+=======+=====================+================================================+
| Surge | :math:`(1, 0, 0)`   | :math:`(0, 0, 0)`                              |
+-------+---------------------+------------------------------------------------+
| Sway  | :math:`(0, 1, 0)`   | :math:`(0, 0, 0)`                              |
+-------+---------------------+------------------------------------------------+
| Heave | :math:`(0, 0, 1)`   | :math:`(0, 0, 0)`                              |
+-------+---------------------+------------------------------------------------+
| Roll  | :math:`(0, -z, y)`  | :math:`(0, 0, 0)`                              |
+-------+---------------------+------------------------------------------------+
| Pitch | :math:`(-z, 0, x)`  | :math:`(0, 0, 1)`                              |
+-------+---------------------+------------------------------------------------+
| Yaw   | :math:`(y, -x, 0)`  | :math:`(0, -1, 0)`                             |
+-------+---------------------+------------------------------------------------+

In other words, the supplementary term is zero except for pitch and yaw.

3. **Gradient of potential in pressure:** The equation relating the potential to the pressure is updated as follows

.. math::
   p = -\rho \left( -i \omega_e \phi + U \frac{\partial \phi}{\partial x} \right)

Similarly the relationship between the potential and the free surface elevation reads

.. math::
   \eta = -\frac{1}{g} \left( -i \omega_e \phi + U \frac{\partial \phi}{\partial x} \right)


The computation of :math:`\frac{\partial \phi}{\partial x}` makes the problems with forward speed typically 50\% slower that problems without.



The overall workflow with forward speed thus looks as follows.

.. mermaid::
    :caption: A simplified flowchart of the internals of Capytaine solver **with forward speed**, where red boxes are the supplementary steps introduced by forward speed.

    graph TD
          h[Water depth] --> gf(Assembling matrices);
          omega[Wave frequency ω] --> doppler(Doppler shift);
          fs[Forward speed U] --> doppler(Doppler shift);
          fs --> un;
          doppler -- Encounter frequency --> gf(Assembling matrices)
          m[Mesh] --> gf;
          gf -- K matrix --> ls(Linear solver);
          dof[Degree of freedom] --> un(Normal velocity on hull);
          un --  RHS of linear problem --> ls;
          gf -- S matrix --> mvp(Matrix vector product);
          ls -- sources distribution σ --> mvp;
          ls -- sources distribution σ --> grad;
          gf -- extended K matrix --> grad;
          mvp -- potential distribution Φ --> int("Integrate pressure on mesh");
          grad(Matrix vector product) -- gradient of Φ --> int;
          m --> int;
          fs --> int;
          int --> f["Hydrodynamic forces"]

          classDef input fill:#DDDDDD,color:#333333,stroke:#444444
          classDef step fill:#88BBBB,color:#003333,stroke:#226666
          classDef newstep fill:#FFAAAA,color:#550000,stroke:#113939
          classDef output fill:#FFE3AA,color:#553900,stroke:#AA8439
          class fs,omega,m,h,dof input
          class doppler,un,grad newstep
          class gf,ls,mvp,int step
          class f output


Post-processing
===============

Forces on body surfaces
-----------------------

Forces acting on body surfaces are computed by integration of the pressure field.

.. math:: F_i = \int_\Gamma p(x) \, n(x) \cdot \delta\!r_i(x) \, dx = j \omega \rho \int_\Gamma \Phi(x) \, n(x) \cdot \delta\!r_i(x) \, dx

where :math:`p = j \omega \rho \Phi` stands for the complex-valued pressure fields in frequency-domain, :math:`n` is the normal vector on the hull :math:`\Gamma` (oriented towards the fluid in Capytaine, see :doc:`../user_manual/conventions`) and :math:`\delta\!r_i` is the local displacement of the hull of the degree of freedom :math:`i`.

For a single rigid body, the degrees of freedom reads:

+---------+----------------------------------------------------------------+
| Dof     | Local hull displacement                                        |
+=========+================================================================+
| Surge   | :math:`\delta\!r(x) = (1, 0, 0)`                               |
+---------+----------------------------------------------------------------+
| Sway    | :math:`\delta\!r(x) = (0, 1, 0)`                               |
+---------+----------------------------------------------------------------+
| Heave   | :math:`\delta\!r(x) = (0, 0, 1)`                               |
+---------+----------------------------------------------------------------+
| Roll    | :math:`\delta\!r(x) = (1, 0, 0) \times (x-x_0, y-y_0, z-z_0)`  |
+---------+----------------------------------------------------------------+
| Pitch   | :math:`\delta\!r(x) = (0, 1, 0) \times (x-x_0, y-y_0, z-z_0)`  |
+---------+----------------------------------------------------------------+
| Yaw     | :math:`\delta\!r(x) = (0, 0, 1) \times (x-x_0, y-y_0, z-z_0)`  |
+---------+----------------------------------------------------------------+

where :math:`(x_0, y_0, z_0)` is the rotation center and :math:`\times` denotes the cross product.


The potential field can be decomposed into three contributions, and so does the resulting force:

1. The Froude-Krylov forces :math:`F_{FK}`, from the integration of the
   incident wave field pressure (incoming plane waves). In Capytaine, the
   incident wave pressure can be retrieved with the
   :func:`~capytaine.bem.airy_waves.airy_wave_pressure` function.
2. The diffraction forces :math:`F_{D}`, from the integration of the diffracted
   wave field (all bodies held fixed).
3. The radiation forces :math:`F_{R}`, which is itself a linear combination of
   the forces exerted by the fluid on the body in response to a motion of each
   degree of freedom.

The component :math:`i` of the radiation force :math:`F_{R}` is further rewritten as

.. math:: F_{R, i} = \sum_k \left[\omega^2 A_{ik} + j \omega B_{ik}\right] X_k

where :math:`A_{ik}` is the added mass matrix, :math:`B_{ik}` is the radiation
damping matrix and :math:`X_k` is the amplitude of the motion of the body along
the degree of freedom :math:`k`.

In other words, one has

.. math::
   A_{ik} & = \frac{1}{\omega^2} \Re \left[ j \omega \rho \int_\Gamma \Phi_k(x) \, n(x) \cdot \delta \! r_i(x) \, dx \right] \\
          & = - \frac{\rho}{\omega} \int_\Gamma \Im [\Phi_k(x)] \, n(x) \cdot \delta \! r_i(x) \, dx

and

.. math::
   B_{ik} & = \frac{1}{\omega} \Im \left[ j \omega \rho \int_\Gamma \Phi_k(x) \, n(x) \cdot \delta \! r_i(x) \, dx \right] \\
          & = \rho \int_\Gamma \Re [\Phi_k(x)] \, n(x) \cdot \delta \! r_i(x) \, dx

where :math:`\Phi_k` is the potential field computed with the normal velocity on the hull :math:`\frac{\partial \Phi_k}{\partial n} = -j \omega \delta \! r_k \cdot n`.
In Capytaine's wording, the degree of freedom :math:`k` defining the normal velocity on the hull is called ``radiating_dof``, while the degree of freedom :math:`i` used in the integration of the force is the ``influenced_dof``.

.. note::
   From Green second identity

   .. math:: \int_\Gamma \left[ \Phi_i \frac{\partial \Phi_k}{\partial n} - \Phi_k \frac{\partial \Phi_i}{\partial n}\right] dx = 0

   one has, when using the definition of the normal velocity of the radiation problem above,

   .. math:: \iint_{\Gamma} \Phi_i \; \delta\!r_k \cdot n = \iint_{\Gamma} \Phi_k \; \delta\!r_i \cdot n

   from which we can deduce the symmetry of the added mass matrix and the radiation dampings matrix.


.. note::
   As an alternative to :math:`\frac{\partial \Phi_k}{\partial n} = -j \omega
   \delta \! r_k \cdot n`, some software such as the version 1 of Capytaine use
   :math:`\frac{\partial \tilde \Phi_k}{\partial n} = \delta \! r_k \cdot n`,
   that is :math:`\tilde \Phi_k = \frac{\Phi_k}{-j \omega}`.

   It leads to the following definition of the added mass and radiation damping

   .. math::
      A_{ik} & = \frac{1}{\omega^2} \Re \left[ j \omega \rho \int_\Gamma (- j \omega \tilde \Phi_k(x)) \, n(x) \cdot \delta \! r_j(x) \, dx \right] \\
             & = \rho \int_\Gamma \Re [\tilde \Phi_k(x)] \, n(x) \cdot \delta \! r_j(x) \, dx

   and

   .. math::
      B_{ik} & = \frac{1}{\omega} \Im \left[ j \omega \rho \int_\Gamma (- j \omega \tilde \Phi_k(x)) \, n(x) \cdot \delta \! r_j(x) \, dx \right] \\
             & = \rho \omega \int_\Gamma \Im [\tilde \Phi_k(x)] \, n(x) \cdot \delta \! r_j(x) \, dx

   This form is convenient since the all the :math:`\omega` in the expression
   of the added mass disappears, which make it possible to compute the value of
   the added mass at frequency such as zero or infinity.

   However, the implementation of :math:`\tilde \Phi` in version 1 of Capytaine
   was not consistent with the use of :math:`\Phi` for diffraction problem and
   it was easy to forget the missing :math:`-j\omega` for some post-processing
   of :math:`\tilde \Phi` for radiation problems.

   In version 2.0 of Capytaine, :math:`\Phi` is used everywhere instead of
   :math:`\tilde \Phi`.
   Since version 2.1, another method has been implemented to take into account
   the cancelling of the :math:`\omega` in the expression of the added mass
   allowing to compute the added mass at zero and infinite frequency.


Dynamic coupling and impedance
------------------------------
Consider a body or a system of bodies. The general linear equation of motion can be expressed in time domain as

.. math:: M_{ij} \ddot{x}_j + C_{ij} \dot{x}_j + K_{ij} x_j = F_i,

and in frequency domain, with the assumed time dependence :math:`x(t) = \mathrm{Re} \left( X e^{-j \omega t} \right)`,

.. math:: \left[-\omega^2M_{ij} - j \omega C_{ij} + K_{ij}\right] X_j = F_i,

where :math:`M_{ij}` is the inertia matrix, accounting for the mass distribution, :math:`C_{ij}` is the mechanical damping matrix, :math:`K_{ij}` is the stiffness matrix which comprises mechanical and hydrostatic effects, and :math:`F_i` are generic external forces.

.. note:: The hydrostatic contribution to matrix :math:`K_{ij}` accounts for a variation of hydrostatic force in direction :math:`i` due to a unit motion in direction :math:`j`. It is a geometric property of the body.

As seen above, forces :math:`F_i` can be decomposed as

.. math:: F_i = F_{FK, i} + F_{D, i} + F_{R, i}

The full system becomes

.. math:: \left[-\omega^2 (M_{ij} + A_{ij}) - j \omega (C_{ij} + B_{ij}) + K_{ij}\right] X_j = F_{FK, i} + F_{D, i}

that is

.. math:: H X = F_{ex}

where :math:`H` denotes the following transfer function matrix

.. math:: H_{ij} = \left[-\omega^2 (M_{ij} + A_{ij}) - j \omega (C_{ij} + B_{ij}) + K_{ij}\right]

and :math:`F_{ex}` denotes the excitation force.

.. math:: F_{ex, i} = F_{FK, i} + F_{D, i}.

The oscillation amplitude is obtained by solving the complex-valued linear system.

.. note:: Matrices :math:`A_{ij}` and :math:`B_{ij}` depend on :math:`\omega`, and so does :math:`H_{ij}` and :math:`X_j`.

Free surface elevation
----------------------

The potential at the reference surface :math:`z = 0` can be connected to the free surface elevation by the dynamic condition

.. math:: \dfrac{\partial \phi}{\partial t} = - g \eta

which, in frequency domain, is

.. math:: \eta = \dfrac{j \omega}{g} \Phi

For a fully coupled problem (bodies free to oscillate, i.e. diffraction and radiation combined), the free surface elevation can be computed as

.. math:: \eta = \eta_{\text{incident}} + \eta_{\text{diffracted}} + \sum_i \eta_{\text{radiated}, i} X_i.


Far-field coefficients
----------------------

TODO
