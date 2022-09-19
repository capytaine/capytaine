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
      \sqrt{R} \left( \frac{\partial \Phi}{\partial R} - i m_0 \right) \left( \Phi - Phi_0 \right)
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
       \Phi_0 = - i \frac{g}{\omega} \frac{\cosh (m_0 (z+h))}{\cosh (m_0 h)} e^{i m_0 (x \cos \beta + y \sin \beta)}

    in finite depth, where the wave number :math:`m_0` is defined by the dispersion relation :math:`\omega^2 = m_0 g \tanh (m_0 h)`, and by

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
   \nabla^2_x G(\xi, x) = \delta(\xi - x), \qquad \forall x,

associated with the boundary condition :eq:`bc_fs` and :eq:`bc_bottom`, where :math:`\xi` is a given point in the domain and :math:`\delta` is the Dirac distribution.

With the help of this Green function :math:`G`, the potential of the surface of the floating body :math:`\Gamma` can be rewritten as a function of a source distribution :math:`\sigma`:

.. math::
   \Phi(x) = \iint_\Gamma \sigma(y) G(x, y) \, \mathrm{dS}(y).
   :label: continuous_source_formulation

.. note:: There is a typo in this equation in [BD15]_.

The integral on the other boundaries of the domain is zero due to the properties of the Green function.

The differentiation of :eq:`continuous_source_formulation` leads to the following equation [Del87]_:

.. math::
   (u \cdot n)(x) = \frac{\sigma(x)}{2} + \iint_\Gamma \sigma(y) \, (\nabla_x G(x, y) \cdot n) \, \mathrm{dS}(y).
   :label: diff_continuous_source_formulation

where :math:`n` is the normal vector on the floating body surface :math:`\Gamma`.

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

Intoducting the dimensionless variables :math:`r = k \sqrt{(\xi_1 - x_1)^2 + (\xi_2 - x_2)^2}` and :math:`z = k (x_3 + \xi_3)`, this term reads

.. math::
    \mathcal{G}(r, z) & = \frac{1}{\sqrt{r^2 + z^2}} + \frac{2}{\pi} \Re \left( \int^{\pi/2}_{-\pi/2}  J(\zeta(r, z, \theta)) \, \mathrm{d} \theta \right) \\
    & \qquad \qquad \qquad \qquad + 2 i \Re \left( \int^{\pi/2}_{-\pi/2} e^{\zeta (r, z, \theta)} \, \mathrm{d} \theta \right)
    :label: green_function_inf_depth_xie

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

The first term of :eq:`green_function_inf_depth_xie` is actually a Rankine-type singularity similar to the first term of :eq:`green_function_inf_depth`, except that one of the point has been reflected through the free surface.

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
    & \qquad \qquad \qquad \qquad + 2 i k \Re \left( \int^{\pi/2}_{-\pi/2} e^{\zeta (r, z, \theta)} \, \mathrm{d} \theta \right)
    :label: green_function_inf_depth_del

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
   \frac{\partial r}{\partial x_1} & = k \frac{x_1 - \xi_1}{r} \\
   \frac{\partial r}{\partial x_2} & = k \frac{x_2 - \xi_2}{r} \\
   \frac{\partial z}{\partial x_3} & = k

and, using the identity :math:`J'(\zeta) = J(\zeta) - 1/\zeta`,

.. math::
   \frac{\partial \mathcal{G}}{\partial r} = & - \frac{r}{(r^2 + z^2)^{3/2}} + \frac{2}{\pi} \Re \left( \int_{-\pi/2}^{\pi/2} i \cos(\theta) \left( J(\zeta) - \frac{1}{\zeta} \right) \, \mathrm{d}\theta \right) \\
   & \qquad \qquad \qquad \qquad + 2 i \Re \left( \int^{\pi/2}_{-\pi/2} i \cos(\theta) e^{\zeta} \, \mathrm{d} \theta \right)

and

.. math::
   \frac{\partial \mathcal{G}}{\partial z} = & - \frac{z}{(r^2 + z^2)^{3/2}} + \frac{2}{\pi} \Re \left( \int_{-\pi/2}^{\pi/2} \left( J(\zeta) - \frac{1}{\zeta} \right) \, \mathrm{d}\theta \right) \\
    & \qquad \qquad \qquad \qquad + 2 i \Re \left( \int^{\pi/2}_{-\pi/2} i \cos(\theta) e^{\zeta } \, \mathrm{d} \theta \right) \\

that is, using :numref:`Lemma {number} <integrate_one_over_zeta>`

.. math::
   \frac{\partial \mathcal{G}}{\partial z} = \mathcal{G}(r, z) + \frac{1}{\sqrt{r^2 + z^2}} - \frac{z}{(r^2 + z^2)^{3/2}}
   :label: green_function_inf_depth_dGdz


.. ..note:: There seems to be a typo in the term of :eq:`green_function_inf_depth_deriv_2` in [Del89]_ and [BD15]_.

.. note::
    The derivative of :math:`G` with respect to :math:`x_1` and :math:`x_2` are antisymmetric in the sense of

    .. math::
       :nowrap:

        \[
        \frac{\partial G}{\partial x_1} (\xi, x) = - \frac{\partial G}{\partial x_1}(x, \xi).
        \]

    Its derivative with respect to :math:`x_3` can be decomposed into an antisymmetric term and a symmetric term.


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

Delhommeau's method is based on expression :eq:`green_function_inf_depth_del` of the Green function.
This expression of the Green function and its derivative require the evaluation of the following real-valued integrals:

.. math::
    D_1(r, z) & = \Re \left( \int^{\pi/2}_{-\pi/2} - i \cos(\theta) \left( J(\zeta) - \frac{1}{\zeta} \right) \, \mathrm{d} \theta \right) \\
    D_2(r, z) & = \Re \left( \int^{\pi/2}_{-\pi/2} - i \cos(\theta) e^{\zeta} \, \mathrm{d} \theta \right) \\
    Z_1(r, z) & = \Re \left( \int^{\pi/2}_{-\pi/2} \left( J(\zeta) - \frac{1}{\zeta} \right) \, \mathrm{d} \theta \right) \\
    Z_2(r, z) & = \Re \left( \int^{\pi/2}_{-\pi/2} e^{\zeta} \, \mathrm{d} \theta \right)


To limit the computational cost of the evaluation of these integrals, they are precomputed for selected values of :math:`r` and :math:`z` and stored in a table.
When evaluating the Green function, the values of the integrals are retrieved by interpolating the values in the tables.

For large values of :math:`r` and :math:`z`, these integrals are asymptotically approximated by the following expressions:

.. math::
      D_1(r, z) & \simeq \pi e^z \sqrt{\frac{2\pi}{r}} \left(\cos(r - \pi/4) - \frac{1}{2r} \sin(r-\pi/4) \right) - \pi \frac{r}{(r^2 + z^2)^{3/2}} \\
      D_2(r, z) & \simeq e^z \sqrt{\frac{2\pi}{r}} \left( \sin(r - \pi/4) + \frac{1}{2r} \cos(r - \pi/4) \right) \\
      Z_1(r, z) & \simeq - \pi e^z \sqrt{\frac{2\pi}{r}} \sin(r - \pi/4) + \pi \frac{z}{(r^2 + z^2)^{3/2}} \\
      Z_2(r, z) & \simeq e^z \sqrt{\frac{2\pi}{r}} \cos(r - \pi/4)


Incorporating these asymptotic approximation in the expression of the Green function, one gets:

.. math::
    \mathcal{G}(r, z) \simeq & -\frac{1}{\sqrt{r^2 + z^2}} - 2 k e^z \sqrt{\frac{2\pi}{r}} \left(\sin(r - \pi/4) - i\cos(r - \pi/4)\right) \\
   & \qquad\qquad\qquad\qquad + 2 k \frac{z}{(r^2 + z^2)^{3/2}}


Xie's variant
~~~~~~~~~~~~~

A slight variant is presented in [X18]_. The authors noticed that the
interpolation of the integral :math:`Z_1` can be inaccurate due to the
singularity :math:`\frac{1}{\zeta}`.
Hence, they proposed to use :eq:`green_function_inf_depth_xie` and to tabulate the integral

.. math::
    \widetilde{Z_1}(r, z) = \Re \left( \int^{\pi/2}_{-\pi/2} J(\zeta) \, \mathrm{d} \theta \right)

By using :numref:`Lemma {number} <integrate_one_over_zeta>`, one has

.. math::
   Z_1 = \widetilde{Z_1} - \frac{\pi}{\sqrt{r^2 + z^2}}

Both the original Delhommeau's method and Xie's variant are implemented in Capytaine.

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

The equations :eq:`continuous_source_formulation` and :eq:`diff_continuous_source_formulation` can be discretized using a collocation method.
Considering a mesh of the surface of the floating body :math:`\Gamma = \cup_i \Gamma_i`:

.. math::
   \Phi_i   & = \Phi(x_i), \\
   \sigma_i & = \sigma(x_i), \\
   u_i      & = (u \cdot n)(x_i) \\
   S_{ij}   & = \iint_{\Gamma_j} G(x_i, y) \mathrm{dS}(y), \\
   V_{ij}   & = \iint_{\Gamma_j} \nabla_{x_i} G(x_i, y) \cdot n \, \mathrm{dS}(y),

where for all :math:`i`, :math:`x_i` is the center of the face :math:`\Gamma_i`.
Each element of the matrices :math:`S` and :math:`V` can be seen as the interaction between two faces of the mesh.

The matrices :math:`S` and :math:`V` relates the vectors :math:`\Phi`, :math:`u` and :math:`\sigma` through the following approximations of :eq:`continuous_source_formulation` and :eq:`diff_continuous_source_formulation`:

.. math::
   \Phi = S \sigma, \qquad u = \left( \frac{\mathbb{I}}{2} + V \right) \sigma.
   :label: discrete_BEM_problem

The resolution of the discrete problem with Nemoh consists of two main steps:

1. The evaluation of the coefficients of the complex-valued matrices :math:`S` and :math:`V`
2. The resolution of the complex-valued linear problem :math:`K \sigma = \left( \frac{\mathbb{I}}{2} + V \right) \sigma = u`.

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

        classDef input fill:#FFAAAA,color:#550000,stroke:#113939
        classDef step fill:#88BBBB,color:#003333,stroke:#226666
        classDef output fill:#FFE3AA,color:#553900,stroke:#AA8439
        class ω,m,un,h input
        class gf,ls,mvp,int step
        class f output


Post-processing
===============

Forces on body surfaces
-----------------------

Forces acting on body surfaces are computed by integration of the pressure field. They can be decomposed into three contributions:

1. The Froude-Krylov forces :math:`F_{FK, i}`, from the integration of the incident wave field pressure (incoming plane waves); :math:`i` denotes the i-th degree of freedom
2. The diffraction forces :math:`F_{D, i}`, from the integration of the diffracted wave field (all bodies held fixed)
3. The radiation forces :math:`F_{R, ij}`, from the result of the radiation problem with radiating degree of freedom :math:`j` and influenced degree of freedom :math:`i`

Dynamic coupling and impedance
------------------------------
Consider a body or a system of bodies. The general linear equation of motion can be expressed in time domain as 

.. math:: M_{ij} \ddot{x}_j + C_{ij} \dot{x}_j + K_{ij} x_j = F_i,

and in frequency domain, with the assumed time dependence :math:`x(t) = \mathrm{Re} \left( X e^{-j \omega t} \right)`,

.. math:: \left[-\omega^2M_{ij} - j \omega C_{ij} + K_{ij}\right] X_j = F_i,

where :math:`M_{ij}` is the inertia matrix, accounting for the mass distribution, :math:`C_{ij}` is the mechanical damping matrix, :math:`K_{ij}` is the stiffness matrix which comprises mechanical and hydrostatic effects, and :math:`F_i` are generic external forces.

.. note:: The hydrostatic contribution to matrix :math:`K_{ij}` accounts for a variation of hydrostatic force in direction :math:`i` due to a unit motion in direction :math:`j`. It is a geometric property of the body.

Forces :math:`F_i` can be decomposed as

.. math:: F_i = F_{FK, i} + F_{D, i} + F_{R, ij}

and :math:`F_{R, ij}` can be further rewritten as 

.. math:: F_{R, ij} = \left[\omega^2 A_{ij} + j\omega B_{ij}\right] X_j

where :math:`A_{ij}` is the added mass matrix and :math:`B_{ij}` is the radiation damping matrix; these properties are thus obtained from the real and imaginary parts of the radiation force. The full system becomes

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

.. math:: \eta = \eta_{\text{incident}} + \eta_{\text{diffracted}} -j \omega \sum_i \eta_{\text{radiated}, i}   X_i

where factor :math:`-j \omega` transforms :math:`\eta_{\text{radiated}, i}` from the radiated wave field corresponding to unit oscillation velocity to the field corresponding to unit oscillation amplitude.


Far-field coefficients
----------------------

TODO

