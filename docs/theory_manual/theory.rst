Nemoh theory manual
===================

.. warning::
   Work in progress...

Hypotheses
----------

1. the fluid is inviscid
2. the fluid is incompressible: :math:`\nabla \cdot u = 0` with :math:`u` the flow velocity. 
3. the flow is irrotational: :math:`\nabla \times u = 0`.
4. the wave amplitude is small with respect to the wavelength.
5. the amplitude of the body motion is small with respect to its dimension.
6. the sea bottom is flat. The water depth is denoted :math:`h`.

Linear boundary value problem
-----------------------------

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

Radiation problem:
    For the radiation problem, the normal velocity on the body surface corresponds to the motion of the body along one of its degree of freedom.
    The resolution of the Laplace problem allows to derive the added mass and the radiation damping associated with this degree of freedom.

Diffraction problem:
    For the diffraction problem, the velocity on the floating body is given by the velocity of Airy's wave field (see below).
    Once the problem has been solved, the linear Froude-Krylov force is computed by the integration of the pressure (:math:`p = i \rho \omega \Phi`) on the floating body.

The incoming wave fields is given by

.. math::
   \Phi_0 = - i \frac{g}{\omega} \frac{\cosh (m_0 (z+h))}{\cosh (m_0 h)} e^{i m_0 (x \cos \beta + y \sin \beta)}

in finite depth, where the wave number :math:`m_0` is defined by the dispersion relation :math:`\omega^2 = m_0 g \tanh (m_0 h)`, and by

.. math::
   \Phi_0 = - i \frac{g}{\omega} e^{k_0 z} e^{i k_0 (x \cos \beta + y \sin \beta)}

in infinite depth, where the wave number :math:`k_0` is defined by :math:`\omega^2 = k_0 g`.


Integral problem
----------------

The partial differential equation can be rewritten as a boundary integral problem.
Let us introduce the Green function :math:`G(\xi, \cdot)`, which is solution of the partial differential equation:

.. math::
   \nabla^2 G(\xi, x) = \delta(\xi - x), \qquad \forall x,

associated with the boundary condition :eq:`bc_fs` and :eq:`bc_bottom`, where :math:`\xi` is a given point in the domain and :math:`\delta` is the Dirac distribution.

With the help of this Green function :math:`G`, the potential of the surface of the floating body :math:`\Gamma` can be rewritten as a function of a source distribution :math:`\sigma`:

.. math::
   \Phi(x) = \iint_\Gamma \sigma(y) G(x, y) \, \mathrm{dS}(y).
   :label: continuous_source_formulation

The integral on the other boundaries of the domain is zero due to the properties of the Green function.

The differentiation of :eq:`continuous_source_formulation` leads to the following equation \cite{delhommeau_problemes_1987}:

.. math::
   (u \cdot n)(x) = \frac{\sigma(x)}{2} + \iint_\Gamma \sigma(y) \, (\nabla_x G(x, y) \cdot n) \, \mathrm{dS}(y).
   :label: diff_continuous_source_formulation

where :math:`n` is the normal vector on the floating body surface :math:`\Gamma`.

Expression for the Green function
---------------------------------

The Green function can be written under the following form 

.. math::
   G(\xi, x) = -\frac{1}{4\pi |x - \xi|} + g_{(h, \omega)}\left((x_1 - \xi_1)^2 +(x_2 - \xi_2)^2, x_3, \xi_3\right)
   :label: green_function

where :math:`x = (x_1, x_2, x_3)` and :math:`\xi = (\xi_1, \xi_2, \xi_3)`. (The coordinate with index :math:`3` is the vertical coordinate.) \cite{delhommeau_problemes_1987}

The first term is the classical 3D Green function for the Laplace equation.
The second term :math:`g_{(h, \omega)}` is complex-valued.
It is introduced to satisfy the boundary conditions :eq:`bc_fs` and :eq:`bc_bottom`.
It depends on the water depth :math:`h` and the wave frequency :math:`\omega`.
We refer to \cite{babarit_theoretical_2015, delhommeau_problemes_1987} for details on its derivation and its evaluation.
(Note that in these references :math:`G` is split into three to four terms.
For the sake of simplicity, several of these terms have been merged into $g$ in this paper.) 

The first term of :eq:`green_function` is invariant under all rotations and translations, whereas the second term is invariant under isometric transformations that don't change the vertical coordinate (reflection across a vertical plane, rotation around a vertical axis, translation following an horizontal vector).


Discretization
--------------

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

1. The evaluation of the coefficients of the complex-valued matrices $S$ and $V$
2. The resolution of the complex-valued linear problem :math:`\left( \frac{\mathbb{I}}{2} + V \right) \sigma = u`.

Once :math:`\sigma` has been computed, :math:`\Phi` can be easily deduced.
Then other magnitudes such as the Froude-Krylov forces or the added mass can be derived.
