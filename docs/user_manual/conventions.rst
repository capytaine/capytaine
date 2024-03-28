==========================================
Conventions and differences to other codes
==========================================

.. note::
    To be expanded.

Unlike most other codes, angles (such as the incoming wave direction) are given in radians in Capytaine.

With respect to WAMIT
---------------------

In this section, the index :math:`W` denotes a magnitude in `WAMIT <https://www.wamit.com/>`_ convention. Other magnitudes use Capytaine convention.

Time dependency
~~~~~~~~~~~~~~~

In Capytaine, the complex-valued amplitudes (phasors) are defined with the convention

.. math::
   x(t) = \Re ( X e^{-i \omega t}).
   :label: time_convention_in_capytaine

It is unlike WAMIT in which the following convention is used

.. math::
   x(t) = \Re ( X_W e^{+ i \omega t})

Incoming potential
~~~~~~~~~~~~~~~~~~

In deep water, the potential of the incoming undisturbed velocity field reads in Capytaine

.. math::
   \Phi_0 = -i \frac{g}{\omega} e^{k z} e^{i k (x \cos \beta + y \sin \beta)},
   :label: incoming_waves_in_capytaine

whereas in WAMIT (equation 15.4 of WAMIT Manual v7.2), it reads

.. math::
   \Phi_{0, W} = - \overline{\Phi_0} = i \frac{g}{\omega} e^{k z} e^{- i k (x \cos \beta + y \sin \beta)},

and similarly in finite depth.

It follows that the incoming velocity field from Capytaine :math:`u_0 = \nabla \Phi_0` is related to the incoming velocity field from WAMIT :math:`u_{0, W} = \nabla \Phi_{0, W}` as

.. math::
   u_0 = \overline{u_{0, W}}.

Then the corresponding Froude-Krylov force and diffraction force (also called scattering force in WAMIT), as well as their sum the excitation force, are complex-conjugate between Capytaine and WAMIT:

.. math::
   F_e = \overline{F_{e, W}}

Normal vector
~~~~~~~~~~~~~

In its theory manual and related publications, WAMIT considers a normal vector on the hull that is oriented towards the inside of the floating body.
In Capytaine, the normal vector on the hull is oriented towards the outside of the floating body.


Green function and source distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In its theory manual and related publications, WAMIT defines the Green function as

.. math::
   G_W(x, \xi) = \frac{1}{|x - \xi|} + ...

whereas in Capytaine, the Green function is defined as

.. math::
   G(x, \xi) = \frac{-1}{4\pi} \left( \frac{1}{|x - \xi|} + ... \right) = \frac{- G_W(x, \xi)}{4\pi}

Similarly, the source distribution follows a different convention than in Capytaine:

.. math::
   \sigma_W = -\frac{\sigma}{4 \pi}.

Together with the convention on the normal vector mentionned above, it results in a slightly different expression for the boundary integral expression.


With respect to Nemoh and Aquadyn
---------------------------------

Capytaine mostly follows the same conventions as `Nemoh <https://gitlab.com/lheea/Nemoh>`_, which are also the same as in Aquadyn.
The main exception is the phase angle of the excitation force in Nemoh and Capytaine is the opposite of the one in Aquadyn.


With respect to HAMS
--------------------

`HAMS <https://github.com/YingyiLiu/HAMS>`_ follows the same conventions :eq:`time_convention_in_capytaine` and :eq:`incoming_waves_in_capytaine` as Capytaine, but in its documentation follows the same convention as WAMIT for normal vectors and Green function.


With respect to OrcaWave
------------------------

OrcaWave follows the same conventions as WAMIT.


With respect to Hydrostar
-------------------------

From (Donatini et al., 2022) [D22]_: 

* Hydrostar provided phases as a phase lead, while Capytaine outputs the phase lag.

* The phase of the incident waves is set such that the maximum is reached at :math:`t=0` at different points: in Capytaine, the reference point is :math:`(0, 0, 0)`, whereas in Hydrostar it is the center of buoyancy of the floating body by default.
