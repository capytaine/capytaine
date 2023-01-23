==========================================
Conventions and differences to other codes
==========================================

.. note::
    To be expanded.

Unlike most other codes, angles (such as the incoming wave direction) are given in radians in Capytaine.

With respect to WAMIT
---------------------

Time dependancy
~~~~~~~~~~~~~~~

In Capytaine, the complex-valued amplitudes (phasors) are defined with the convention :math:`x(t) = \Re ( X e^{-i \omega t})`.
It is unlike WAMIT in which the convention :math:`x(t) = \Re ( X e^{+ i \omega t})` is used.

Incoming potential
~~~~~~~~~~~~~~~~~~

In deep water, the potential of the incoming undisturbed velocity field reads in Capytaine

.. math::
   \Phi_0 = -i \frac{g}{\omega} e^{k z} e^{i k (x \cos \beta + y \sin \beta)},

whereas in WAMIT (equation 15.4 of WAMIT Manual v7.2), it reads

.. math::
   \Phi_{0, W} = - \overline{\Phi_0} = i \frac{g}{\omega} e^{k z} e^{- i k (x \cos \beta + y \sin \beta),

and similarly in finite depth.

It follows that the incoming velocity field from Capytaine :math:`u_0 = \nabla \Phi_0` is related to the incoming velocity field from WAMIT :math:`u_{0, W} = \nabla \Phi_{0, W}` as

.. math::
   u_0 = \overline{u_{0, W}}


With respect to Nemoh and Aquadyn
---------------------------------

Capytaine mostly follows the same conventions as Nemoh, which are also the same as in Aquadyn.
The main exception is the phase angle of the excitation force in Nemoh and Capytaine is the opposite of the one in Aquadyn.
