==============================
Green function in finite depth
==============================

Green function's Hankel form
----------------------------

Denoting by :math:`h` the water depth, the Green function can be written as (eq (2.33) from [Del87]_)

.. math::
   & -4 \pi G(x, \xi, \omega, h) = \frac{1}{|x - \xi|} + \frac{1}{|x - S_h(\xi)|}  \\
                            & \quad + 2 \int_0^\infty F(\kappa) e^{- 2 \kappa h} \cosh (\kappa (x_3 + h)) \cosh(\kappa (\xi_3 + h)) J_0 \left(\kappa r \right) d \kappa
   :label: green_function_fin_depth

where :math:`r = \sqrt{(\xi_1 - x_1)^2 + (\xi_2 - x_2)^2}` and :math:`S_h` denotes the symmetry across the horizontal plane :math:`z=h`:

.. math::
   S_h(\xi) = (\xi_1, \xi_2, -\xi_3 - 2h).

.. Note that :math:`S_h(S_h(x)) = x` and :math:`| S_h(x) - S_h(\xi) | = | x - \xi |`.

And :math:`F` is defined as (eq. (2.36) from [Del87]_):

.. math::
   F(\kappa) = \frac{\left(\kappa + \frac{\omega^2}{g}\right) e^{\kappa h}}{\kappa \sinh(\kappa h) - \frac{\omega^2}{g} \cosh(\kappa h)}

As in infinite depth, the integrand of :eq:`green_function_fin_depth` is singular for :math:`\kappa = k`, where :math:`k` the wavenumber solution of :math:`\omega^2/g = k \tanh(kh)`.
Replacing :math:`\omega^2/g` in the definition of :math:`F`:

.. math::
   F(\kappa) & = \frac{(\kappa + k \tanh(kh)) e^{\kappa h}}{\kappa \sinh(\kappa h) - k \tanh(kh) \cosh(\kappa h)} \\
            & = \frac{(\kappa + k \tanh(kh))(1 + \tanh(\kappa h))}{\kappa \tanh(\kappa h) - k \tanh(kh)} \\


Reduced variables
-----------------

We have by definition of :math:`\cosh`:

.. math::
   & \cosh (\kappa (x_3 + h)) \cosh(\kappa (\xi_3 + h)) \\
   & \qquad = \frac{1}{4} \left( e^{\kappa(x_3 + h)} + e^{\kappa(- x_3 - h)}\right)\left(e^{\kappa(\xi_3 + h)} + e^{\kappa(- \xi_3 - h)} \right) \\
   & \qquad = \frac{1}{4} \left( e^{\kappa(x_3 + \xi_3 +2h)} + e^{\kappa(x_3 - \xi_3)} + e^{\kappa(- x_3 + \xi_3)} + e^{\kappa(-x_3 - \xi_3 - 2h)} \right)

hence

.. math::
   & -4 \pi G(x, \xi, \omega, h) = \frac{1}{|x - \xi|} + \frac{1}{|x - S_h(\xi)|}  \\
   & \quad + \frac{1}{2} \int_0^\infty F(\kappa) e^{- 2 \kappa h}
   \left( e^{\kappa(x_3 + \xi_3 +2h)} + e^{\kappa(x_3 - \xi_3)} + e^{\kappa(- x_3 + \xi_3)} + e^{\kappa(-x_3 - \xi_3 - 2h)} \right)
   J_0 \left(\kappa r \right) d \kappa

Following Newman, we notice that the Green function can be written as two calls to a single function of three variables:

.. math::
   -4 \pi G(x, \xi, \omega, h) = L(r, x_3 - \xi_3, h) + L(r, 2h + x_3 + \xi_3, h)

where

.. math::
   L(r, z, h) & = \frac{1}{\sqrt{r^2 + z^2}} + \mathcal{L}(r, z, h) \\
   \mathcal{L}(r, z, h) & = \frac{1}{2} \int_0^{\infty} F(\kappa) e^{- 2 \kappa h} \left( e^{\kappa z} + e^{- \kappa z } \right) J_0(\kappa r) \mathrm{d} \kappa \\
                        & = \frac{1}{2} \int_0^{\infty} F(\kappa) e^{\kappa(z - 2 h)} J_0(\kappa r) \mathrm{d} \kappa + \frac{1}{2} \int_0^{\infty} F(\kappa) e^{- \kappa(z + 2 h)} J_0(\kappa r) \mathrm{d} \kappa

Since :math:`x_3 \in [-h, 0]` and :math:`\xi_3 \in [-h, 0]`, we have :math:`x_3 - \xi_3 \in [-h, h]` and :math:`2h + x_3 + \xi_3 \in [0, 2h]`, hence :math:`\mathcal{L}` should be defined for :math:`z \in [-h, 2h]`.
We kept a factor :math:`e^{-2\kappa h}` outside of :math:`F`, such that the exponent of the exponential is negative, in order to use Lipschitz-Hankel integral later.

Singularities extraction
------------------------

We extract the simple :math:`\kappa = k` singularity in :math:`F`:

.. math::
   (\kappa - k) F(\kappa) & = \frac{\kappa - k}{\kappa \tanh(\kappa h) - k \tanh(kh)} (\kappa + k \tanh(kh))(1 + \tanh(\kappa h))\\
   & \longrightarrow_{\kappa \rightarrow k} \frac{1}{\tanh(kh) + kh(1 - \tanh^2(kh))} (k + k \tanh(kh))(1 + \tanh(kh)) \\
   & \qquad = \frac{k (1 + \tanh(kh))^2}{\tanh(kh) + kh(1 - \tanh^2(kh))}

such that

.. math::

   F(\kappa) & = \frac{1}{\kappa - k} \frac{k (1 + \tanh(kh))^2}{\tanh(kh) + kh(1 - \tanh^2(kh))} + F_1(\kappa) \\
             & = \frac{A(kh)}{\kappa h - k h} + F_1(\kappa)

with :math:`F_1` a bounded function (on :math:`[0, \infty)` for :math:`kh > 0`) and :math:`A` is defined as:

.. math::
   A(kh) = \frac{k h \left(1 + \tanh(k h) \right)^2}{\tanh(kh) + kh \left(1-\tanh^2(kh)\right)}

also written in [Del87]_ and in the code as

.. math::
   A(kh) = \frac{(k h)^2 +  \left( kh \tanh(k h) \right)^2}{kh \tanh(kh) + kh^2 - \left(kh\tanh(kh)\right)^2}

Revisiting :math:`\mathcal{L}`:

.. math::
   \mathcal{L}(r, z, h) & = \frac{1}{2} \int_0^{\infty} \left(\frac{A(kh)}{\kappa h - k h} + F_1(\kappa)\right) e^{- 2 \kappa h} \left( e^{\kappa z} + e^{- \kappa z } \right) J_0(\kappa r) \mathrm{d} \kappa \\
                        & = \frac{A(kh)}{2h} \int_0^{\infty} \frac{e^{- 2 \kappa h} \left( e^{\kappa z} + e^{- \kappa z } \right)}{\kappa - k} J_0(\kappa r) \mathrm{d} \kappa \\
                           & \qquad + \frac{1}{2} \int_0^\infty F_1(\kappa) e^{- 2 \kappa h} \left( e^{\kappa z} + e^{- \kappa z } \right) J_0(\kappa r) \mathrm{d} \kappa \\

Terms of the form :math:`\int_0^\infty \frac{e^\kappa Z}{\kappa - k} J_0(\kappa r) d\kappa` have been seen in :doc:`inf_depth_green_function` as infinite-depth free surface term, and the same numerical evaluation method can be used.

Namely,

.. math::
   & \frac{A(kh)}{2h} \int_0^{\infty} \frac{e^{- 2 \kappa h} \left( e^{\kappa z} + e^{- \kappa z } \right)}{\kappa - k} J_0(\kappa r) \mathrm{d} \kappa \\
   &  \qquad = \frac{A(kh)}{2h} \frac{1}{2} \left(\mathcal{G}^+(r, z - 2h) +  \mathcal{G}^+(r, -(z + 2h))\right)


For the evaluation of the remaining term involving :math:`F_1`, the strategy chosen by Guével, Daubisse and Delhommeau is to approximate it by Prony decomposition, that is as a sum of exponentials, in order to employ the other kind of integral that was evaluated in :doc:`inf_depth_green_function`, namely the Lipschitz integrals.

.. prf:lemma::
   Hankel transform of Rankine terms

   That is the application of the Lipschitz integral presented in :doc:`inf_depth_green_function`:

   .. math::
      \frac{1}{|x - \xi|} = \int_0^\infty e^{- \kappa |x_3 - \xi_3|} J_0(\kappa r) d\kappa

   .. math::
      \frac{1}{|x - S_h(\xi)|} = \frac{1}{|S_h(x) - \xi|} = \int_0^\infty e^{- \kappa |x_3 - (-\xi_3 - 2h)|} J_0(\kappa r) d\kappa

   .. math::
      \frac{1}{|x - S_0(\xi)|} = \frac{1}{|S_0(x) - \xi|} = \int_0^\infty e^{\kappa (x_3 - (-\xi_3))} J_0(\kappa r) d\kappa

   .. math::
      \frac{1}{|x - S_0(S_h(\xi))|} = \frac{1}{|S_h(S_0(x)) - \xi|} = \int_0^\infty e^{\kappa (x_3 - (\xi_3 + 2h))} J_0(\kappa r) d\kappa

   .. math::
      \frac{1}{|x - S_h(S_0(\xi))|} = \frac{1}{|S_0(S_h(x)) - \xi|} = \int_0^\infty e^{-\kappa (x_3 - (\xi_3 - 2h))} J_0(\kappa r) d\kappa

   .. math::
      & \frac{1}{|x - S_h(S_0(S_h(\xi)))|} = \frac{1}{|S_h(S_0(S_h(x))) - \xi|} \\
       & \qquad \qquad = \frac{1}{|x - S_{2h}(\xi)|} = \frac{1}{|S_{2h}(x) - \xi|} = \int_0^\infty e^{-\kappa (x_3 - (-\xi_3 - 4h))} J_0(\kappa r) d\kappa

   but also

   .. math::
      \forall \lambda > 0, \qquad \frac{1}{|x - S_0(S_{\lambda h}(\xi))|} = \int_0^\infty e^{\kappa (x_3 - (\xi_3 + 2 \lambda h))} J_0(\kappa r) d\kappa

   .. math::
      \forall \lambda > 0, \qquad \frac{1}{|x - S_{\lambda h}(\xi)|} = \int_0^\infty e^{\kappa (x_3 + \xi_3 + 2 \lambda h))} J_0(\kappa r) d\kappa


We look for a Prony decomposition of :math:`F_1` of the following form:

.. math::
   F_1(\kappa) \simeq 2 + \sum_{i=1}^N a_i e^{-\lambda_i h \kappa}

with :math:`\forall i > 1, \lambda_i > 0` to always recover :math:`\lim_{\kappa \rightarrow \infty} F_1(\kappa) = 2`.
Note that :math:`F_1` depends on the water depth and the frequency, and so does :math:`a_i` and :math:`\lambda_i`. They are recomputed for each new :math:`(k, h)`.
The number of terms :math:`N` is chosen big enough to ensure a good accuracy of the approximation.

Then,

.. math::
   & \frac{1}{2} \int_0^\infty F_1(\kappa) e^{- 2 \kappa h} \left( e^{\kappa z} + e^{- \kappa z } \right) J_0(\kappa r) \mathrm{d} \kappa \\
   & \qquad \qquad \simeq \int_0^\infty e^{\kappa (z - 2h)} J_0(kr) d\kappa + \int_0^\infty e^{-\kappa (z + 2h)} J_0(kr) d\kappa \\
   &  \qquad \qquad \qquad + \sum_{i=1}^N \frac{a_i}{2} \int_0^\infty e^{\kappa (z - 2h - \lambda_i h)} J_0(kr) d\kappa
   + \sum_{i=1}^N \frac{a_i}{2} \int_0^\infty e^{- \kappa (z + 2h + \lambda_i h)} J_0(kr) d\kappa \\
   & \qquad \qquad = \frac{1}{\sqrt{r^2 + (z - 2h)^2}} + \frac{1}{\sqrt{r^2 + (z + 2h)^2}} \\
   &  \qquad \qquad \qquad + \sum_{i=1}^N \frac{a_i}{2} \frac{1}{\sqrt{r^2 + (z- 2h - \lambda_i h)^2}}
      + \sum_{i=1}^N \frac{a_i}{2} \frac{1}{\sqrt{r^2 + (z + 2h + \lambda_i h)^2}}


So :math:`\mathcal{L}` reads:

.. math::
   \mathcal{L}(r, z, h) & \simeq \frac{A(kh)}{4h} \left(\mathcal{G}^+(r, z - 2h) +  \mathcal{G}^+(r, -(z + 2h))\right) \\
   & \qquad + \frac{1}{\sqrt{r^2 + (z - 2h)^2}} + \frac{1}{\sqrt{r^2 + (z + 2h)^2}} \\
   &  \qquad \qquad + \sum_{i=1}^N \frac{a_i}{2} \frac{1}{\sqrt{r^2 + (z- 2h - \lambda_i h)^2}}
      + \sum_{i=1}^N \frac{a_i}{2} \frac{1}{\sqrt{r^2 + (z + 2h + \lambda_i h)^2}}

And the full unabrigde Green function reads:

..
   -4 \pi G(x, \xi, h)
   & = L(r, x_3 - \xi_3, h) + L(r, 2h + x_3 + \xi_3, h) \\
   & = \frac{1}{\sqrt{r^2 + (x_3-\xi_3)^2}} + \mathcal{L}(r, x_3 - \xi_3) \\
   & \qquad + \frac{1}{\sqrt{r^2 + (2h + x_3 + \xi_3)^2}} + \mathcal{L}(r, 2h + x_3 + \xi_3) \\

.. math::
   -4 \pi G(x, \xi, h)
   & = \frac{1}{\sqrt{r^2 + (x_3-\xi_3)^2}} + \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 - 2h)^2}} + \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 + 2h)^2}} \\
   & \qquad + \frac{A(kh)}{4h} \mathcal{G}^+(r, x_3 - \xi_3 - 2h) + \frac{A(kh)}{4h} \mathcal{G}^+(r, -(x_3 - \xi_3 + 2h)) \\
   &  \qquad \qquad + \sum_{i=1}^N \frac{a_i}{2} \frac{1}{\sqrt{r^2 + (x_3 - \xi_3- 2h - \lambda_i h)^2}}
      + \sum_{i=1}^N \frac{a_i}{2} \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 + 2h + \lambda_i h)^2}} \\
   & + \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 2h)^2}} + \frac{1}{\sqrt{r^2 + (x_3 + \xi_3)^2}} + \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 4h)^2}} \\
   & \qquad + \frac{A(kh)}{4h} \mathcal{G}^+(r, x_3 + \xi_3) + \frac{A(kh)}{4h} \mathcal{G}^+(r, -(x_3 + \xi_3 + 4h)) \\
   &  \qquad \qquad + \sum_{i=1}^N \frac{a_i}{2} \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 - \lambda_i h)^2}}
      + \sum_{i=1}^N \frac{a_i}{2} \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 4h + \lambda_i h)^2}}


Notice it contains six Rankine terms that appears in all expressions of the finite depth Green function, :math:`4N` other Rankine terms due to the approximation we used for the residual term, and four infinite-depth wave terms.
