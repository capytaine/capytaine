==============================
Green function in finite depth
==============================

.. contents:: Contents

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

.. note::

   .. math::
      \frac{1}{|x - \xi|} = \frac{1}{\sqrt{r^2 + (x_3 - \xi_3)^2}}
      :label: rankine_term

   .. math::
       \frac{1}{|x - S_h(\xi)|} & = \frac{1}{\sqrt{r^2 + (x_3 - (- \xi_3 - 2h))^2}} \\
      & = \frac{1}{\sqrt{r^2 + ((- x_3 - 2h) - \xi_3)^2}} = \frac{1}{|S_h(x) - \xi|}
      :label: rankine_term_sea_bottom_reflection


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
   :label: green_function_fin_depth_after_cosh_expansion


Following Newman [N85]_, we notice that the Green function can be written as two calls to a single function of three variables:

.. math::
   -4 \pi G(x, \xi, \omega, h) = L(r, x_3 - \xi_3, h) + L(r, 2h + x_3 + \xi_3, h)

where

.. math::
   L(r, z, h) & = \frac{1}{\sqrt{r^2 + z^2}} + \mathcal{L}(r, z, h) \\
   \mathcal{L}(r, z, h) & = \frac{1}{2} \int_0^{\infty} F(\kappa) e^{- 2 \kappa h} \left( e^{\kappa z} + e^{- \kappa z } \right) J_0(\kappa r) \mathrm{d} \kappa \\
                        & = \frac{1}{2} \int_0^{\infty} F(\kappa) e^{\kappa(z - 2 h)} J_0(\kappa r) \mathrm{d} \kappa + \frac{1}{2} \int_0^{\infty} F(\kappa) e^{- \kappa(z + 2 h)} J_0(\kappa r) \mathrm{d} \kappa

Since :math:`x_3 \in [-h, 0]` and :math:`\xi_3 \in [-h, 0]`, we have :math:`x_3 - \xi_3 \in [-h, h]` and :math:`2h + x_3 + \xi_3 \in [0, 2h]`, hence :math:`\mathcal{L}` should be defined for :math:`z \in [-h, 2h]`.
We kept a factor :math:`e^{-2\kappa h}` outside of :math:`F`, such that the exponent of the exponential is negative, in order to use Lipschitz-Hankel integral later.

.. .. note::
   Alternative decompositions (to be completed)

   Chen [C93]_ writes instead the Green function as follows (up to the non-dimensionnalization of the parameters):

   .. math::
      -4 \pi G(x, \xi, \omega, h) = \frac{1}{|x - \xi|} + \frac{1}{|x - S_h(\xi)|} + G_M(r, x_3 - \xi_3, h) + G_P(r, x_3 + \xi_3, h)

   where in our notations

   .. math::
      G_M(r, z, h) & = \mathcal{L}(r, z, h) \\
      G_P(r, z, h) & = \mathcal{L}(r, 2h + z, h)

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

or

.. math::
   A(kh) = \frac{(1 + \tanh(kh))^2}{1 - \tanh(kh)^2 + \frac{\tanh(kh)}{kh}}

Revisiting :math:`\mathcal{L}`:

.. math::
   \mathcal{L}(r, z, h) & = \frac{1}{2} \int_0^{\infty} \left(\frac{A(kh)}{\kappa h - k h} + F_1(\kappa)\right) e^{- 2 \kappa h} \left( e^{\kappa z} + e^{- \kappa z } \right) J_0(\kappa r) \mathrm{d} \kappa \\
                        & = \frac{A(kh)}{2h} \int_0^{\infty} \frac{e^{- 2 \kappa h} \left( e^{\kappa z} + e^{- \kappa z } \right)}{\kappa - k} J_0(\kappa r) \mathrm{d} \kappa \\
                           & \qquad + \frac{1}{2} \int_0^\infty F_1(\kappa) e^{- 2 \kappa h} \left( e^{\kappa z} + e^{- \kappa z } \right) J_0(\kappa r) \mathrm{d} \kappa \\

Terms of the form :math:`\int_0^\infty \frac{e^{\kappa Z}}{\kappa - k} J_0(\kappa r) d\kappa` have been seen in :doc:`inf_depth_green_function` as infinite-depth free surface term, and the same numerical evaluation method can be used.

Namely,

.. math::
   & \frac{A(kh)}{2h} \int_0^{\infty} \frac{e^{- 2 \kappa h} \left( e^{\kappa z} + e^{- \kappa z } \right)}{\kappa - k} J_0(\kappa r) \mathrm{d} \kappa \\
   &  \qquad = \frac{A(kh)}{2h} \frac{1}{2} \left(\mathcal{G}^+(r, z - 2h) +  \mathcal{G}^+(r, -(z + 2h))\right)


For the evaluation of the remaining term involving :math:`F_1`, the strategy chosen by Guével, Daubisse and Delhommeau is to approximate it by Prony decomposition, that is as a sum of exponentials, in order to employ the other kind of integral that was evaluated in :doc:`inf_depth_green_function`, namely the Lipschitz integrals.



We look for a Prony decomposition of :math:`F_1` of the following form:

.. math::
   F_1(\kappa) \simeq 2 + \sum_{i=1}^N a_i e^{-\lambda_i h \kappa}

with :math:`\forall i \ge 1, \lambda_i > 0` to always recover :math:`\lim_{\kappa \rightarrow \infty} F_1(\kappa) = 2`.
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
   :label: green_function_fin_depth_full_expression


Notice it contains six Rankine terms that appears in all expressions of the finite depth Green function, :math:`4N` other Rankine terms due to the approximation we used for the residual term, and four infinite-depth wave terms.

The Rankine terms are real-valued. The four infinite-depth Green functions :math:`\mathcal{G}^+` are complex-valued.

.. note::
   *Interpreting the Rankine terms as reflections*

   Besides :eq:`rankine_term` and :eq:`rankine_term_sea_bottom_reflection`, the Green function contains the following Rankine terms:

   .. math::
      \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 - 2h)^2}} = \frac{1}{|x - S_0(S_h(\xi))|} = \frac{1}{|S_h(S_0(x)) - \xi|}

   .. math::
      \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 + 2h)^2}} = \frac{1}{|x - S_h(S_0(\xi))|} = \frac{1}{|S_0(S_h(x)) - \xi|}

   .. math::
      \frac{1}{\sqrt{r^2 + (x_3 + \xi_3)^2}} = \frac{1}{|x - S_0(\xi)|} = \frac{1}{|S_0(x) - \xi|}
      :label: rankine_term_free_surface_reflection

   .. math::
      \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 4h)^2}} & = \frac{1}{|x - S_h(S_0(S_h(\xi)))|} = \frac{1}{|S_h(S_0(S_h(x))) - \xi|} \\
       & = \frac{1}{|x - S_{2h}(\xi)|} = \frac{1}{|S_{2h}(x) - \xi|}

   but also

   .. math::
      \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 - 2h - \lambda_i h)^2}} = \frac{1}{|x - S_0(S_{(1 + \lambda_i/2)h}(\xi))|}

   .. math::
      \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 + 2h + \lambda_i h)^2}} = \frac{1}{|x - S_{(1+\lambda_i/2)h}(S_0(\xi))|}

   .. math::
      \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 - \lambda_i h)^2}} = \frac{1}{|x - S_{(-\lambda_i/2) h}(\xi)|}

   .. math::
      \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 4h + \lambda_i h)^2}} = \frac{1}{|x - S_{(2+\lambda_i/2)h}(\xi)|}


   such that most terms of the Green function can be interpreted as Rankine term with :math:`6 + 4N` reflexions of :math:`\xi` (or :math:`x`).


.. note::
   *Singularities*

   Since :math:`\lambda_i > 0`, all the reflected terms are outside of the fluid domain and do not need a special treatment when integrating the Green function on a panel to avoid a singularity.

   The only terms that should be taken care of are :eq:`rankine_term` (always), :eq:`rankine_term_sea_bottom_reflection` (when the panel is on the sea bottom) and :eq:`rankine_term_free_surface_reflection` (when the panel is on the free surface).

   Beside the singularities of Rankine term, the infinite-depth wave terms :math:`\mathcal{G}^+` have a logarithmic singularity, but only the :math:`\mathcal{G}^+(r, x_3 + \xi_3)` term may reach its singularity.

.. _prony-decomposition-evaluation:

Evaluation of the Prony decomposition
-------------------------------------

For each :math:`k` and :math:`h`, the coefficients :math:`a_i` and :math:`\lambda_i` need to be precomputed.

The function :math:`\kappa h \mapsto F(\kappa h) - \frac{A(kh)}{\kappa h - k h} - 2` is evaluated for a range of values of :math:`\kappa h` in :math:`[-0.1, 20]`.
Looking for an approximation as small as possible, we start by evaluating the function at a few reference points and fit a Prony decomposition with few exponentials. If the resulting accuracy computed with more reference evaluations of the function is not sufficient, the number of points and the number of exponentials is increased (until a maximal value of 30 exponential terms).

Because of the singularity, the expression above cannot be evaluated when :math:`\kappa = k`, even if the function is mathematically smooth at this point. In Nemoh and in the legacy Fortran code of Capytaine, the function is replaced by a polynomial interpolation in :math:`[k - \max(0.1, 0.1k), k + \max(0.1, 0.1k)]` to avoid this.
In practice we noticed that, in double floating point precision, the effect of the singularity on the computational accuracy is restricted to a much smaller interval. Instead of patching the function, we add some small noise on the choice of the reference points ensure that among the several tries with different resolutions, most of them will not hit the singularity.

A remaining issue with this method is the second singularity of :math:`F`: the function has another pole of opposite sign in :math:`\kappa = -k`. When :math:`kh` is small, this singularity interferes with the evaluation of the Prony decomposition (see also GH635_). For this reason, it is currently not possible to compute finite depth problems in Capytaine with :math:`kh < 0.1`.

.. _GH635: https://github.com/capytaine/capytaine/issues/635#issuecomment-2729948600

Asymptotics
-----------

Deep water asymptotics
~~~~~~~~~~~~~~~~~~~~~~

When :math:`h \rightarrow \infty`, we have

.. math::
   F(\kappa) \sim_{h \rightarrow \infty} 2 \frac{\kappa + k}{\kappa - k}

and :eq:`green_function_fin_depth_after_cosh_expansion` becomes:

.. math::
   -4 \pi G(x, \xi, \omega, h) = \frac{1}{|x - \xi|} + \int_0^\infty \frac{\kappa + k}{\kappa - k} e^{\kappa(x_3 + \xi_3)} J_0 \left(\kappa r \right) d \kappa

recovering the expression found in :doc:`inf_depth_green_function`.

Alternatively, showing that

.. math::
   A(kh) \sim_{h \rightarrow \infty} 4 kh

and

.. math::
   F_1(\kappa) \sim_{h \rightarrow \infty} 2

hence

.. math::
   \forall i \ge 1, \quad a_i(kh) \sim_{h \rightarrow \infty} 0

makes :eq:`green_function_fin_depth_full_expression` simplify as follows

.. math::
   -4 \pi G(x, \xi, h)
   & = \frac{1}{\sqrt{r^2 + (x_3-\xi_3)^2}} + \cancel{\frac{1}{\sqrt{r^2 + (x_3 - \xi_3 - 2h)^2}}} + \cancel{\frac{1}{\sqrt{r^2 + (x_3 - \xi_3 + 2h)^2}}} \\
   & \qquad + \cancel{\frac{A(kh)}{4h} \mathcal{G}^+(r, x_3 - \xi_3 - 2h)} + \cancel{\frac{A(kh)}{4h} \mathcal{G}^+(r, -(x_3 - \xi_3 + 2h))} \\
   &  \qquad \qquad + \sum_{i=1}^N \cancel{\frac{a_i}{2} \frac{1}{\sqrt{r^2 + (x_3 - \xi_3- 2h - \lambda_i h)^2}}}
      + \sum_{i=1}^N \cancel{\frac{a_i}{2} \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 + 2h + \lambda_i h)^2}}} \\
   & + \cancel{\frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 2h)^2}}} + \frac{1}{\sqrt{r^2 + (x_3 + \xi_3)^2}} + \cancel{\frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 4h)^2}}} \\
   & \qquad + \frac{A(kh)}{4h} \mathcal{G}^+(r, x_3 + \xi_3) + \cancel{\frac{A(kh)}{4h} \mathcal{G}^+(r, -(x_3 + \xi_3 + 4h))} \\
   &  \qquad \qquad + \sum_{i=1}^N \cancel{\frac{a_i}{2} \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 - \lambda_i h)^2}}}
      + \sum_{i=1}^N \frac{a_i}{2} \cancel{\frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 4h + \lambda_i h)^2}}}

hence

.. math::
   -4 \pi G(x, \xi, h)
   = \frac{1}{\sqrt{r^2 + (x_3-\xi_3)^2}} + \frac{1}{\sqrt{r^2 + (x_3 + \xi_3)^2}} + k \mathcal{G}^+(r, x_3 + \xi_3)

where we also used the fact that

.. math::
   \mathcal{G}^+(r, z) \rightarrow_{z\rightarrow -\infty} 0

Infinite frequency asymptotics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When :math:`k \rightarrow \infty`, we have

.. math::
   F(\kappa) \sim - (1 + \tanh(\kappa h)) = - 2 \frac{1}{1 + e^{-2 \kappa h}} = - 2 - 2 \sum_{i=1}^\infty (-1)^i e^{-2 i \kappa h}

There is no :math:`\kappa - k` singularity directly visible anymore at the high-frequency asymptotics, so no infinite-depth wave terms, only Rankine terms are left.

.. note::

   .. math::
      \lim_{k \rightarrow \infty} \lim_{\kappa \rightarrow \infty } F(\kappa) = 2 \neq
      \lim_{\kappa \rightarrow \infty} \lim_{k \rightarrow \infty } F(\kappa) = -2

   because the :math:`\kappa - k` singularity is still around and messing with us.

A Prony decomposition of :math:`F` can be written as

.. math::
   a_i = - 2 (-1)^i, \qquad \lambda_i = 2i

and the Green function reads

.. math::
   & -4 \pi G(x, \xi, h)  \\
   & = \frac{1}{\sqrt{r^2 + (x_3-\xi_3)^2}} - \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 - 2h)^2}} - \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 + 2h)^2}} \\
   &  \qquad \qquad - \sum_{i=1}^\infty (-1)^i \frac{1}{\sqrt{r^2 + (x_3 - \xi_3- 2(1+i)h)^2}}
      - \sum_{i=1}^\infty (-1)^i \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 + 2( 1 + i)h)^2}} \\
   & + \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 2h)^2}} - \frac{1}{\sqrt{r^2 + (x_3 + \xi_3)^2}} - \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 4h)^2}} \\
   &  \qquad \qquad - \sum_{i=1}^\infty (-1)^i \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 - 2 i h)^2}}
      - \sum_{i=1}^\infty (-1)^i \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 2 (2 + i) h )^2}}

or equivalently

.. .. math::
   -4 \pi G(x, \xi, h)
   & = \frac{1}{\sqrt{r^2 + (x_3-\xi_3)^2}} + \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 2h)^2}} \\
   &  \qquad \qquad - \sum_{i=0}^\infty (-1)^i \frac{1}{\sqrt{r^2 + (x_3 - \xi_3- 2(1+i)h)^2}}
      - \sum_{i=0}^\infty (-1)^i \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 + 2( 1 + i)h)^2}} \\
   &  \qquad \qquad - \sum_{i=0}^\infty (-1)^i \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 - 2 i h)^2}}
      - \sum_{i=0}^\infty (-1)^i \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 2 (2 + i) h )^2}}

.. math::
   & -4 \pi G(x, \xi, h) \\
   & = - \sum_{i=0}^\infty (-1)^i \frac{1}{\sqrt{r^2 + (x_3 - \xi_3- 2(1+i)h)^2}}
      - \sum_{i=0}^\infty (-1)^{i+1} \frac{1}{\sqrt{r^2 + (x_3 - \xi_3 + 2 i h)^2}} \\
   &  \qquad - \sum_{i=0}^\infty (-1)^i \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 - 2 i h)^2}}
      - \sum_{i=0}^\infty (-1)^{i+1} \frac{1}{\sqrt{r^2 + (x_3 + \xi_3 + 2 (1 + i) h )^2}}

That is an infinity of identical (up to the sign) reflected Rankine terms meant to enforce the boundary conditions

.. math::
   \Phi = 0, \text{ on } z = 0, \qquad \frac{\partial \Phi}{\partial z} = 0, \text{ on } z = -h.

Zero frequency asymptotics
~~~~~~~~~~~~~~~~~~~~~~~~~~

As discussed in :ref:`prony-decomposition-evaluation`, the function :math:`F` has a second singularity in :math:`\kappa+k` cancelling out the :math:`\kappa-k` singularity for :math:`k = 0`.

.. math::
   F(\kappa) \sim \frac{1 + \tanh(\kappa h)}{\tanh(\kappa h)} = 2 + 2 \sum_{i=1}^\infty e^{-2 i \kappa h}

which would result in a infinite sum of Rankine kernels, similarly to the infinite-frequency case, ensuring the zero-frequency boundary conditions:

.. math::
   \frac{\partial \Phi}{\partial z} = 0, \text{ on } z = 0, \qquad \frac{\partial \Phi}{\partial z} = 0, \text{ on } z = -h.

This is currently not implemented in Capytaine.

Gradient of the Green function
------------------------------

We already know from :doc:`inf_depth_green_function` how to compute the derivative of all the terms constituing the finite depth Green function.
