==========================================================
Capytaine: a Python-based linear potential flow BEM solver
==========================================================

Capytaine is a Python package for the simulation of the interaction between water waves and floating bodies in frequency domain.

It is built around a full rewrite of the open source Boundary Element Method (BEM) solver Nemoh_ for the linear potential flow wave theory.

This documentation is for version |release| (released |today|).

The documentations of older versions are available there: `v2.1 <https://capytaine.github.io/v2.1>`_, `v2.2.1 <https://capytaine.github.io/v2.2.1>`_, `v2.3.1 <https://capytaine.github.io/v2.3.1>`_.

.. toctree::
   :maxdepth: 1

   features.rst

Private support, custom developments and training can be provided by `Mews Labs <https://www.mews-labs.com/>`_ (contact@mews-labs.com).

For free support, you can ask questions on the public `Github discussion page <https://github.com/capytaine/capytaine/discussions>`_ or as `Github issues <https://github.com/capytaine/capytaine/issues/>`_.
Please do not contact the developers directly by email, unless you are looking for private paid support.

Contributions are welcome!
Please report bugs and suggest improvements as `Github issues <https://github.com/capytaine/capytaine/issues/>`_.

.. raw:: html

    <div style="float: right;">
    <video src="_static/front_page_animation.webm" loop autoplay muted width="320" height="240">
    </video>
    </div>

Documentation
=============

.. toctree::
   :maxdepth: 2

   user_manual/index.rst

.. toctree::
   :maxdepth: 2

   examples/index.rst

.. toctree::
   :maxdepth: 2

   developer_manual/index.rst

.. toctree::
   :maxdepth: 2

   theory_manual/index.rst

.. toctree::
   :maxdepth: 1

   citing.rst
   changelog.rst

Source code
===========

Available on `Github <https://github.com/capytaine/capytaine>`_

License
=======

Capytaine is developed by Matthieu Ancellin with the welcome help of `several contributors <https://github.com/capytaine/capytaine/graphs/contributors>`_.

Since April 2022, the development of Capytaine is funded by the Alliance for Sustainable Energy, LLC, Managing and Operating Contractor for the National Renewable Energy Laboratory (NREL) for the U.S. Department of Energy.

Since April 2025, the development of Capytaine is also funded by `Mews Labs <https://www.mews-labs.com/>`_ and BPI France.

From April 2017 to March 2019, the development of Capytaine at University College Dublin (UCD) was funded by Science Foundation Ireland (SFI) under Marine Renewable Energy Ireland (MaREI), the SFI Centre for Marine Renewable Energy Research.

It is based on version 2 of `Nemoh <https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp>`_, which has been developed by Gérard Delhommeau, Aurélien Babarit et al., (École Centrale de Nantes) and was distributed under the Apache License 2.0.

Since version 3, Capytaine is licensed under the Apache License, Version 2.0.
You may obtain a copy of the License at  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

This documentation is licensed under the `Creative Commons Attribution-ShareAlike 4.0 International License`_ |CCBYSA|.

.. |CCBYSA| image:: https://i.creativecommons.org/l/by-sa/4.0/80x15.png
.. _`Creative Commons Attribution-ShareAlike 4.0 International License`: http://creativecommons.org/licenses/by-sa/4.0/

The `boat mesh`_ in the animation above is in the public domain.

.. _`boat mesh`: https://opengameart.org/content/low-poly-pirate-ship


.. Indices and tables
   ------------------
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
