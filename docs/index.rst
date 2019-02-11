===============================================
Capytaine: a Python-based distribution of Nemoh
===============================================

Capytaine is a Python package for the simulation of the interaction between water waves and floating bodies in frequency domain.

It is build around a full rewrite of the open source Boundary Element Method (BEM) solver Nemoh_ for the linear potential flow wave theory.

Latest release: |release| (|today|).

.. toctree::
   :maxdepth: 1

   features.rst

Contributions are welcome!
Please report bugs and suggest improvements as `Github issues <https://github.com/mancellin/capytaine/issues/>`_.

.. raw:: html

    <div style="float: right;">
    <video src="_static/front_page_animation.ogv" loop autoplay width="320" height="240">
    </video>
    </div>

Documentation
=============

.. toctree::
   :maxdepth: 2

   user_manual/index.rst

.. toctree::
   :maxdepth: 2

   developer_manual//index.rst

.. toctree::
   :maxdepth: 2

   theory_manual/index.rst

.. toctree::
   :maxdepth: 1

   citing.rst
   changelog.rst

Source code
===========

Available on `Github <https://github.com/mancellin/capytaine>`_

License
=======

.. raw:: html

    <div style="float: right;">
    <a href="http://www.ucd.ie/">
    <img src="_static/UCD-logo-nobg-mini.png" height="100" alt="UCD">
    </a>
    <a href="http://www.marei.ie/">
    <img src="_static/MaREI-Centre-logo-mini.png" height="100" alt="MaREI">
    </a>
    </div>

Capytaine has been developped by `Matthieu Ancellin`_ (University College Dublin) and is distributed under the terms of the GNU General Public License (GPL) v3.0. See the ``LICENSE`` file in the `code repository`_.

.. _`code repository`: https://github.com/mancellin/capytaine
.. _`Matthieu Ancellin`: http://ancell.in

This documentation is licensed under the `Creative Commons Attribution-ShareAlike 4.0 International License`_ |CCBYSA|.

.. |CCBYSA| image:: https://i.creativecommons.org/l/by-sa/4.0/80x15.png
.. _`Creative Commons Attribution-ShareAlike 4.0 International License`: http://creativecommons.org/licenses/by-sa/4.0/

Capytaine is a fork of Nemoh_, which has been developed by Gérard Delhommeau, Aurélien Babarit et al., (École Centrale de Nantes) and is distributed under the Apache License 2.0.

Capytaine includes code from meshmagick_ by François Rongère (École Centrale de Nantes), licensed under the GNU General Public License (GPL).

.. _meshmagick: https://github.com/LHEEA/meshmagick
.. _Nemoh: https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp

The `boat mesh`_ in the animation above is in the public domain.

.. _`boat mesh`: https://opengameart.org/content/low-poly-pirate-ship

.. Indices and tables
   ------------------
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
