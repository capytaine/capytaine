======================
Installation for users
======================

Capytaine is available on Windows, MacOS [#]_ and Linux.

.. [#] For the latest informations on the arm64 architectures (Apple M1), see https://github.com/capytaine/capytaine/issues/190

Capytaine requires Python 3.6 or higher.
Thus it is compatible with `all currently supported version of Python <https://devguide.python.org/versions/>`_.

With Conda
----------

The easiest way to install Capytaine is the precompiled package available on Conda_.
Download and install the `Anaconda distribution`_ or its lightweight counterparts Miniconda_ and Miniforge_.

.. _Conda: https://conda.io
.. _`Anaconda distribution`: https://www.anaconda.com/download/
.. _Miniconda: https://conda.io/miniconda.html
.. _Miniforge: https://github.com/conda-forge/miniforge

Once Conda has been installed, you might want to `create a dedicated environment <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`_.
Capytaine's package is available in the `conda-forge` channel and can be installed with the following command ::

    conda install -c conda-forge capytaine

The required dependencies should be installed automatically.


Optional dependencies
---------------------

Optional dependencies can be manually installed.
They are nice to have but not necessary for Capytaine's main features.

+------------+------------------------------------------+------------------------------+
| Name       | Example installation command             | Usage                        |
+============+==========================================+==============================+
| matplotlib | :code:`conda install matplotlib`         | Used in several examples     |
|            |                                          | in the documentation and     |
|            |                                          | the cookbook                 |
+------------+------------------------------------------+------------------------------+
| vtk        | :code:`conda install -c conda-forge vtk` | For 3D visualization         |
+------------+------------------------------------------+------------------------------+
| joblib     | :code:`conda install joblib`             | For parallel resolution      |
+------------+------------------------------------------+------------------------------+
| quadpy     | :code:`pip install quadpy`               | For higher order quadratures |
|            |                                          | (experimental)               |
+------------+------------------------------------------+------------------------------+


With Pip
--------

The package is available on PyPI, although only as a source distribution.
That means that you'll need a Fortran compiler [#]_ in order to install the package.
If you do, you can install Capytaine as::

    pip install capytaine

If you can't install a compiler, it is recommended to use Conda instead.

.. [#] For example, on Ubuntu or Debian: :code:`sudo apt install gfortran`.
       On macOS, see `for instance these instructions <https://github.com/capytaine/capytaine/issues/115#issuecomment-1143987636>`_.

