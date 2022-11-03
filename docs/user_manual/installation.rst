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
| meshio     | :code:`pip install meshio`               | To load more mesh formats    |
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


With Docker
-----------

The following command will create a Docker image called `my_capytaine_image` based on Ubuntu 22.04 with the version 1.5 of Capytaine::

    docker build -t my_capytaine_image https://github.com/capytaine/capytaine.git#v1.5

Use for instance the following command to open an IPython shell in which Capytaine can be imported::

    docker run -it my_capytaine_image ipython3

Or the following command to make the current directory accessible from the Docker image and run the file :code:`my_script.py`::

    docker run -it -v $(pwd):/home/user my_capytaine_image python3 my_scipt.py
