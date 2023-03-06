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
Download and install the `Anaconda distribution`_ or one of its lightweight counterparts Miniconda_ and Miniforge_.

.. _Conda: https://conda.io
.. _`Anaconda distribution`: https://www.anaconda.com/download/
.. _Miniconda: https://conda.io/miniconda.html
.. _Miniforge: https://github.com/conda-forge/miniforge

Once Conda has been installed, you can install Capytaine from the `conda-forge` channel.
It is recommended to do the installation into a `dedicated virtual environment <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`_ (here arbitrarily named ``capytaine_env``)::

    conda create --name capytaine_env --channel conda-forge capytaine

Then activate the environment to use it on the command line with::

    conda activate capytaine_env
    
or set it in the project configuration of your IDE (for instance see `the documentation of PyCharm <https://www.jetbrains.com/help/pycharm/conda-support-creating-conda-virtual-environment.html>`_ or the `documentation of Spyder <https://github.com/spyder-ide/spyder/wiki/Working-with-packages-and-environments-in-Spyder#working-with-other-environments-and-python-installations>`_).

Alternatively, Capytaine can be installed in an existing environment with the following command::

    conda install --channel conda-forge capytaine

You can check which version of Capytaine has been installed by opening a Python shell and running::

    import capytaine; print(capytaine.__version__)

Optional dependencies
---------------------

All the required dependencies should be installed automatically with the Conda package.
More optional dependencies can be manually installed.
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

The following command will create a Docker image based on Ubuntu 22.04 with the version v1.5 of Capytaine::

    docker build -t capytaine:v1.5 https://github.com/capytaine/capytaine.git#v1.5

Replace :code:`v1.5` by :code:`master` to download instead the latest development version.
Use the following command to open an IPython shell in which Capytaine can be imported::

    docker run -it capytaine:v1.5 ipython3

Or the following command to make the current directory accessible from the Docker image and run the file :code:`my_script.py` from the current directory::

    docker run -it -v $(pwd):/home/user capytaine:v1.5 python3 my_scipt.py

Note that graphical displays (matplotlib, vtk, ...) might require a complex setup to work from the Docker image.
