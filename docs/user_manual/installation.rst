======================
Installation for users
======================

Capytaine is available on Windows, MacOS [#]_ and Linux.

.. [#] For the latest informations on the Apple arm64 architectures, see https://github.com/capytaine/capytaine/issues/190

Capytaine requires Python 3.7 or higher.
It is compatible with `all currently supported version of Python <https://devguide.python.org/versions/>`_.


On a cloud platform
-------------------

For a quick try of Capytaine without installing anything on your computer, you can use an online Python-based computing environment such as `CoCalc <https://cocalc.com/>`_, on which Capytaine is already installed by default, or `Google Colab <https://colab.research.google.com/>`_.
On such a `Jupyter <https://jupyter.org/>`_-based environment, Capytaine can be installed by running the following command in a cell::

    %pip install capytaine

Then run the following line to check that the latest version of Capytaine has been installed::

    import capytaine as cpt; print(cpt.__version__)

You may need to restart the computing environment (kernel) of the notebook for the installation to be effective.

All the core feature of Capytaine are accessible from such a Jupyter-based environment, except for some 3D visualization tools.


As a standalone executable
--------------------------

This is a work in progress. See :pull:`324` for progress updates.


Installing with pip package manager
-----------------------------------

Since version 2.0, Capytaine is available as precompiled package on all platform on `PyPI <https://pypi.org/project/capytaine/>`_, the package registry used by the ``pip`` command. After installing a Python interpreter, run the following command line in a terminal to install Capytaine and its dependencies::

    python -m pip install capytaine

Then run the following line to check that the latest version of Capytaine has been installed::

    python -c 'import capytaine as cpt; print(cpt.__version__)'

You might want to use a `virtual environment <https://docs.python.org/3/library/venv.html>`_ to install Capytaine independently of your other Python packages and avoid any risk of dependency conflict.

The package can also be installed by other modern PyPI-based Python package managers, such as PDM_ or poetry_.

.. _PDM: https://pdm.fming.dev
.. _poetry: https://python-poetry.org

Installing with Conda package manager
-------------------------------------

Capytaine is also available in the Anaconda package repository, that can be accessed with the `Anaconda distribution`_ or one of its lightweight counterparts Miniconda_ and Miniforge_.

.. _Conda: https://conda.io
.. _`Anaconda distribution`: https://www.anaconda.com/download/
.. _Miniconda: https://conda.io/miniconda.html
.. _Miniforge: https://github.com/conda-forge/miniforge
.. _Mamba: https://mamba.readthedocs.io/en/latest/

.. note::
    If you experience very long processing time when installing a package with ``conda``, you might want to `install the libmamba solver with ``conda`` <https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community>`_ or `fully replace ``conda`` with Mamba_.

Once Conda has been installed, you can install Capytaine from the `conda-forge` channel.
It is recommended to do the installation into a `dedicated virtual environment <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments>`_ (here arbitrarily named ``capytaine_env``)::

    conda create --name capytaine_env --channel conda-forge capytaine

Then activate the environment to use it on the command line with::

    conda activate capytaine_env

or set it in the project configuration of your IDE (for instance see `the documentation of PyCharm <https://www.jetbrains.com/help/pycharm/conda-support-creating-conda-virtual-environment.html>`_, `the documentation of VSCode <https://code.visualstudio.com/docs/python/environments#_working-with-python-interpreters>`_ or the `documentation of Spyder <https://github.com/spyder-ide/spyder/wiki/Working-with-packages-and-environments-in-Spyder#working-with-other-environments-and-python-installations>`_).

Alternatively, Capytaine can be installed in an existing environment with the following command::

    conda install --channel conda-forge capytaine

You can check which version of Capytaine has been installed by running the following command line::

    python -c 'import capytaine as cpt; print(cpt.__version__)'

The latest version is currently |version|.


Optional dependencies
---------------------

All the required dependencies should be installed automatically when installing with ``pip`` or ``conda``.
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

After creating the Conda environment containing Capytaine, you can add more packages to this environment by activating it with ``conda activate`` and then using the ``conda install`` or ``pip install`` commands.
However, it is often more efficient to specify the packages you'd like in your environment from the start when creating it, such as in the following example::

    conda create --name capy_and_other_env --channel conda-forge capytaine jupyter matplotlib vtk


With Docker
-----------

The following command will create a Docker image based on Ubuntu 22.04 with the version v2.0 of Capytaine::

    docker build -t capytaine:v2.0 https://github.com/capytaine/capytaine.git#v2.0

Replace :code:`v2.0` by :code:`master` to download instead the latest development version.
Use the following command to open an IPython shell in which Capytaine can be imported::

    docker run -it capytaine:v2.0 ipython3

Or the following command to make the current directory accessible from the Docker image and run the file :code:`my_script.py` from the current directory::

    docker run -it -v $(pwd):/home/user capytaine:v2.0 python3 my_scipt.py

Note that graphical displays (matplotlib, vtk, ...) might require a complex setup to work from the Docker image.
