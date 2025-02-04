======================
Installation for users
======================

Capytaine requires Python 3.8 or higher.
It is regularly tested on Windows, MacOS and Linux using `all currently supported version of Python <https://devguide.python.org/versions/>`_

Precompiled packages are distributed on `PyPI <https://pypi.org/project/capytaine/>`_ and `Conda-forge <https://conda-forge.org/>`_ for Windows, MacOS (both Intel and ARM processors) and Linux for all supported versions of Python.
However, it might take a few weeks before precompiled packages for the latest version of Python (3.13 at the time of writing) are available and you might prefer using an earlier version of Python.


On a cloud platform
-------------------

For a quick try of Capytaine without installing anything on your computer, you can use an online Python-based computing environment such as `CoCalc <https://cocalc.com/>`_, on which Capytaine is already installed by default, or `Google Colab <https://colab.research.google.com/>`_.
On such a `Jupyter <https://jupyter.org/>`_-based environment, Capytaine can be installed by running the following command in a cell::

    %pip install capytaine

Then run the following line to check which version of Capytaine has been installed::

    import capytaine as cpt; print(cpt.__version__)

You may need to restart the computing environment (kernel) of the notebook for the installation to be effective.

All the core feature of Capytaine are accessible from such a Jupyter-based environment, except for some 3D visualization tools.


Locally with the `uv` package manager
-------------------------------------

As of 2025, the best compromise between ease-of-use, speed and flexibility to install Capytaine on your machine is to use the `uv <https://docs.astral.sh/uv/>`_ package manager.
Once you have installed `uv`, run the following command to run a file script::

    uv run --with capytaine --script path/to/my_script.py

`uv` will take care of installing Python and all required dependencies on the fly.
You can start an interactive console with Capytaine available as follows::

    uv run --with capytaine --with ipython ipython

Or a Matlab-like development environment with::

    uv run --with capytaine --with spyder spyder

Or the Jupyter notebook interface::

    uv run --with capytaine --with jupyter jupyter lab

Execute the following Python code to check that Capytaine is correctly installed::

    import capytaine as cpt; print(cpt.__version__)

More optional dependencies can be specified, as well as specific versions of capytaine::

    uv run --with "capytaine==2.1" --with matplotlib --with meshio --script path/to/my_script.py

The following optional dependencies can be used together with Capytaine.

+------------+---------------------------------------------------------------------------------+
| Name       | Usage                                                                           |
+============+=================================================================================+
| matplotlib | To plot graphs. Used in several examples in the documentation and the cookbook. |
+------------+---------------------------------------------------------------------------------+
| vtk        | For 3D visualization                                                            |
+------------+---------------------------------------------------------------------------------+
| meshio     | To load more mesh formats                                                       |
+------------+---------------------------------------------------------------------------------+
| netcdf4    | To export in NetCDF4 format                                                     |
+------------+---------------------------------------------------------------------------------+
| joblib     | For parallel resolution                                                         |
+------------+---------------------------------------------------------------------------------+

You can ask UV to install them all at the same time as Capytaine with the following command::

    uv run --with "capytaine[optional]" --script path/to/my_script.py


As a standalone executable
--------------------------

An experimental distribution of Capytaine bundled with a full Python distribution in a single executable file can be found at `<https://github.com/capytaine/capytaine-standalone>`_.
Please refer to the instruction on that page for download and usage.

The standalone executable is the simplest way to use Capytaine locally, although it has some limitations, such a longer startup time and the current lack of interactive Matplotlib figures.

You can check the bundled version of Capytaine with the following command::

    .\ipython-with-capytaine-windows.exe -c 'print(cpt.__version__)'

(or the corresponding file name on other platforms than Windows).

Installing with pip package manager
-----------------------------------

Since version 2.0, Capytaine is available as precompiled package on all platforms on `PyPI <https://pypi.org/project/capytaine/>`_, the package registry used by the ``pip`` command. After installing a Python interpreter, run the following command line in a terminal to install Capytaine and its dependencies::

    python -m pip install capytaine

Then run the following line to check that the latest version of Capytaine has been installed::

    python -c 'import capytaine as cpt; print(cpt.__version__)'

You might want to use a `virtual environment <https://docs.python.org/3/library/venv.html>`_ to install Capytaine independently of your other Python packages and avoid any risk of dependency conflict.
If you are using an IDE, you can install Capytaine in a virtual environment using the graphical interface such as `in PyCharm <https://www.jetbrains.com/help/pycharm/creating-virtual-environment.html>`_ or `in VSCode <https://code.visualstudio.com/docs/python/environments#_creating-environments>`_.

The package can also be installed by other modern PyPI-based Python package managers, such as UV_ (see above), PDM_ or poetry_.

.. _UV: https://docs.astral.sh/uv/
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

It is often more efficient to specify all the optional packages you'd like in your environment from the start when creating it, such as in the following example::

    conda create --name capy_and_other_env --channel conda-forge capytaine jupyter matplotlib vtk


More build recipes
------------------

More advanced build recipes for Capytaine are available in the dedicated repository `https://github.com/capytaine/capytaine-extra-build-recipes`_.
In particular, build recipes for Docker and Guix might be useful for reproducible computations.
