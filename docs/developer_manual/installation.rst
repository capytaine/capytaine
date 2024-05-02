===========================
Installation for developers
===========================

Development environment
-----------------------

To work on the source of Capytaine, it is recommended to use Conda_ (or the alternative implementation with the same user interface Mamba_) as package manager and virtual environment manager.
Other Python package manager such as PDM_ can also be used.
However Conda/Mamba can simplify the installation of some non-Python tools such as the Fortran compiler and will thus be used as example for the rest of this page.

Please check also the section on Conda in the :doc:`user installation instructions </user_manual/installation>`.

.. _Conda: https://conda.io
.. _Mamba: https://mamba.readthedocs.io/en/latest/
.. _PDM: https://pdm.fming.dev/latest/

Let us first `create a new virtual environment <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_ with::

    conda create --name capy_dev python pip
    conda activate capy_dev

By default, Conda will install the latest version of Python.
Capytaine requires Python 3.7 or higher, and is compatible with `all currently supported version of Python <https://devguide.python.org/versions/>`_.

Getting the source code
-----------------------

If ``git`` is not available in your environment, you can install it through ``conda``::

    conda install -c conda-forge git

Then the source code can be downloaded using ``git`` with::

    git clone https://github.com/capytaine/capytaine
    cd capytaine

Alternatively, the source code can be directly downloaded from Github web interface.

Getting a Fortran compiler
--------------------------

Several options are available to get a Fortran compiler.
Please choose below the most relevant to your case.

.. collapse:: GFortran compiler on Linux or Windows Subsystem for Linux

    You can install ``gfortran`` with the package
    manager of your distribution. For instance on Debian or Ubuntu::

        sudo apt install gfortran

    Alternatively, ``gfortran`` is also available from the ``conda-forge``
    channel of the Conda package repository::

        conda install -c conda-forge gfortran


.. collapse:: GFortran compiler on macOS

    You can install ``gfortran`` via `Homebrew`_::

        brew install gcc

    Make sure that the compilers installed by Homebrew are in you path (e.g.,
    :code:`which gcc`); this can be accomplished by adding the relevant
    directories to your path::

        export PATH="/usr/local/bin:$PATH"

    or through the use of aliases, e.g.,::

        alias gcc=/usr/local/bin/gcc-10

.. _`Homebrew`: https://brew.sh


.. collapse:: GFortran on Windows

   The GNU toolchain, including ``gfortran`` can be installed with the help of ``conda``::

        conda install -c conda-forge m2w64-toolchain


.. collapse:: Intel compiler on Windows

    Microsoft Visual Studio is required for linking the Fortran binaries

        * https://visualstudio.microsoft.com/downloads/
        * During installation check the box to include :code:`Desktop development with C++`

    Intel oneAPI HPC toolkit is required for compiling the Fortran binaries (you do not need the base kit)

        * https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html
        * Install to the default file location

    Create a ``LIB`` environment variable to point towards the intel directory for compiler ``.lib`` files

        * If oneAPI is installed to the default location, assign the LIB user variable a value of::

            C:\Program Files (x86)\Intel\oneAPI\compiler\2022.1.0\windows\compiler\lib\intel64_win

        * If oneAPI is installed to a different location then adjust the path above as necessary

    Test if your Fortran compiler was installed correctly by entering :code:`ifort` on your command line


Compiling and installing the code
---------------------------------

The ``Makefile`` file in Capytaine repository contains short forms for the most common commands required to build Capytaine.

If ``make`` is not available in your environment, you can install it through ``conda``::

    conda install -c conda-forge make

To compile the code and install it in the current environment, you can run::

    make install

which is just synonym of::

    pip install .

You can check that the package is installed by running::

    python -c 'import capytaine as cpt; print(cpt.__version__)'

or by checking the complete list of packages installed in the current environment with::

    conda list

.. note::

    If you have an error of the form::

        ModuleNotFoundError:: No module named 'capytaine.green_functions.libs.Delhommeau_float64'

    when importing Capytaine, it may be because the Python interpreter is
    trying to load the content of the local directory ``capytaine`` (containing
    only the source code) and not the actual compiled package.

    Running ``python`` from any other directory on your system should fix the
    issue, since there won't be a local ``capytaine`` directory to confuse the
    module importer.

    Alternatively, recent versions of Python (>=3.11) have the ``-P`` option
    which will disable the loading of the local ``capytaine`` directory.

When using ``make install``, you will need to re-run the installation
for any change made to the code to take effect in the installed version. For
development, it is more convenient to use instead::

    make develop

Then all change made to the source code should automatically affects the
installed package. (You may need to restard you Python interpreter.)

Testing
-------

To check that the installed packaged is working fine, you can run the test suite with Pytest.
If Pytest is not already install, it can be done with::

    pip install pytest

Then run the following command from the root of Capytaine repository to test the code::

    python -m pytest

Alternatively, `Nox <https://nox.thea.codes>`_ can be used to set up an isolated environment and build and test the code.
After installing Nox, use::

    nox -s build_and_test_on_locked_env

to test the current source code in an environment with fixed versions of Capytaine's dependencies, whereas::

    nox -s build_and_test_on_latest_env

will test the code in an environment using the latest available version of Capytaine's dependencies.


Building the documentation
--------------------------

TODO

Contributing
------------

For instructions about how to help with the development of Capytaine, see the `contributing guide`_.

.. _`contributing guide`: https://github.com/capytaine/capytaine/blob/master/CONTRIBUTING.md
