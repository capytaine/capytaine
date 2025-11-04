===========================
Installation for developers
===========================

The development of Capytaine mostly happens on Linux or Windows Subsystem for Linux (WSL).
You can find below some instructions to compile Capytaine on macOS or Windows.
These instruction are not updated regularly and might be out-of-date.
Please open an issue on Github if you need help.

Previous versions of this documentation recommended the use of Conda_ (or Mamba_) to develop Capytaine.
This is not really necessary on Linux, where it might be simpler and faster to use a PyPI-based workflow.
Conda might nonetheless be useful on Windows in order to install development tools such as `git`, `just` or `gfortran`.

.. _Conda: https://conda.io
.. _Mamba: https://mamba.readthedocs.io/en/latest/


Installation for development using `pip` on Linux or WSL
--------------------------------------------------------

As of February 2025, the build backend used by Capytaine since version 2.0 (meson-python_) is not fully compatible with some modern tools such as `uv` for development (see https://github.com/astral-sh/uv/issues/10214).

.. _meson-python: https://mesonbuild.com/meson-python/index.html

You might need to have Python installed as well as common development tools such as ``git``, ``gfortran`` and ``just``.
On Ubuntu or Debian, this can be done with::

    sudo apt install python3 python3-pip python3-venv python-is-python3 git gfortran just

Not that ``just`` is only available in this way in recent versions of Ubuntu (24.04 and newer) and Debian (13 and newer).

Get the source code from Github using ``git``::

    git clone https://github.com/capytaine/capytaine
    cd capytaine

If you wish to contribute to the project, you might want to create a fork of the project and clone it using the SSH interface (see Github's documentation for more details).
Let us create a virtual environment in which the development version of Capytaine and its dependencies will be installed, such as::

    python -m venv /tmp/capy_dev_venv

Note that the build-system used by Capytaine does not like the use of virtual environment in the same directory as the source code (https://github.com/capytaine/capytaine/issues/396).

Activate the virtual environment with::

    source /tmp/capy_dev_venv/bin/activate  # with bash shell, change accordingly for other shell

To prepare the development environment, we'll install the required dependencies.
Assuming you have a recent enough version of ``pip``, this can be done with::

    pip install --group editable_install

and then build Capytaine in editable mode::

    pip install --no-build-isolation --editable .

The above two commands can also be executed by typing::

    just editable_install

As long as the virtual environment is activated, every import of Capytaine in Python will use the development version in this directory.
All change made to the source code should automatically affects the
installed package. (You may need to restart you Python interpreter, but
rerunning `pip install` should not be necessary.)

Call for instance the following line to check that Capytaine has been correctly installed::

    python -c 'import capytaine as cpt; print(cpt.__version__)'


Installation for development using Conda
----------------------------------------

If you are for instance on Windows, it is recommended to use Conda_ (or the alternative implementation with the same user interface Mamba_) as package manager and virtual environment manager to install Capytaine for development.
Let us first `create a new virtual environment <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_ with::

    conda create --name capy_dev python pip
    conda activate capy_dev

By default, Conda will install the latest version of Python.

If ``git`` is not available in your environment, you can install it through ``conda``::

    conda install -c conda-forge git

Then the source code can be downloaded using ``git`` with::

    git clone https://github.com/capytaine/capytaine
    cd capytaine

Alternatively, the source code can be directly downloaded from Github web interface.

Several options are available to get a Fortran compiler.
Please choose below the most relevant to your case.

.. tabs::
   .. tab:: Linux or WSL

       You can install ``gfortran`` with the package
       manager of your distribution. For instance on Debian or Ubuntu::

           sudo apt install gfortran

       Alternatively, ``gfortran`` is also available from the ``conda-forge``
       channel of the Conda package repository::

           conda install -c conda-forge gfortran


   .. tab:: macOS

       You can install ``gfortran`` via `Homebrew`_::

           brew install gcc

       Make sure that the compilers installed by Homebrew are in you path (e.g.,
       :code:`which gcc`); this can be accomplished by adding the relevant
       directories to your path::

           export PATH="/usr/local/bin:$PATH"

       or through the use of aliases, e.g.,::

           alias gcc=/usr/local/bin/gcc-10

   .. _`Homebrew`: https://brew.sh


   .. tab:: Windows (GFortran)

      The GNU toolchain, including ``gfortran`` can be installed with the help of ``conda``::

           conda install -c conda-forge m2w64-toolchain


   .. tab:: Windows (Intel)

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


Once you have a Fortran compiler installed, the same instructions as above can be used to install the Python dependencies of Capytaine::

    pip install --group editable_install

and then build Capytaine in editable mode::

    pip install --no-build-isolation --editable .

If ``just`` is not available in your environment, you can install it through ``conda``::

    conda install -c conda-forge just

and simply use the following line to install Capytaine in editable mode in your conda environment::

    just editable_install

You can check that the package is installed by running::

    python -c 'import capytaine as cpt; print(cpt.__version__)'


Building the documentation
--------------------------

In a ``pip`` or ``conda`` virtual environment (which can be the same as above or a different one), install Capytaine with the extra dependencies ``docs``::

    pip install --group docs
    pip install .[optional]

Then run the ``make`` command in the ``docs/`` directory::

    make --directory="./docs/"

and the documentation will be built in the ``docs/_build`` directory.

Alternatively, use::

    just build_docs
