===========================
Installation for developers
===========================

On Linux, MacOS, or Windows Subsystem for Linux (WSL)
-----------------------------------------------------

It is recommended to use a `conda environment <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_, for instance with::

    conde create --name capy_dev python=3.9 numpy=1.20 scipy pip
    conda activate capy_dev

Capytaine requires **Python 3.6** or higher.
The latest version is tested on Python 3.7 and 3.9, and Numpy 1.16 and 1.20.

You'll also need a Fortran compiler:

* **On Linux or WSL,** you can install :code:`gfortran` with the package manager of your distribution (e.g. on Debian or Ubuntu: :code:`sudo apt install gfortran`).

* **On macOS,** you can install the required compilers via `Homebrew`_. Make sure that
  the compilers installed by Homebrew are in you path (e.g., :code:`which gcc`);
  this can be accomplished by adding the relevant directories to your path::

  	export PATH="/usr/local/bin:$PATH"

  or through the use of aliases, e.g.,::

  	alias gcc=/usr/local/bin/gcc-10

.. _`Homebrew`: https://brew.sh

Then, download the source code from Github web interface or using ``git`` with::

    git clone --recurse-submodules https://github.com/capytaine/capytaine
    cd capytaine

To compile the code, install all optional dependencies, and put it in your path::

    make develop

The test suite can be run with::

    make test

If you need to recompile::

    make clean
    make develop

For instructions about how to help with the development of Capytaine, see the `contributing guide`_.

.. _`contributing guide`: https://github.com/capytaine/capytaine/blob/master/CONTRIBUTING.md
