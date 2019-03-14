===========================
Installation for developers
===========================

Capytaine requires **Python 3.6** or higher.
It has been successfully tested on Python 3.6 and 3.7, and Numpy 1.15 and 1.16.

It is recommended to use a `conda environment`_.

.. _`conda environment`: https://conda.io/docs/user-guide/tasks/manage-environments.html

Ensure that Numpy is installed in your enviroment::

    conda install numpy

You'll also need a Fortran compiler:

* On Linux, you can install `gfortran` with the package manager of your distribution (e.g. `sudo apt install gfortran`).

* On Windows, the code can be compiled with MinGW.
  Add the directory with the `gfortran` binary to your path. For instance with `set PATH=C:\\mingw-w64\\x86_64-7.2.0-posix-seh-rt_v5-rev1\\mingw64\\bin;%PATH%`.
  You should also let Python know about the compiler by creating a file with the following two lines::

    echo [build]
    echo compiler=mingw32

  as `C:\\path\\to\\anaconda\\Lib\\distutils\\distutils.cfg`.

Then, download the source code from Github web interface or using ``git`` with::

    git clone https://github.com/mancellin/capytaine

In the main directory of Capytaine (where ``setup.py`` is located), run the following command to compile the Fortran code::

    python setup.py build_ext --inplace

Re-run this command later to recompile the code if you change one of the Fortran files.

Add the current directory to the Python path with the following command::

    pip install -e .

(This last two commands are included into the ``Makefile`` of the main directory.
They can be run by ``make develop``.
The command ``make clean`` deletes the binaries.)
