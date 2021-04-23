===========================
Installation for developers
===========================

Capytaine requires **Python 3.6** or higher.
It has been successfully tested on Python 3.6, 3.7, and 3.9, and Numpy 1.15, 1.16, and 1.20.

It is recommended to use a `conda environment`_.

.. _`conda environment`: https://conda.io/docs/user-guide/tasks/manage-environments.html

Ensure that Numpy is installed in your enviroment::

    conda install numpy

You'll also need a Fortran compiler:

* **On Linux,** you can install :code:`gfortran` with the package manager of your distribution (e.g. :code:`sudo apt install gfortran`).

* **On Windows,** the code can be compiled with MinGW.
  Add the directory with the `gfortran` binary to your path. For instance with :code:`set PATH=C:\\mingw-w64\\x86_64-7.2.0-posix-seh-rt_v5-rev1\\mingw64\\bin;%PATH%`.
  You should also let Python know about the compiler by creating a file with the following two lines::

    echo [build]
    echo compiler=mingw32

  as `C:\\path\\to\\anaconda\\Lib\\distutils\\distutils.cfg`.

* **On macOS,** you can install the required compilers via `Homebrew`_. Make sure that
  the compilers installed by Homebrew are in you path (e.g., :code:`which gcc`); 
  this can be accomplished by adding the relevant directories to your path::

  	export PATH="/usr/local/bin:$PATH"

  or through the use of aliases, e.g.,::
  
  	alias gcc=/usr/local/bin/gcc-10
  
.. _`Homebrew`: https://brew.sh

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

You will also likely want to install Capytaine's optional dependencies::

	conda install matplotlib vtk pytest quadpy
	pip install pygmsh 
