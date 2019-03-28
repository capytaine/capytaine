======================
Installation for users
======================

.. warning::
    The binary package for Windows is currently broken, sorry.
    (Please take a look at this `Github issue <https://github.com/mancellin/capytaine/issues/1>`_ if you can help.)

    However, the Linux version works fine on the *Windows Subsystem for Linux*.
    Until the native Windows version is fixed, it is the recommended way to install Capytaine on Windows.

With conda on Linux or Mac
--------------------------

The easiest way to install Capytaine is the precompiled package available on Conda_.
Download and install the `Anaconda distribution`_ or its lightweight counterpart Miniconda_.

.. _Conda: https://conda.io
.. _`Anaconda distribution`: https://www.anaconda.com/download/
.. _Miniconda: https://conda.io/miniconda.html

Capytaine requires **Python 3.6** or higher.
It has been successfully tested on Python 3.6 and 3.7, and Numpy 1.15 and 1.16.

Once Conda has been installed, run the following command in a terminal to install Capytaine::

    conda install -c conda-forge capytaine

All the necessary code from Nemoh and Meshmagick is already included into Capytaine and all the other required dependencies should be automatically installed.

In case you get an error message because of a missing module, you can try to install them all manually::

    conda install "numpy>=1.15,<1.17" libgfortran=3 attrs scipy matplotlib pandas xarray vtk 


With Pip
--------

The package is available on PyPI, although only as a source distribution.
That means that you'll nead a Fortran compiler [#]_ in order to install the package.
If you do, you can install Capytaine as::

    pip install numpy
    pip install capytaine

If you can't install a compiler, it is recommanded to use Conda instead.

.. [#] For example, on Ubuntu or Debian: :code:`sudo apt install gfortran`.

