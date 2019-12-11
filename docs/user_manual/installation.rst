======================
Installation for users
======================

With Conda
----------

The easiest way to install Capytaine is the precompiled package available on Conda_.
Download and install the `Anaconda distribution`_ or its lightweight counterpart Miniconda_.

.. _Conda: https://conda.io
.. _`Anaconda distribution`: https://www.anaconda.com/download/
.. _Miniconda: https://conda.io/miniconda.html

Capytaine requires **Python 3.6** or higher.
It has been successfully tested on Python 3.6 and 3.7, and Numpy 1.15 and 1.16.

Once Conda has been installed, run the following command in a terminal to install Capytaine::

    conda install -c conda-forge capytaine

Some optional dependencies are not automatically installed with Capytaine.
Install them with::

    conda install matplotlib
    conda install vtk

They are plotting and display libraries that are useful but not necessary for Capytaine's main purpose.
Some users have had troubles installing them, so they are not required anymore since version 1.2.
Note that some examples in the documentation might be assuming you have them installed.

With Pip
--------

The package is available on PyPI, although only as a source distribution.
That means that you'll nead a Fortran compiler [#]_ in order to install the package.
If you do, you can install Capytaine as::

    pip install numpy
    pip install capytaine

If you want, you can then install the optional dependencies::

    pip install matplotlib vtk

If you can't install a compiler, it is recommanded to use Conda instead.

.. [#] For example, on Ubuntu or Debian: :code:`sudo apt install gfortran`.

