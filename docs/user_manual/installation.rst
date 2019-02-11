======================
Installation for users
======================


With conda
----------

On Windows
~~~~~~~~~~

The binary Windows version is currently broken, sorry.
Please take a look at this `Github issue <https://github.com/mancellin/capytaine/issues/1>`_ if you can help.

On Linux
~~~~~~~~

The easiest way to install the package is through Conda_.
Download and install the `Anaconda distribution`_ or its lightweight counterpart Miniconda_.
Capytaine is developped and tested with **Python 3.6**.
It should also work fine on newer versions of Python.

.. _Conda: https://conda.io
.. _`Anaconda distribution`: https://www.anaconda.com/download/
.. _Miniconda: https://conda.io/miniconda.html

Once Conda has been installed, run the following command in a terminal to install Capytaine::

    conda install -c mancellin capytaine

All the necessary code from Nemoh and Meshmagick is already included into Capytaine and all the other required dependencies should be automatically installed.

In case you get an error message because of a missing module, you can try to run::

    conda install attrs numpy scipy pandas xarray matplotlib vtk libgfortran=3

to install them all manually.

On Mac
~~~~~~

If you are Mac-savvy and would like to help me building and testing Mac binaries, please let me know!


With pip
--------

Binaries are not available on PyPI (:code:`pip install`) at the moment.
If you need them, please let me know.
