Installation for users
======================

On Windows or Linux
-------------------

The easiest way to install the package is through Conda_.
Download and install the `Anaconda distribution`_ or its lightweight counterpart Miniconda_.
Capytaine requires **Python 3.6**.

.. _Conda: https://conda.io
.. _`Anaconda distribution`: https://www.anaconda.com/download/
.. _Miniconda: https://conda.io/miniconda.html

Once Conda has been installed, run the following command in a terminal to install a dependency:

::

    conda install -c mancellin meshmagick

This is a fork of Meshmagick_ updated for Python 3 and Capytaine.

.. _Meshmagick: https://github.com/LHEEA/meshmagick

Then the code itself can be installed by the following command:

::

    conda install -c mancellin capytaine

All the required dependencies should be automatically installed with the previous command.
In case you get an error message because of a missing module, you can try running

::

    conda install attrs numpy scipy pandas xarray matplotlib vtk

to install them all manually.

All the necessary code from Nemoh is already included into Capytaine.

On Mac
------

If you are Mac-savvy and would like to help me building and testing Mac binaries, please let me know!

