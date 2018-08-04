===========================
Installation for developers
===========================

It is recommended to use a `conda environment`_.
Capytaine requires **Python 3.6**.

Install the dependencies of Capytaine.

.. _`conda environment`: https://conda.io/docs/user-guide/tasks/manage-environments.html

::

    conda install numpy
    conda install -c mancellin meshmagick

Download the source code from Github web interface or using ``git`` with

::

    git clone https://github.com/mancellin/capytaine

In the main directory of Capytaine (where ``setup.py`` is located), run the following command to compile the Fortran code.

::

    python setup.py build_ext --inplace

Re-run this command later to recompile the code if you change one of the Fortran files.

Add the current directory to the Python path with the following command.

::

    pip install -e .

(This last two commands are included into the ``Makefile`` of the main directory. They can be run by ``make develop``. The command ``make clean`` delete the binaries.)
