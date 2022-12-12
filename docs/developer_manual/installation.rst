===========================
Installation for developers
===========================

On Linux, MacOS, or Windows Subsystem for Linux (WSL)
-----------------------------------------------------

It is recommended to use a `conda environment <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_, for instance with::

    conda create --name capy_dev python
    conda activate capy_dev

By default, conda will install the latest version of Python.
Capytaine requires Python 3.6 and is compatible with `all currently supported version of Python <https://devguide.python.org/versions/>`_.


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


On Windows
----------

Microsoft Visual Studio is required for linking the Fortran binaries

    * https://visualstudio.microsoft.com/downloads/
    * During installation check the box to include :code:`Desktop development with C++`

Intel oneAPI HPC toolkit is required for compiling the Fortran binaries (you do not need the base kit)

    * https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html
    * Install to the default file location
    
Create a **"LIB"** environment variable to point towards the intel directory for compiler :code:`.lib` files

    * If oneAPI is installed to the default location, assign the LIB user variable a value of:
    
        :code:`C:\Program Files (x86)\Intel\oneAPI\compiler\2022.1.0\windows\compiler\lib\intel64_win`
    
    * If oneAPI is installed to a different location then adjust the path above as necessary

Test if your Fortran compiler was installed correctly by entering :code:`ifort` on your command line

Open the anaconda powershell and create a new Python environment (by default, with the latest version of Python) for Capytaine-related development (e.g. :code:`capy_dev`)::
    
    conda create --name capy_dev python
    conda activate capy_dev
        
Clone the Capytaine repo to your preferred location (e.g. "C:/code/")::
        
    cd C:/code/
    git clone --recurse-submodules https://github.com/capytaine/capytaine.git
        
Install Capytaine as a developer!::
    
    cd capytaine
    pip install -e .

Be sure to check setup.py => install_requires = [...] to ensure that your environment has all required packages installed. You can check your environment's packages using:::

    conda list
        
If any packages are missing simply install them using:::
    
    pip install <package name>


Contributing
------------

For instructions about how to help with the development of Capytaine, see the `contributing guide`_.

.. _`contributing guide`: https://github.com/capytaine/capytaine/blob/master/CONTRIBUTING.md
