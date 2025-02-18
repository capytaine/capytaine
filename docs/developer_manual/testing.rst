
=================
Testing Capytaine
=================

Running the test suite with Pytest
----------------------------------

To check that the installed packaged is working fine, you can run the test suite with `Pytest`_.
If Pytest is not already installed, it can be done with e.g.::

    pip install pytest

Then run the following command from the root of Capytaine repository to test the code::

    python -m pytest

Note that some tests are skipped if the required optional dependencies are not
installed, hence it is advisable to install all of Capytaine optional
dependencies (see :doc:`../user_manual/installation.rst`) to fully test
the code.

Some useful features of Pytest include:

- stopping at the first failure::

   python -m pytest --maxfail=1

- running all test starting with the ones that failed in the previous run::

   python -m pytest --failed-first


.. _`Pytest`: https://docs.pytest.org/


Testing in isolated environments with Nox
-----------------------------------------

`Nox`_ can be used to set up an isolated environment
and build and run the unit test suite as well as some of the example code from
the cookbook (``docs/user_manual/examples/`` directory).
The test protocols are defined in the ``noxfile.py`` file at the root of the
repository.

.. _`Nox`: https://nox.thea.codes

The ``build_and_test_on_locked_env`` Nox session is used to test the current
source code in an environment with fixed versions of Capytaine's dependencies.
This is meant to test changes in Capytaine without interferences from possible
changes in dependencies.
Two environments are predefined for this test, one is older and meant to be
used with the oldest version of Python supported by Capytaine (Python 3.8 at
the time of writing), while the other is more recent and is meant to be used
with a recent version of Python (Python 3.12 at the time of writing).
Their lockfiles can be found in the ``pytest/envs/`` directory.
Assuming you have `uv <https://docs.astral.sh/uv/>`_ installed, you can them with::

    uv run --with nox --python 3.8 nox -s build_and_test_on_locked_env
    uv run --with nox --python 3.12 nox -s build_and_test_on_locked_env

This Nox session is the main part of the
``.github/workflows/test_new_commits.yaml`` Github Actions workflow that is run
at each new commit and pull request.

Alternatively, the Nox sessions ``build_and_test_on_latest_env`` and
``editable_build_and_test_on_latest_env`` uses the latest dependencies
available on PyPI to test Capytaine, respectively by installing it as a normal
package or by installing it in development mode as described above on this
page.
You can run it locally with, e.g.::

    nox -s build_and_test_on_latest_env

assuming you have ``nox`` installed.
If you have ``uv`` installed, you can prefix the above command with ``uv run
--with nox`` to install ``nox`` on-the-fly and run the tests.

It is the main part of the
``.github/workflows/test_with_latest_dependencies.yaml`` Github Actions workflow
that is run twice a month on the main branch of the Github repository.

Finally the session ``build_and_test_on_nightly_builds`` fetches yet-unreleased
versions of Capytaine dependencies and run the same tests. It is mostly meant
to anticipate breaking changes in the dependencies, such as the Numpy 1 to
Numpy 2 transition.


Testing the Fortran core
------------------------

The core Fortran routine of Capytaine located at
``capytaine/green_functions/libDelhommeau/`` includes some example meant to
demonstrate its usage on its own outside of Capytaine.
See ``capytaine/green_functions/libDelhommeau/Makefile`` for details.
The good compilation of this standalone Fortran code should be tested regularly.
It is done in CI for each new commit in ``.github/workflows/test_new_commits.yaml``.


While testing new developments in the Fortran core it is also useful to test
the compilation of Capytaine's Fortran core AND its Python binding, but without
the full Capytaine test suite.
This can be done with Capytaine's main ``Makefile`` as::

   make test_fortran_compilation


Sanity checks with Pre-commit
-----------------------------

Capytaine repository includes a config file for `pre-commit`_.
It is meant to catch common mistakes before creating a new commit.
While not as important as the functional tests described above, it is
recommended to install pre-commit and set it up in Capytaine repository while
developping Capytaine.

.. _`pre-commit`: https://pre-commit.com/


Testing in CI with Github Actions
---------------------------------

The directory ``.github/workflows/`` contains the definition of the Github
Actions that are run to test Capytaine in CI.
The workflows ``test_new_commits.yaml`` and
``test_with_latest_dependencies.yaml`` are mostly thin wrapper around the Nox
sessions described above.
The former is run on each new commit, while the latter is run periodically.
Follow the instruction above in the section about Nox to run them locally.
Both only run on Linux, testing other platforms is only done with the
``build_wheel.yaml`` workflow described below.

The workflow ``build_docs.yaml`` builds the documentation and push it to
`https://capytaine.org/master/`, whenever a commit modifies the sources in the
``docs/`` directory.

The workflow ``build_wheels.yaml`` builds the wheels (precompiled packages) for
all supported platform that are distributed on PyPI.
It is run when a new version is released.
It is also run periodically to check that the current main branch of Capytaine
can be built and run without issue on all platforms.
