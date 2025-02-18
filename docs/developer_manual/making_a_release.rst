================
Making a release
================

This page is a memo for the Capytaine developer(s) detailling the process of doing a new release.

Last commit of released version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Open a pull request similar to https://github.com/capytaine/capytaine/pull/613 to update:
- the version number in ``capytaine/__about__.py``,
- the changelog in ``docs/changelog.rst``,
- the link to previous documentations in ``docs/index.rst``.


Dry-run packaging
~~~~~~~~~~~~~~~~~

- The ``build_wheel.yaml`` Github Action has a ``workflow_dispatch`` event trigger, which means it can be run manually from Github at https://github.com/capytaine/capytaine/actions/workflows/build_wheels.yaml
Run it on version to be released.


- To dry-run the conda packaging, open a pull-request named ``release-vx.y`` on https://github.com/conda-forge/capytaine-feedstock and update the ``recipe/meta.yaml`` to use::

  source:
    git_url: https://github.com/capytaine/capytaine
    git_rev: master
    git_depth: 1

to build the version of the master branch.
The conda-forge CI will run the compilation.
DO NOT merge the PR yet.


Github release and actual packaging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Make a new release on Github https://github.com/capytaine/capytaine/releases titled "Capytaine vX.Y", creating a new tag ``vX.Y``.

- Update the PR opened in https://github.com/conda-forge/capytaine-feedstock to use the released tarball as source::

  source:
    url: https://github.com/capytaine/capytaine/archive/refs/tags/v{{ version }}.tar.gz
    sha256: {{ hash }}

If the CI is all green, the PR can be merged and the package will be uploaded to the conda-forge repository in the following hours.

- Build the source distribution, using e.g. locally ::

   uv run --with build python -m build --sdist

- The Github Action https://github.com/capytaine/capytaine/actions/workflows/build_wheels.yaml should run automatically on the tagged commit. Download the wheels on the page of the action.

- Upload the sdist and the wheels using `twine <https://twine.readthedocs.io>`_ ::

   uv run --with twine twine upload dist/*

You'll need the access token to PyPI.


Update documentation
~~~~~~~~~~~~~~~~~~~~

In the repository of the website https://github.com/capytaine/capytaine.github.io, use the latest version of the documentation available in the directory ``master`` to make the new documentation and update the symlink ``stable``::

   mv master v2.1
   rm stable
   ln -s v2.1 stable
   cp v2.0/_static/front_page_animation.webm v2.1/_static/

In the new version documentation (`v2.1` in this example), remove the banner on top of the pages ::

   sed -i '/new-doc-warning-banner/d' $(rg -l 'new-doc')

In the former stable documentation (`v2.0` in this example), add a banner on top of the pages just after the ``<body>`` tag::

   sed -i '/<body>/a   <div id="old-doc-warning-banner" style="width:100%; background-color:#F0E68C; text-align: center;">This page is part of the documentation of an old version of Capytaine. <a href="https://capytaine.org/stable/">Latest stable version is available here.</a></div>' $(rg -l "<body>")

Commit and push to https://github.com/capytaine/capytaine.github.io.


Archival on Zenodo
~~~~~~~~~~~~~~~~~~

Zenodo record at https://zenodo.org/records/14178807 should be created automatically, but you might want to review it.


Other distributions
~~~~~~~~~~~~~~~~~~~

To release a new version of https://github.com/capytaine/capytaine-standalone:
- Update pyproject.toml to change the version of Capytaine (as well as other bundled packages if necessary)
- Run `pdm lock` to update the lock file.
- Test locally with `pdm install && pdm bundle`
- Test the CI by opening a pull request or by triggering manually on Github
- Do the release. The packages are automatically linked to the release on Github.

You also might want to update https://github.com/capytaine/capytaine-extra-build-recipes


First commit of new development branch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Update version number on the development branch to ``x.y.dev``, as in https://github.com/capytaine/capytaine/pull/616.
