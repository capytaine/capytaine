#!/usr/bin/env python
# coding: utf-8

import os

from setuptools import dist
dist.Distribution().fetch_build_eggs(['numpy'])

from numpy.distutils.core import Extension, setup

########################
#  Fortran extensions  #
########################

Delhommeau_source = [
        "capytaine/green_functions/libDelhommeau/src/constants.f90",
        "capytaine/green_functions/libDelhommeau/src/Delhommeau_integrals.f90",
        "capytaine/green_functions/libDelhommeau/src/old_Prony_decomposition.f90",
        "capytaine/green_functions/libDelhommeau/src/Green_Rankine.f90",
        "capytaine/green_functions/libDelhommeau/src/Green_wave.f90",
        "capytaine/green_functions/libDelhommeau/src/matrices.f90",
    ]

Delhommeau_extension = Extension(
    name="capytaine.green_functions.Delhommeau_f90",
    sources=Delhommeau_source,
    extra_f90_compile_args=['-O2', '-fopenmp', '-cpp'],
    extra_link_args=['-fopenmp'],
    # # Uncomment the following lines to get more verbose output from f2py.
    # define_macros=[
    #     ('F2PY_REPORT_ATEXIT', 1),
    #     ('F2PY_REPORT_ON_ARRAY_COPY', 1),
    # ],
)

XieDelhommeau_extension = Extension(
    name="capytaine.green_functions.XieDelhommeau_f90",
    sources=Delhommeau_source,
    extra_f90_compile_args=['-O2', '-fopenmp', '-cpp', '-DXIE_CORRECTION'],
    extra_link_args=['-fopenmp'],
    # # Uncomment the following lines to get more verbose output from f2py.
    # define_macros=[
    #     ('F2PY_REPORT_ATEXIT', 1),
    #     ('F2PY_REPORT_ON_ARRAY_COPY', 1),
    # ],
)


########################################################
#  Read version number and other info in __about__.py  #
########################################################

base_dir = os.path.dirname(__file__)
src_dir = os.path.join(base_dir, "capytaine")

about = {}
with open(os.path.join(src_dir, "__about__.py")) as f:
    exec(f.read(), about)


##########
#  Main  #
##########

if __name__ == "__main__":
    setup(name=about["__title__"],
          version=about["__version__"],
          description=about["__description__"],
          author=about["__author__"],
          license=about["__license__"],
          url=about["__uri__"],
          packages=[
              'capytaine',
              'capytaine.meshes',
              'capytaine.matrices',
              'capytaine.bodies',
              'capytaine.bodies.predefined',
              'capytaine.bem',
              'capytaine.green_functions',
              'capytaine.post_pro',
              'capytaine.ui',
              'capytaine.ui.vtk',
              'capytaine.io',
              'capytaine.tools',
          ],
          install_requires=[
              'numpy',
              'scipy',
              'pandas>=1.3',
              'xarray',
          ],
          extras_require={
            'develop': [
              'pytest',
              'hypothesis',
              'ipython',
              'matplotlib',
              'vtk',
              'meshio',
              'pygmsh',
              'gmsh',
              'quadpy',
              'bemio @ git+https://github.com/michaelcdevin/bemio.git@master-python3#egg=bemio',
              'sphinx',
              'sphinxcontrib-proof',
            ],
            'ci': [
              'pytest',
              'hypothesis',
                ],
            'extra': [
              'ipython',
              'matplotlib',
              'vtk',
              'meshio',
              'pygmsh',
              'gmsh',
              'quadpy',
              'bemio @ git+https://github.com/michaelcdevin/bemio.git@master-python3#egg=bemio',
            ]
          },
          entry_points={
              'console_scripts': [
                  'capytaine=capytaine.ui.cli:main',
              ],
          },
          ext_modules=[
              Delhommeau_extension,
              XieDelhommeau_extension,
          ],
          )
