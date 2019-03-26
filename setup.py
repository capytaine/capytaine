#!/usr/bin/env python
# coding: utf-8

import os
import re

from numpy.distutils.core import Extension, setup

#######################
#  Fortran extension  #
#######################

NemohCore = Extension(
    name="capytaine.bem.NemohCore",
    sources=[
        "capytaine/bem/NemohCore/constants.f90",
        "capytaine/bem/NemohCore/Green_Rankine.f90",
        "capytaine/bem/NemohCore/Initialize_Green_wave.f90",
        "capytaine/bem/NemohCore/Green_wave.f90",
        "capytaine/bem/NemohCore/old_Prony_decomposition.f90",
    ],
    extra_compile_args=['-O2', '-fopenmp'],
    extra_f90_compile_args=['-O2', '-fopenmp'],
    extra_link_args=['-fopenmp'],
    # # Uncomment the following lines to get more verbose output from f2py.
    # define_macros=[
    #     ('F2PY_REPORT_ATEXIT', 1),
    #     ('F2PY_REPORT_ON_ARRAY_COPY', 1),
    # ],
)

#########################
#  Read version number  #
#########################

here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    with open(os.path.join(here, *parts), 'r') as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


if __name__ == "__main__":

    ##########
    #  Main  #
    ##########

    setup(name='capytaine',
          version=find_version('capytaine', '__init__.py'),
          description='a Python-based linear potential flow solver',
          url='http://github.com/mancellin/capytaine',
          author='Matthieu Ancellin',
          author_email='matthieu.ancellin@ucd.ie',
          license='GPL-3.0',
          packages=[
              'capytaine',
              'capytaine.meshes',
              'capytaine.matrices',
              'capytaine.bodies',
              'capytaine.bodies.predefined',
              'capytaine.bem',
              'capytaine.post_pro',
              'capytaine.ui',
              'capytaine.ui.vtk',
              'capytaine.io',
          ],
          install_requires=[
              'attrs',
              'numpy',
              'scipy',
              'pandas',
              'xarray',
              'matplotlib',
              'vtk',
          ],
          entry_points={
              'console_scripts': [
                  'capytaine=capytaine.ui.cli:main',
              ],
          },
          ext_modules=[NemohCore],
          )
