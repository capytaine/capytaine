#!/usr/bin/env python
# coding: utf-8

from numpy.distutils.core import Extension, setup

VERSION = '0.4a'

Green = Extension(
    name="capytaine._Green",
    sources=[
        "capytaine/NemohCore/constants.f90",
        "capytaine/NemohCore/Green_1.f90",
        "capytaine/NemohCore/Initialize_Green_2.f90",
        "capytaine/NemohCore/Green_2.f90",
        "capytaine/NemohCore/old_Prony_decomposition.f90",
    ],
    # # Uncomment the following lines to get more verbose output from f2py.
    # define_macros=[
    #     ('F2PY_REPORT_ATEXIT', 1),
    #     ('F2PY_REPORT_ON_ARRAY_COPY', 1),
    # ],
)

Wavenumber = Extension(
    name="capytaine._Wavenumber",
    sources=[
        "capytaine/NemohCore/Wavenumber.f90"
    ]
)

if __name__ == "__main__":
    setup(name='capytaine',
          version=VERSION,
          description='Nemoh python wrapper',
          url='http://github.com/mancellin/capytaine',
          author='Matthieu Ancellin',
          author_email='matthieu.ancellin@ucd.ie',
          license='GPLv3',
          packages=['capytaine', 'capytaine.tools', 'capytaine.geometric_bodies'],
          install_requires=[
              'attrs',
              'numpy',
              'scipy',
              'xarray',
              'meshmagick',
          ],
          entry_points={
              'console_scripts': [
                  'capytaine=capytaine.cli:main',
              ],
          },
          ext_modules=[
              Green,
              Wavenumber
          ],
          )
