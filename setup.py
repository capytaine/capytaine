#!/usr/bin/env python
# coding: utf-8

from numpy.distutils.core import Extension, setup

VERSION = '0.5.1'

NemohCore = Extension(
    name="capytaine.NemohCore",
    sources=[
        "capytaine/NemohCore/constants.f90",
        "capytaine/NemohCore/Green_Rankine.f90",
        "capytaine/NemohCore/Initialize_Green_wave.f90",
        "capytaine/NemohCore/Green_wave.f90",
        "capytaine/NemohCore/old_Prony_decomposition.f90",
    ],
    extra_compile_args=['-fopenmp'],
    extra_link_args=['-fopenmp'],
    # # Uncomment the following lines to get more verbose output from f2py.
    # define_macros=[
    #     ('F2PY_REPORT_ATEXIT', 1),
    #     ('F2PY_REPORT_ON_ARRAY_COPY', 1),
    # ],
)

if __name__ == "__main__":
    setup(name='capytaine',
          version=VERSION,
          description='Nemoh python wrapper',
          url='http://github.com/mancellin/capytaine',
          author='Matthieu Ancellin',
          author_email='matthieu.ancellin@ucd.ie',
          license='GPLv3',
          packages=[
              'capytaine',
              'capytaine.mesh',
              'capytaine.geometric_bodies',
              'capytaine.tools',
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
