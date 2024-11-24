# Capytaine: a linear potential flow BEM solver with Python.

![CI status](https://github.com/capytaine/capytaine/actions/workflows/test_new_commits.yaml/badge.svg?event=push)
![CI status](https://github.com/capytaine/capytaine/actions/workflows/test_with_latest_dependencies.yaml/badge.svg)


Capytaine is Python package for the simulation of the interaction between water waves and floating bodies in frequency domain.
It is built around a full rewrite of the open source Boundary Element Method (BEM) solver Nemoh for the linear potential flow wave theory.

## Installation

[![PyPI](https://img.shields.io/pypi/v/capytaine)](https://pypi.org/project/capytaine)
[![Conda-forge](https://img.shields.io/conda/vn/conda-forge/capytaine)](https://github.com/conda-forge/capytaine-feedstock)

Packages for Windows, macOS and Linux are available on PyPI:

```bash
pip install capytaine
```
and Conda-forge

```bash
conda install -c conda-forge capytaine
```

or as a standalone executable (with some drawbacks) at https://github.com/capytaine/capytaine-standalone.

## Documentation

[https://capytaine.github.io/](https://capytaine.github.io/)

[![DOI](http://joss.theoj.org/papers/10.21105/joss.01341/status.svg)](https://doi.org/10.21105/joss.01341)

## License

Copyright (C) 2017-2024, Matthieu Ancellin

Since April 2022, the development of Capytaine is funded by the Alliance for Sustainable Energy, LLC, Managing and Operating Contractor for the National Renewable Energy Laboratory (NREL) for the U.S. Department of Energy.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

It is based on version 2 of [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp), which has been developed by Gérard Delhommeau, Aurélien Babarit et al., (École Centrale de Nantes) and was distributed under the Apache License 2.0.

Some core Fortran routines of Capytaine coming from Nemoh version 2 are also available under the Apache License 2.0. They can be found in the [`capytaine/green_functions/libDelhommeau`](https://github.com/capytaine/capytaine/tree/master/capytaine/green_functions/libDelhommeau) directory of Capytaine's repository.

Capytaine includes code from [meshmagick](https://github.com/LHEEA/meshmagick/) by François Rongère (École Centrale de Nantes), licensed under the GNU General Public License (GPL).
