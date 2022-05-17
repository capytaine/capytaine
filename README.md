# Capytaine: a linear potential flow BEM solver with Python.

[![DOI](http://joss.theoj.org/papers/10.21105/joss.01341/status.svg)](https://doi.org/10.21105/joss.01341)

Capytaine is Python package for the simulation of the interaction between water waves and floating bodies in frequency domain.
It is built around a full rewrite of the open source Boundary Element Method (BEM) solver Nemoh for the linear potential flow wave theory.

## Installation

On Windows, macOS and Linux, using the [Conda package manager](https://www.anaconda.com/distribution/):

```bash
conda install -c conda-forge capytaine
```

## Documentation

[https://ancell.in/capytaine/latest/](https://ancell.in/capytaine/latest/)

## License

Copyright (C) 2017-2022, Matthieu Ancellin

Since April 2022, the development of Capytaine is funded by the Alliance for Sustainable Energy, LLC, Managing and Operating Contractor for the National Renewable Energy Laboratory (NREL) for the U.S. Department of Energy.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

It is based on [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp), which has been developed by Gérard Delhommeau, Aurélien Babarit et al., (École Centrale de Nantes) and is distributed under the Apache License 2.0.

It includes code from [meshmagick](https://github.com/LHEEA/meshmagick/) by François Rongère (École
Centrale de Nantes), licensed under the GNU General Public License (GPL).
