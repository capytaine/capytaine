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

Copyright (C) 2017-2026, Matthieu Ancellin and the Capytaine developers.

Since April 2022, the development of Capytaine is funded by the Alliance for Sustainable Energy, LLC, Managing and Operating Contractor for the National Renewable Energy Laboratory (NREL) for the U.S. Department of Energy.

Since April 2025, the development of Capytaine is also funded by [Mews Labs](https://www.mews-labs.com/) and BPI France.

It is based on version 2 of [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp), which has been developed by Gérard Delhommeau, Aurélien Babarit et al., (École Centrale de Nantes) and was distributed under the Apache License 2.0.

Since version 3, Capytaine is licensed under the Apache License, Version 2.0.
You may obtain a copy of the License at  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
