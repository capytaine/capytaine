# Capytaine: Nemoh's Python interface

Capytaine is a BEM solver prototype for the potential flow linear wave theory, written in Python and Fortran.

It is based on
* [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp)'s Fortran core routines for the computation of the Green function,
* [meshmagick](https://github.com/LHEEA/meshmagick) for the manipulation of meshes,
* and various tools from the [Python scientific stack](https://scipy.org/).

It's main purpose is to give developers and researchers an easier access to the inside gears of Nemoh.

For users, it can be used in a similar way as Nemoh's Matlab interface. Due to some experimental optimizations, it is significantly faster than Nemoh 2.0. See `example/cylinder.py` for an example of configuration and usage.

## Installation

[An installation guide is available in the wiki.](https://github.com/mancellin/capytaine/wiki/Installation)
Capytaine is still under development, please report any bug or lacking feature.

## License

Capytaine is developed by Matthieu Ancellin (University College Dublin).
It is licensed under the GNU General Public License v3.0.

It includes some code from [Nemoh](https://lheea.ec-nantes.fr/logiciels-et-brevets/nemoh-presentation-192863.kjsp), which has been developed by Gérard Delhommeau, Aurélien Babarit et al., (École Centrale de Nantes) and is distributed under the Apache License 2.0.
