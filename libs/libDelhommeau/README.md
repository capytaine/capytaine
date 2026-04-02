# LibDelhommeau

A Fortran 90 library for the computation of the Green function for linear potential flow BEM problems, using Delhommeau's method as in Aquaplus, [Nemoh](https://lheea.ec-nantes.fr/valorisation/logiciels-et-brevets/nemoh-presentation) and [Capytaine](https://joss.theoj.org/papers/10.21105/joss.01341).

This Fortran repository contains the source code of the compiled extensions of the Python software [Capytaine](https://github.com/capytaine/capytaine).
One of the reason to have it in a separate repository is because libDelhommeau is distributed with a more permissive license than Capytaine.
Another reason to decouple it from Capytaine is to support reuse of the code in other software.

However, unless you are an experienced user, you might be more interested in using the full [Capytaine](https://github.com/capytaine/capytaine) software which is easier to use and better documented than this library.

## Authors
Matthieu Ancellin, based on the work of Gérard Delhommeau, Aurélien Babarit, and other contributors to the aforementioned codes.

## License
Apache license 2.0
