---
title: 'Capytaine: a Python-based linear potential flow solver'
tags:
  - python
  - fortran
  - hydrodynamics
  - potential flow
  - water wave
  - wave energy
authors:
 - name: Matthieu Ancellin
   orcid: 0000-0002-0316-3230
   affiliation: 1
 - name: Frédéric Dias
   orcid: 0000-0002-5123-4929
   affiliation: 1
affiliations:
 - name: UCD School of Mathematics and Statistics, University College Dublin, MaREI Centre, Ireland
   index: 1
date: 19 March 2019
bibliography: paper.bib
nocite: '@folley_numerical_2016; @babarit_ocean_2017; @matplotlib'
---

# Context

One of the first phases of the design of any ship or floating structure is its simulation with a linear potential flow model.
This model gives a description of the interaction between the ocean waves and the floating body [@newman_marine_1977].

Assuming an incompressible, inviscid and irrotational flow, the fluid motion can be modeled by a scalar field $\phi$, called velocity potential and such that $u = \nabla \phi$.
The conservation of the mass of the fluid can be written as a Laplace problem $\nabla^2 \phi = 0$ .
Assuming small motions of the floating body and small wave heights, the boundary condition at the free surface can be linearized.
The resulting problem is fully linear and can be solved in the frequency domain.

The problem can be formulated as a boundary integral problem and solved by the Boundary Element Method (BEM).
In contrast with other BEM problems, the linear potential flow solvers use an expression of the Green's function which accounts for the boundary conditions at the free surface and the sea bottom.
Thus, only the wet surface of the floating body needs to be discretized.
The Green's function can be computed in several ways [@xie_comparison_2018], which are implemented in several software, most of them commercial and closed-source. 

In 2014, the École Centrale de Nantes (ECN) released the code *Nemoh* as an open-source software [@babarit_theoretical_2015].
This solver is based on in-house and commercial codes developed in the ECN since the 1970's [@delhommeau_problemes_1987; @delhommeau_amelioration_1989; @delhommeau_seakeeping_1993].
Since its open-source release, it has been used in numerous applications, in particular by the community of marine renewable energies [@penalbaretes_using_2017].

However, the capabilities of the code are limited in comparison with its commercial counterparts.
Developing new features is difficult since the core of the implementation dating back from the 1970's is poorly readable and the documentation is scarce.
The goal of the present work is the modernization of this code in order to have an open-source linear potential flow solver that is easier to maintain and develop.

# The `capytaine` Python package

The present work is a complete modernization of the Nemoh code, including a refactoring of the core routines and a Python user interface. 

The core routines have been kept in Fortran, although a more modern coding style has been used.
Dead code has been removed, no global variables are used anymore, more meaningful names have been used for functions and variables and a lot of comments have been added.

This new version is meant to be used within the scientific Python ecosystem.
The integration of the Fortran core with Python has been done with F2PY [@f2py].
The naive linear solver included in Nemoh has been replaced by a call to the state-of-the-art solvers from the Numpy and Scipy libraries [@numpy; @scipy].
The Intel MKL library can thus be used to solve the linear system in parallel with OpenMP.
The rest of the code is the independent computation of the coefficients of a matrix and has been straightforwardly parallelized with OpenMP, making most of the code parallel.

The code has also been integrated with other Python libraries.
Reading, writing and transforming meshes is done by the integration of the `meshmagick` [@meshmagick] library.
The default output format is a `xarray.Dataset` [@hoyer_xarray_2017] that can be saved in the standard NetCDF format.
`VTK` [@vtk] is used for the 3D visualization of the meshes and the results.
Testing of the code is done with the `pytest` [@pytest] library.

Efforts have been made to follow best practices for the development of scientific software.
For instance, the default output object includes the inputs and the settings of the solver in order to ease the reproducibility of the results.
However, the present version is a step back with respect to Nemoh regarding the durability of the code.
Indeed, the original code is a single Fortran project with no dependency.
Thus it requires far less maintenance than the present code which is built on top of several layers of fast evolving languages and software.

This more modular version of the code allows the experimentation of new techniques to improve its performance.
Hierarchical matrices have been implemented in the code and their combination with the use of local symmetries of the mesh is being tested [@ancellin_using_2018; @ancellin_hierarchical_2019].
It might have been more efficient to try to merge the code into an existing BEM solver instead of redeveloping advanced BEM techniques.
The present work is nonetheless a step towards a fully modular code, in which the current Green's function evaluation routines could be used in another BEM solver (e.g. [@alouges_fem_2018]), and conversely other methods for the computation of the Green's function [@xie_comparison_2018] could be plugged in to the present user interface.

# Acknowledgment

This work is funded by Science Foundation Ireland (SFI) under Marine Renewable Energy Ireland (MaREI), the SFI Centre for Marine Renewable Energy Research (grant 12/RC/2302).

# References
