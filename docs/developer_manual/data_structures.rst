===============
Data structures
===============

This page discuss some of the implementation details of how data are stored internally in Capytaine.

Simple mesh
-----------

A mesh can be seen as a list of faces.
Data on faces (such as the area of each face) are stored in array of length the number of faces.
The ordering in the array should be consistent for all fields defined on this mesh.

Symmetric mesh
--------------

A reflection symmetric mesh is stored as two halves.
An array storing a data field on the mesh is organised such that the first half of the array is the reference half of the mesh and the second half of the array is the reflected half of the mesh.
Taking advantage of the symmetries in the resolution requires the matrices to have a block circulant structure that would be too inconvenient to implement if the reference faces and the symmetric faces were ordered in a different way.


.. raw::

    |  half 1    |   half 2   |


Body with lid
-------------

Assume a FloatingBody defined with a hull mesh and a lid mesh:

- the dofs of the body are defined only on the hull
- the input Neumann boundary condition of the LinearPotentialFlowProblem, as well as the source and potential outputs fields are defined on the hull and the lid

In general, the lid are the last panels of the boundary condition, but this will not be the case in the following more complicated cases, so to access only the panels on the hull, the attribute ``hull_masks``, stored as a boolean array, should be used.

.. raw::

               |  hull      |  lid  |
    dof:       |------------|
    bc:        |--------------------|

    hull_mask: |TTTTTTTTTTTTTFFFFFFF|


Body with lid and symmetry
--------------------------

For a single body with a lid and symmetry, the last panels are always the symmetric part, so the data about the lid panels are distributed in the array:

.. raw::

               |  hull      |   lid   |  sym hull  | sym lid |
    dof:       |------------|         |------------|
    bc:        |---------------------------------------------|

    hull_mask: |TTTTTTTTTTTTTFFFFFFFFFTTTTTTTTTTTTTTFFFFFFFFF|


Multibody
---------

For multibodies, the dofs in ``Multibody.dofs`` are stored in arrays of the length of the total number of panels of the mesh of all bodies, possibly using the ``DofOnSubmesh`` class.


.. raw::

               |  hull 1    |   hull2    |
    dof body 1:|------------|            |
    dof body 2:             |------------|
    dof:       |-------------------------|
    bodymask 1:|TTTTTTTTTTTTTFFFFFFFFFFFF|
    bodymask 2:|FFFFFFFFFFFFFTTTTTTTTTTTT|

    bc:        |-------------------------|



Multibody with symmetry
-----------------------

For multiple bodies composed of symmetric bodies


.. raw::
               |  hull 1    |   hull2    | sym hull 1 | sym hull 2 |
    dof body 1:|------------|            |------------|
    dof body 2:             |------------|            |------------|
    dof:       |---------------------------------------------------|
    bodymask 1:|TTTTTTTTTTTTTFFFFFFFFFFFFFTTTTTTTTTTTTFFFFFFFFFFFFF|
    bodymask 2:|FFFFFFFFFFFFFTTTTTTTTTTTTFFFFFFFFFFFFFTTTTTTTTTTTTT|

    bc:        |---------------------------------------------------|


Multibody with lid
------------------

.. raw::
               |  hull 1    |  lid 1 |  hull2    | lid 2 |
    dof body 1:|------------|
    dof body 2:                      |-----------|
    dof:       |------------|        |-----------|
    bodymask 1:|TTTTTTTTTTTT|        |FFFFFFFFFFF|
    bodymask 2:|FFFFFFFFFFFF|        |TTTTTTTTTTT|

    bc:        |-----------------------------------------|
    hull_mask: |TTTTTTTTTTTTTFFFFFFFFTTTTTTTTTTTTTFFFFFFF|


 The body from which each lid face originates is not tracked, as there is no use for it.


Multibody with lid and symmetry
-------------------------------

.. raw::
               |  hull 1    | hull 2     | lid 1     | lid 2     | sym hull 1 | sym hull 2 | sym lid 1 | sym lid 2 |
    dof body 1:|------------|                                    |------------|
    dof body 2:             |------------|                                    |------------|
    dof:       |-------------------------|                       |-------------------------|
    bodymask 1:|TTTTTTTTTTTTFFFFFFFFFFFFF|                       |TTTTTTTTTTTTFFFFFFFFFFFFF|
    bodymask 1:|FFFFFFFFFFFFTTTTTTTTTTTTT|                       |FFFFFFFFFFFFTTTTTTTTTTTTT|

    bc:        |---------------------------------------------------------------------------------------------------|
    hull_mask: |TTTTTTTTTTTTTTTTTTTTTTTTTFFFFFFFFFFFFFFFFFFFFFFFFTTTTTTTTTTTTTTTTTTTTTTTTTTTFFFFFFFFFFFFFFFFFFFFFFF|
