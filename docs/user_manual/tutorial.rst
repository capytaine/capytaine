========
Tutorial
========

Main concepts
=============

**Mesh**:
    The mesh of a floating body in its averaged position. It is stored as a :code:`Mesh` class
    from meshmagick.

    The mesh is defined as a list of vertices (a vertex is a triplet of real-valued coordinates)
    and a list of faces (a face is a quadruplet of indices of vertices). By default, faces are
    assumed to be quadrangular. Triangular faces are supported as a quadrangle with two identical
    vertices.

    The :code:`Mesh` class also stores some data computed from the vertices and the faces: the
    faces normals, the faces centers and the faces areas.

**Dof**:
    A degree of freedom (or dof) defines a small motion or a deformation of the floating body
    around its averaged position. It is stored as a vector at the center of each faces of the mesh.

    .. Rigid-body dofs can be generated with the :code:`add_translation_dof` and
       :code:`add_rotation_dof` methods.

**FloatingBody**:
    A :code:`FloatingBody` is simply the reunion of a :code:`Mesh` and some degrees of freedom.

**LinearPotentialFlowProblem**:
    A problem is a collection of several parameters: a :code:`FloatingBody`, the wave frequency
    :math:`\omega`,the water depth :math:`h`, the water density :math:`\rho` and the gravity
    acceleration :math:`g`.
    
    The abstract class :code:`LinearPotentialFlowProblem` has two child classes:
    :code:`RadiationProblem` (that requires also the name of the dof that is radiating) and
    :code:`DiffractionProblem` (that requires the angle of the incoming wave field :math:`\beta`).

    Most of the parameters are optionals. A default value is used when they are not provided.

**Solver**
    The core of the code. It has a :code:`solve` method that takes a
    :code:`LinearPotentialFlowProblem` as input and returns a :code:`LinearPotentialFlowResult`.

**LinearPotentialFlowResult**:
    The class storing the results is similar to the class storing a problem, with some
    supplementary data such as :code:`result.added_masses` and :code:`result.radiation_dampings`
    for radiation problems and :code:`result.forces` for diffraction problems.

Step-by-step examples
=====================

TODO
