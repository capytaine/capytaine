========
Tutorial
========

Main concepts
=============

:class:`~capytaine.meshes.meshes.Mesh`
    The mesh of a floating body in its averaged position. It is stored as a
    instance of the :code:`Mesh` class.

    The mesh is defined as a list of vertices (a vertex is a triplet of real-valued coordinates)
    and a list of faces (a face is a quadruplet of indices of vertices). By default, faces are
    assumed to be quadrangular. Triangular faces are supported as quadrangles with two identical
    vertices.

    The :code:`Mesh` class also stores some data computed from the vertices and the faces such as
    the faces normals, the faces centers and the faces areas.

**Dof**
    A degree of freedom (or dof) defines a small motion or deformation of the floating body
    around its average position. It is stored as a vector at the center of each faces of the mesh.

    Degrees of freedom appears in two forms in the code:
    :code:`radiating_dof` denotes an actual motion of the body, whereas
    :code:`influenced_dof` denotes a component of a (generalized) force.

    .. note:: For mathematicians in the field of Galerkin Boundary Element Method, the concept
        of degree of freedom might have a different meaning (a basis function of the Galerkin
        decomposition). Here, the degrees of freedom are the physical degrees of freedom of the
        floating body, typically the rigid body translations and rotations.

:class:`~capytaine.bodies.bodies.FloatingBody`
    A :code:`FloatingBody` is mainly the reunion of a :code:`Mesh` and some degrees of freedom.

    The degree of freedom of the body are referred by a name (e.g. `Heave`).
    They should stay in the order in which they have been defined, but `the code
    does not strictly guarantee it <https://github.com/capytaine/capytaine/issues/4>`_.
    Accessing them by name rather than by index should be preferred.

    Beside the mesh and the dofs, some other physical information can be
    stored in a :code:`FloatingBody` instance, such as the mass and the
    position of the center of mass. This information is only required for
    some specific actions (see :doc:`hydrostatics`) and can be left unspecified
    in many cases.

:class:`~capytaine.bem.problems_and_results.LinearPotentialFlowProblem`
    A problem is a collection of several parameters: a :code:`FloatingBody`, the wave angular frequency
    :math:`\omega`, the water depth :math:`h`, the water density :math:`\rho` and the gravity
    acceleration :math:`g`.

    The abstract class :code:`LinearPotentialFlowProblem` has two child classes:
    :class:`~capytaine.bem.problems_and_results.RadiationProblem` (that requires also the name of the dof that is radiating) and
    :class:`~capytaine.bem.problems_and_results.DiffractionProblem` (that requires the angle of the incoming wave field :math:`\beta`).

    Most of the parameters are optional. A default value is used when they are not provided (see the page :doc:`problem_setup`).

:class:`Solver <capytaine.bem.solver.BEMSolver>`
    The core of the code. It has a :meth:`~capytaine.bem.solver.BEMSolver.solve` method that takes a
    :code:`LinearPotentialFlowProblem` as input and returns a :code:`LinearPotentialFlowResult`.
    It calls a class computing the Green function and a class to build the matrices.
    See :doc:`resolution` for details.

:class:`~capytaine.bem.problems_and_results.LinearPotentialFlowResult`
    The class storing the results is similar to the class storing a problem, with some
    supplementary data such as :code:`result.added_masses` and :code:`result.radiation_dampings`
    for radiation problems and :code:`result.forces` for diffraction problems.
    The forces are stored as a dictionary associating the name of a degree of freedom to a value.
    The value is the integral of the force along this degree of freedom.
    For example, to retrieve the components of the force vector on a rigid body in Cartesian coordinates, check the
    value of the force with respect to :code:`Surge`, :code:`Sway` and :code:`Heave`.

Step-by-step example
====================

Launch an interactive Python console such as :code:`ipython`.
All the main features of Capytaine can be loaded with::

    from capytaine import *

Note that Capytaine uses the logging module from Python. Then, you can optionally get some feedback from the code
by initializing the logging module with the following commands::

    import logging
    logging.basicConfig(level=logging.INFO)

Replace :code:`INFO` by :code:`DEBUG` to get more information on everything that is happening
inside the solver. On the other hand, if you set the level to :code:`WARNING`, only important
warnings will be printed out by the solver (this is the default behavior when the logging module
has not been set up). 

Load a mesh
-----------

For this tutorial we will use one of the mesh generators included into Capytaine for simple
geometric shapes::

    sphere = Sphere(radius=1.0, center=(0, 0, -2), name="my buoy")

Users can also import mesh from various file formats as shown in the :doc:`mesh`
section of the documentation. The mesh is stored as a
:class:`~capytaine.mesh.mesh.Mesh` object. You can for instance access of
coordinates of some of the vertices, faces centers or faces normal vectors using
the following syntax::

    sphere.mesh.vertices[:10]  # First ten vertices.
    sphere.mesh.faces_centers[5]  # Center of the sixth face (Python arrays start at 0).
    sphere.mesh.faces_normals[5]  # Normal vector of the sixth face.

If `vtk` has been installed, the mesh can be displayed in 3D using::

    sphere.show()

Defining dofs
-------------

Before solving a diffraction or radiation problem, we need to define the degrees of freedom (dofs) of our
body. It can be done in several ways:

* The manual way: define a list a vectors where each vector is the displacement of the
  body at the center of a face. The example below is the simplest example of a rigid body motion in
  the :math:`x` direction::

    sphere.dofs['Surge'] = [(1, 0, 0) for face in sphere.mesh.faces]

* Helpers functions are available to define rigid body translations and rotations. For instance for
  the motion in the :math:`z` direction, we can use :meth:`FloatingBody.add_translation_dof <capytaine.bodies.bodies.FloatingBody.add_translation_dof>`.
  It can recognize some dof names such as "Surge", "Sway" and "Heave"::

    sphere.add_translation_dof(name="Heave")

  See the documentation of :meth:`FloatingBody.add_rotation_dof <capytaine.bodies.bodies.FloatingBody.add_rotation_dof>` and :meth:`FloatingBody.add_all_rigid_body_dofs <capytaine.bodies.bodies.FloatingBody.add_all_rigid_body_dofs>`.

The degrees of freedoms are stored in the :code:`dofs` dictionary. To access the name of the dofs of a
body, you can use for instance::

    print(sphere.dofs.keys())
    # dict_keys(['Surge', 'Heave'])

Hydrostatics
------------

Capytaine can directly perform some hydrostatic computation for a given mesh. You can get parameters such as volume, wet surface area, waterplane area, center of buoyancy, metacentric radius and height, hydrostatic stiffness and interia mass for any given :code:`FloatingBody`.

Let us give the code some more information about the body::

    sphere.center_of_mass = (0, 0, -2)
    sphere.rotation_center = (0, 0, -2)

The "rotation center" is the point used to define the rotation dofs.
(Due to a current limitation of the hydrostatics methods, the definition of the rotation center is required as soon as there is a rigid body dof — here surge and heave —, even if it is not a rotation dof.)

Each hydrostatic parameter can be computed by a dedicated method::

    print(sphere.volume)
    # 3.82267415555807

    print(sphere.center_of_buoyancy)
    # [-3.58784373e-17 -2.59455034e-17 -2.00000000e+00]

    print(sphere.compute_hydrostatic_stiffness())
    # <xarray.DataArray 'hydrostatic_stiffness' (influenced_dof: 2, radiating_dof: 2)>
    # array([[0.00000000e+00, 0.00000000e+00],
    #        [0.00000000e+00, 2.38246922e-13]])
    # Coordinates:
    #   * influenced_dof  (influenced_dof) <U5 'Surge' 'Heave'
    #   * radiating_dof   (radiating_dof) <U5 'Surge' 'Heave'

    print(sphere.compute_rigid_body_inertia())
    # <xarray.DataArray 'inertia_matrix' (influenced_dof: 2, radiating_dof: 2)>
    # array([[3822.67415556,    0.        ],
    #        [   0.        , 3822.67415556]])
    # Coordinates:
    #   * influenced_dof  (influenced_dof) <U5 'Surge' 'Heave'
    #   * radiating_dof   (radiating_dof) <U5 'Surge' 'Heave'

The matrices here are :math:`2 \times 2` matrices as we have defined only two dofs for our sphere.

You can also use :code:`compute_hydrostatics` method which computes all hydrostatic parameters and returns a :code:`dict` of parameters and values::

    hydrostatics = sphere.compute_hydrostatics()

.. note::
   Before computing hydrostatic parameters, you might want to crop your mesh using `body.keep_immersed_part()`.
   It is not required here since the sphere is fully immersed.
   Cropping is included in the :code:`compute_hydrostatics()` function.


Defining linear potential flow problems.
----------------------------------------

Let us define a radiation problem for the heave of our sphere::

    from numpy import infty
    problem = RadiationProblem(body=sphere, radiating_dof="Heave", omega=1.0, sea_bottom=-infty, g=9.81, rho=1000)

The argument :code:`radiating_dof` must be the name of one of the dofs of the floating body given as the
:code:`body` argument. The wave angular frequency has been set arbitrarily as :math:`\omega = 1 \, \text{rad/s}`.
The water depth is infinite, the gravity acceleration is :math:`g = 9.81 \, \text{m/s}^2` and the water density has
been chosen as :math:`\rho = 1000 \, \text{kg/m}^3`. These last parameters are actually optional.
Since we are using their default value, we could have defined the radiation problem as::

    problem = RadiationProblem(body=sphere, radiating_dof="Heave", omega=1.0)

Some more parameters are automatically computed, such as::

    print(problem.wavenumber)
    # 0.1019367991845056
    print(problem.period)
    # 6.283185307179586

Solve the problem
-----------------

Let us initialize the BEM solver::

    solver = BEMSolver()

Solver settings could have been given at this point, but in this tutorial, we will use the default settings.
Let us now solve the problem we defined earlier::

    result = solver.solve(problem)

The :meth:`~capytaine.bem.solver.BEMSolver.solve` method returns a result object. The result object contains all of the data from
the problem it comes from::

    print(result.omega)
    # 1.0
    print(result.body.name)
    # "my buoy"
    print(result.radiating_dof)
    # "Heave"
    print(result.period)
    # 6.283185307179586

Of course, it also stores some output data. Since we solved a radiation problem, we can now access
the added mass and radiation damping::

    print(result.added_masses)
    # {'Surge': 9.154531598110083e-06, 'Heave': 2207.8423200090374}

The :code:`added_masses` dictionary stores the resulting force on each of the "influenced dofs" of the body.
In this example, the radiating dof is heave and the reaction force in the
:math:`x` direction (:code:`result.added_masses['Surge']`) is negligible with
respect to the one in the :math:`z` direction
(:code:`result.added_masses['Heave']`).

::

    print(result.radiation_dampings)
    # {'Surge': -5.792518686098536e-07, 'Heave': 13.62318484050783}

Gather results in arrays
------------------------

Let us compute the added mass and radiation damping for surge::

    other_problem = RadiationProblem(body=sphere, radiating_dof="Surge", omega=1.0)
    other_result = solver.solve(other_problem)

Note that this second resolution should be faster than the first one. The solver has stored some
intermediate data for this body and will reuse them to solve this other problem.

The results can be gathered together as follow::

    dataset = assemble_dataset([result, other_result])

The new object is a NetCDF-like dataset from the xarray package. It is storing the added mass and
radiation damping from the result objects in an organized way. In our example, it is basically two
2x2 matrices. The matrices can be accessed for instance in the following way::

    dataset['added_mass'].sel(radiating_dof=["Surge", "Heave"], influenced_dof=["Surge", "Heave"], omega=1.0)

You'll probably want to solve problems for a wide range of parameters without
defining each test individually. This can be done with the :code:`fill_dataset`
method of the solver. See :doc:`problem_setup`.

