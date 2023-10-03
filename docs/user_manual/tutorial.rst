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

Launch an interactive Python console such as :code:`ipython` and import the Capytaine package::

    import capytaine as cpt

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

    sphere = cpt.mesh_sphere(radius=1.0, center=(0, 0, -2), name="my sphere")

Users can also import mesh from various file formats as shown in the :doc:`mesh`
section of the documentation. The mesh is stored as a
:class:`~capytaine.mesh.mesh.Mesh` object. You can for instance access of
coordinates of some of the vertices, faces centers or faces normal vectors using
the following syntax::

    sphere.vertices[:10]  # First ten vertices.
    sphere.faces_centers[5]  # Center of the sixth face (Python arrays start at 0).
    sphere.faces_normals[5]  # Normal vector of the sixth face.

If `vtk` has been installed, the mesh can be displayed in 3D using::

    sphere.show()

Defining a floating body
------------------------

Before solving a diffraction or radiation problem, we need to define the degrees of freedom (dofs) of our body.
In Capytaine, this is done by creating a :code:`FloatingBody` object::

    body = cpt.FloatingBody(mesh=sphere,
                            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -2)),
                            center_of_mass=(0, 0, -2))

The new body defined here will have the six degrees of freedom of a rigid body.
The :code:`rotation_center` is used for the definition of the rotation dofs.
The :code:`center_of_mass` is used for some hydrostatics properties but not required for the diffraction-radiation problems.

The degrees of freedoms are stored in the :code:`dofs` dictionary. To access the name of the dofs of a body, you can use for instance::

    print(body.dofs.keys())
    # dict_keys(['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'])

Dofs can also be defined manually, for instance to model a flexible body. For this purpose, one has to define a list a vectors where each vector is the displacement of the body at the center of a face::

    import numpy as np
    body.dofs["x-shear"] = [(np.cos(np.pi*z/2), 0, 0) for x, y, z in sphere.faces_centers]

    print(body.dofs.keys())
    # dict_keys(['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw', 'x-shear'])

Hydrostatics
------------

Capytaine can directly perform some hydrostatic computations. You can get parameters such as volume, wet surface area, waterplane area, center of buoyancy, metacentric radius and height, hydrostatic stiffness and inertia matrix for any given :code:`FloatingBody`::

    hydrostatics = body.compute_hydrostatics(rho=1025.0)

    print(hydrostatics["disp_volume"])
    # 3.82267415555807

    print(hydrostatics["hydrostatic_stiffness"])
    # <xarray.DataArray 'hydrostatic_stiffness' (influenced_dof: 7, radiating_dof: 7)>
    # [...]
    # Coordinates:
    #   * influenced_dof  (influenced_dof) <U7 'Surge' 'Sway' ... 'Yaw' 'x-shear'
    #   * radiating_dof   (radiating_dof) <U7 'Surge' 'Sway' ... 'Yaw' 'x-shear'

    print(hydrostatics["inertia_matrix"])
    # <xarray.DataArray 'inertia_matrix' (influenced_dof: 7, radiating_dof: 7)>
    # [...]
    # Coordinates:
    #   * influenced_dof  (influenced_dof) <U7 'Surge' 'Sway' ... 'Yaw' 'x-shear'
    #   * radiating_dof   (radiating_dof) <U7 'Surge' 'Sway' ... 'Yaw' 'x-shear'

The matrices here are :math:`7 \times 7` matrices as we have defined seven dofs for our sphere.
The matrices are stored as :code:`DataArray` from the `xarray <https://xarray.dev/>`_ package (see below for an example of usage).
Note that the inertia matrix can only be computed for rigid bodies (assuming constant density). The matrix was filled with :code:`NaN` for the generalized dof :code:`x-shear`.


Defining linear potential flow problems.
----------------------------------------

Let us define a radiation problem for the heave of our sphere::

    from numpy import infty
    problem = cpt.RadiationProblem(body=body, radiating_dof="Heave", omega=1.0, water_depth=infty, g=9.81, rho=1000)

The argument :code:`radiating_dof` must be the name of one of the dofs of the floating body given as the
:code:`body` argument. The wave angular frequency has been set arbitrarily as :math:`\omega = 1 \, \text{rad/s}`.
The water depth is infinite, the gravity acceleration is :math:`g = 9.81 \, \text{m/s}^2` and the water density has
been chosen as :math:`\rho = 1000 \, \text{kg/m}^3`. These last parameters are actually optional.
Since we are using their default value, we could have defined the radiation problem as::

    problem = cpt.RadiationProblem(body=body, radiating_dof="Heave", omega=1.0)

Besides, one can give a :code:`period`, a :code:`wavelength` or a :code:`wavenumber` to specify the frequency::

    problem = cpt.RadiationProblem(body=body, radiating_dof="Heave", wavelength=60.0)

Some more parameters are automatically computed, such as::

    print(problem.wavenumber)
    # 0.1019367991845056
    print(problem.period)
    # 6.283185307179586

Capytaine also implement a :code:`DiffractionProblem` class which does not take a :code:`radiating_dof` argument but instead requires a :code:`wave_direction` in radians::

    diffraction_problem = cpt.DiffractionProblem(body=body, wave_direction=np.pi/2, omega=1.0)

Solve the problem
-----------------

Let us initialize the BEM solver::

    solver = cpt.BEMSolver()

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
    # {'Surge': -2.2737367544323206e-13, 'Sway': 2.8421709430404007e-13, 'Heave': 2207.842170942399, 'Roll': -5.3290705182007514e-14, 'Pitch': 1.3289369604763124e-13, 'Yaw': -4.0933040491086855e -15, 'x-shear': 1.3677947663381929e-13}

The :code:`added_masses` dictionary stores the resulting force on each of the "influenced dofs" of the body.
In this example, the radiating dof is heave and the reaction force in the
:math:`x` direction (:code:`result.added_masses['Surge']`) is negligible with
respect to the one in the :math:`z` direction
(:code:`result.added_masses['Heave']`).

::

    print(result.radiation_dampings)
    # {'Surge': -9.103828801926284e-15, 'Sway': 0.0, 'Heave': 13.623038091677039, 'Roll': -1.7486012 637846216e-15, 'Pitch': 3.4937330806172895e-15, 'Yaw': 9.7897500502683e-17, 'x-shear': 7.271960811294775e-15}

The same thing hold for the diffraction problem::

    diffraction_result = solver.solve(diffraction_problem)
    print(diffraction_result.forces)
    # {'Surge': (1.2789769243681803e-13+2.1760371282653068e-13j), 'Sway': (5.9629097739514805-1923.8976950141728j), 'Heave': (-1802.7102076027684-10.957655022820937j), 'Roll': (-0.010382265075485009+4.15060693393701j), 'Pitch': (-7.799316747991725e-14+2.8310687127941492e-14j), 'Yaw': (-1.2477555116385028e-15+7.261633298807041e-14j), 'x-shear': (-8.43769498715119e-14-2.0972112935169207e-13j)}


Gather results in arrays
------------------------

Let us compute the added mass and radiation damping for all the dofs of our body::

    all_radiation_problems = [cpt.RadiationProblem(body=body, radiating_dof=dof, omega=1.0) for dof in body.dofs]
    all_radiation_results = solver.solve_all(all_radiation_problems)

Here, we used :code:`solve_all` instead of :code:`solve` since we are passing a
list of problems and not a single one. Note that this resolution should be
faster than the first one. The solver has stored some intermediate data for
this body at this wave frequency and will reuse it to solve the new problems.

The results can be gathered together as follow::

    dataset = cpt.assemble_dataset([diffraction_result] + all_radiation_results)

The new object is a NetCDF-like dataset from the xarray package. It is storing the added mass and
radiation damping from the result objects in an organized way. In our example, it is basically two
7×7 matrices. The matrices can be accessed for instance in the following way::

    dataset['added_mass'].sel(radiating_dof=["Surge", "Heave"], influenced_dof=["Surge", "Heave"], omega=1.0)

You'll probably want to solve problems for a wide range of parameters without
defining each test individually. This can be done with the :code:`fill_dataset`
method of the solver. See :doc:`problem_setup`.
