=============
Floating body
=============

As described in the :doc:`tutorial`, a floating body is defined as a mesh with degrees of freedom and optionally other properties, such as a mass and a center of mass.

Initialization
--------------

The floating body is set up by initializing the
:class:`~capytaine.bodies.bodies.FloatingBody` class::

    body = cpt.FloatingBody(mesh=mesh, dofs={})

The above example creates a new body without degrees of freedom.

Mesh
~~~~

The ``mesh`` can be a :class:`~capytaine.meshes.meshes.Mesh` object, or any of
the variants defined in Capytaine, such as a
:class:`~capytaine.meshes.symmetric.ReflectionSymmetricMesh`.
Meshes from `meshio` can also be given directly to the ``FloatingBody``
constructor without calling :func:`~capytaine.io.meshio.load_from_meshio`.

Lid mesh
~~~~~~~~

For irregular frequencies removal, a second mesh can be provided when defining
the body. It is meant to be a mesh of a lid of the part of the free surface
that is inside the body. The lid can either on the free surface or slightly
below. This mesh is ignored for many computations (such as hydrostatics) and is
only used when solving a BEM problem.

It is set as in the following example::

    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh)

Once a lid mesh has been defined, it is automatically used for irregular
frequencies removal without any other action from the user.

For irregular frequencies removal, the lid is expected to have normals going
down (towards the fluid, through the body). If the lid appears to be horizontal
with normal going up, Capytaine will switch the direction of its normal
vectors.

Currently, meshes with a symmetry are not supported, in the sense that the
computation will be done without using the symmetries when a lid is added. This
should be improved to support at least vertical symmetry plane in a future
version.

Dofs
~~~~

The degrees of freedom are defined as a Python dictionary associating the name
of each dof to a Numpy array of shape ``(nb_faces, 3)``.
This array stores the displacement vector at the center of each face of the
mesh::

   body = cpt.FloatingBody(
           mesh=mesh,
           dofs={
               "heave": np.array([(0, 0, 1) for x, y, z in mesh.faces_centers]),
               "x-shear": np.array([(np.cos(np.pi*z/2), 0, 0) for x, y, z in mesh.faces_centers])
               },
           )

:meth:`cpt.rigid_body_dofs() <~capytaine.bodies.dofs.rigid_body_dofs>` can
be used to automatically give a body the six degrees of freedom of a rigid
body::

   body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -1)))
   print(body.dofs.keys())
   # dict_keys(['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'])

Other parameters
~~~~~~~~~~~~~~~~

Besides the mandatory ``mesh`` and ``dofs``, a floating body can also be
defined with a ``mass``, given a floating point number and a
``center_of_mass``, given a three coordinates.
These two arguments are required for :doc:`hydrostatics` but not for
first-order wave-structure interaction.
Their only role then is for the definition of the rotation degrees of freedom.
When defining a rotation dof, the code looks for attributes called
:code:`rotation_center`, :code:`center_of_mass` or :code:`geometric_center` (in
that order), and use them to define the rotation axis.
If none of them are define, the rotation is defined around the origin of
the domain :math:`(0, 0, 0)`.

Finally, as the mesh objects, the floating body can be assigned a name.


Display and animation
---------------------

The methods :meth:`~capytaine.bodies.bodies.FloatingBody.show()` and
:meth:`~capytaine.bodies.bodies.FloatingBody.show_matplotlib()` of meshes can
also be used on ``FloatingBody``.

Once a :code:`FloatingBody` with dofs has been defineds, the
:meth:`~capytaine.bodies.bodies.FloatingBody.animate`
method can be used to visualize a given motion of the body::

    anim = body.animate(motion={"Heave": 0.1, "Surge": 0.1j}, loop_duration=1.0)
    anim.run()

The above example will present an interactive animation of the linear combination of heave and surge.

Jupyter notebooks can also include a (non-interactive) video of the animation::

    anim.embed_in_notebook(camera_position=(-1.0, -1.0, 1.0), resolution=(400, 300))


Geometric transformations
-------------------------

All the geometric transformation defined on meshes in :doc:`mesh` can also be
applied to ``FloatingBody``. Beside updating the mesh, they also update the
definition of the degrees of freedom and the center of mass (if relevant).


Multiple bodies
---------------

Multiple bodies problems can be defined by combining several bodies with the ``join_bodies`` method::

    all_bodies = cpt.FloatingBody.join_bodies(body_1, body_2, body_3, body_4)

For two-body problems, the ``+`` operator can also be used::

   two_bodies = body_1 + body_2

But it is not recommended to use it for large number of bodies as it is not
strictly associative (that is ``body_1 + (body_2 + body_3)`` has some internal
differences with ``(body_1 + body_2) + body_3``).

When two floating bodies with dofs are merged, the resulting body inherits from
the dofs of the individual bodies with the new name :code:`body_name__dof_name`::

    print(two_bodies.nb_dofs)
    # 12
    print(two_bodies.dofs.keys())
    # dict_keys(['body_1__Surge', 'body_1__Sway', 'body_1__Heave', 'body_1__Roll', 'body_1__Pitch', 'body_1__Yaw', 'body_2__Surge', 'body_2__Sway', 'body_2__Heave', 'body_2__Roll', 'body_2__Pitch', 'body_2__Yaw'])

Capytaine also include helper functions to create arrays of identical bodies::

    array = body.assemble_regular_array(distance=1.0, nb_bodies=(4, 5))

places copies of the ``body`` on a regular grid of :math:`4 \times 5` with distance between bodies of 1 meter, and::

    locations = np.array([[0.0, 0.0], [1.0, 2.0], [3.0, 4.5], [3.0, -0.5]])
    array = body.assemble_arbitrary_array(locations)

places copies of the ``body`` at the list of locations specified.

.. warning::
   As currently implemented in Capytaine, the multiple bodies are stored as a
   single body with a non-connex mesh and generalized degrees of freedom.
   Hence some information about the individual bodies is lost.
   It includes the center of mass and the center of rotation of the individual
   bodies (although the latter could be recovered indirectly by studying the
   definition of the rotation dof).
   Although it does not affect first order wave-structure interaction, it
   hinders the computation of hydrostatics for multiple rigid bodies and will
   need to be fixed in the future.
