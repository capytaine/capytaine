
Multiple bodies
===============

Multiple bodies problems can be defined by combining several bodies in the :class:`~capytaine.bodies.multibodies.Multibody` class::

    all_bodies = cpt.Multibody([body_1, body_2, body_3, body_4])

The resulting object can be used in the same way as the single :doc:`body`.

For two-body problems, the ``+`` operator can also be used::

   two_bodies = body_1 + body_2

But it is not recommended to use it for large number of bodies as it is not
strictly associative (that is ``body_1 + (body_2 + body_3)`` has some minor
internal differences with ``(body_1 + body_2) + body_3``).

The multi-body object inherits the dofs of the individual bodies with the new name :code:`body_name__dof_name`::

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
