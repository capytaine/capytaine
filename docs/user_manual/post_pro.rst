===============
Post-processing
===============

Free surface elevation
----------------------

To compute the free surface elevation, let us first initialize a
:class:`~capytaine.post_pro.free_srufaces.FreeSurface` object::

    from capytaine import FreeSurface
    fs = FreeSurface(x_range=(-10, 10), nx=10, y_range=(-5, 5), ny=10)

The above code generates a regular free surface mesh of :math:`10 \times 10`
cells. This object can be used by the solver with a results object to get the
free surface elevation::

    fs_elevation = solver.get_free_surface_elevation(result, free_surface)

The output is a numpy array storing the free surface elevation in frequency
domain as a complex number at each point of the free surface (in the present
example an array of shape :code:`(10, 10)`).

The result object should have been computed with the option
:code:`keep_details=True`. The solver does not need to be the one that computed
the result object.

The undisturbed incoming waves (Airy waves) can be computed as follow::

    incoming_waves = fs.incoming_waves(DiffractionProblem(omega=1.0, angle=pi/2))

See the examples in the :doc:`cookbook` for usage in a 3D animation.
