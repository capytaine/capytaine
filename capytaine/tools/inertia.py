#!/usr/bin/env python
# coding: utf-8

# TODO: add the possibility to compute the inertia to an other point than [0, 0, 0]
def eval_plain_mesh_inertias(mesh, rho_medium=1023.):
    """Evaluates the mesh inertia under the assumption of an enclosed volume made of an homogeneous medium of the given density.

    Parameters
    ----------
    rho_medium : float, optional
        The medium density (kg/m**3). Default is 1023 kg.m**3 (salt water)

    Returns
    -------
    RigidBodyInertia
        The mesh inertia instance expressed at origin (0, 0, 0)
    """
    # TODO: allow to specify an other point for inertia matrix expression
    # TODO: manipuler plutot un objet inertia --> creer une classe !
    rho_medium = float(rho_medium)

    volume = mesh.volume
    mass = rho_medium * volume

    integrals = mesh.get_surface_integrals()[6:15]
    sigma_6_8 = integrals[:3]

    normals = mesh.faces_normals.T

    cog = (normals * sigma_6_8).sum(axis=1) / (2*volume)

    sigma9, sigma10, sigma11 = (normals * integrals[3:6]).sum(axis=1)
    sigma12, sigma13, sigma14 = (normals * integrals[6:10]).sum(axis=1)

    xx = rho_medium * (sigma10 + sigma11) / 3.
    yy = rho_medium * (sigma9 + sigma11) / 3.
    zz = rho_medium * (sigma9 + sigma10) / 3.
    xy = rho_medium * sigma12 / 2.
    xz = rho_medium * sigma14 / 2.
    yz = rho_medium * sigma13 / 2.

    return RigidBodyInertia(mass, cog, xx, yy, zz, yz, xz, xy, point=[0, 0, 0])

def eval_shell_mesh_inertias(mesh, rho_medium=7850., thickness=0.02):
    """Evaluates the mesh inertia under the assumption of an enclosed volume made of an homogeneous medium of the
    given density.

    Parameters
    ----------
    rho_medium : float, optional
        The medium density (kg/m**3). Default is 7850 kg/m**3 (Steel density)
    thickness : flaot, optional
        The hull thickness (m). Default is 0.02 m.

    Returns
    -------
    RigidBodyInertia
        The mesh inertia instance expressed at origin (0, 0, 0)
    """
    rho_medium = float(rho_medium)
    thickness = float(thickness)
    surf_density = rho_medium * thickness

    surface = mesh.faces_areas.sum()
    mass = surf_density * surface

    s0, s1, s2, s3, s4, s5, s6, s7, s8 = mesh.get_surface_integrals()[:9].sum(axis=1)

    cog = np.array([s0, s1, s2], dtype=np.float) / mass

    xx = surf_density * (s7 + s8)
    yy = surf_density * (s6 + s8)
    zz = surf_density * (s6 + s7)
    yz = surf_density * s3
    xz = surf_density * s4
    xy = surf_density * s5

    return RigidBodyInertia(mass, cog, xx, yy, zz, yz, xz, xy, point=[0, 0, 0])


