#!/usr/bin/env python
# coding: utf-8
"""Floating bodies to be used in radiation-diffraction problems."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
import copy
from itertools import chain, accumulate, product, zip_longest

import datetime

import numpy as np
import xarray as xr

from capytaine.meshes.geometry import Abstract3DObject, Plane, inplace_transformation
from capytaine.meshes.meshes import Mesh
from capytaine.meshes.symmetric import build_regular_array_of_meshes
from capytaine.meshes.collections import CollectionOfMeshes

LOG = logging.getLogger(__name__)

TRANSLATION_DOFS_DIRECTIONS = {"surge": (1, 0, 0), "sway": (0, 1, 0), "heave": (0, 0, 1)}
ROTATION_DOFS_AXIS = {"roll": (1, 0, 0), "pitch": (0, 1, 0), "yaw": (0, 0, 1)}


class FloatingBody(Abstract3DObject):
    """A floating body described as a mesh and some degrees of freedom.

    The mesh structure is stored as a Mesh from capytaine.mesh.mesh or a
    CollectionOfMeshes from capytaine.mesh.meshes_collection.

    The degrees of freedom (dofs) are stored as a dict associating a name to
    a complex-valued array of shape (nb_faces, 3). To each face of the body
    (as indexed in the mesh) corresponds a complex-valued 3d vector, which
    defines the displacement of the center of the face in frequency domain.

    Parameters
    ----------
    mesh : Mesh or CollectionOfMeshes, optional
        the mesh describing the geometry of the floating body.
        If none is given, a empty one is created.
    dofs : dict, optional
        the degrees of freedom of the body.
        If none is given, a empty dictionary is initialized.
    mass : float or None, optional
        the mass of the body in kilograms.
        Required only for some hydrostatics computation.
        If None, the mass is implicitly assumed to be the mass of displaced water.
    center_of_mass: 3-element array, optional
        the position of the center of mass.
        Required only for some hydrostatics computation.
    name : str, optional
        a name for the body.
        If none is given, the one of the mesh is used.
    """

    def __init__(self, mesh=None, dofs=None, mass=None, center_of_mass=None, name=None):
        if mesh is None:
            mesh = Mesh(name="dummy_mesh")

        if dofs is None:
            dofs = {}

        if name is None:
            name = mesh.name

        assert isinstance(mesh, Mesh) or isinstance(mesh, CollectionOfMeshes)
        self.mesh = mesh
        self.full_body = None
        self.dofs = dofs
        self.mass = mass
        self.center_of_mass = center_of_mass
        self.name = name

        if self.mesh.nb_vertices == 0 or self.mesh.nb_faces == 0:
            LOG.warning(f"New floating body (with empty mesh!): {self.name}.")
        else:
            self.mesh.heal_mesh()
            LOG.info(f"New floating body: {self.name}.")

    @staticmethod
    def from_meshio(mesh, name=None) -> 'FloatingBody':
        """Create a FloatingBody from a meshio mesh object."""

        import meshio
        if not isinstance(mesh, meshio._mesh.Mesh):
            raise TypeError('mesh must be of type meshio._mesh.Mesh, received {:}'.format(type(mesh)))

        if name is None:
            date_str = datetime.datetime.now().strftime('%Y%m%d%H%M%S%f')
            name = 'fb_{:}'.format(date_str)

        def all_faces_as_quads(cells):
            all_faces = []
            if 'quad' in cells:
                all_faces.append(cells['quad'])
            if 'triangle' in cells:
                num_triangles = len(mesh.cells_dict['triangle'])
                LOG.info("Stored {:} triangle faces as quadrilaterals".format(num_triangles))
                triangles_as_quads = np.empty((cells['triangle'].shape[0], 4), dtype=int)
                triangles_as_quads[:, :3] = cells['triangle'][:, :]
                triangles_as_quads[:, 3] = cells['triangle'][:, 2]  # Repeat one node to make a quad
                all_faces.append(triangles_as_quads)
            return np.concatenate(all_faces)

        cpt_mesh = Mesh(vertices=mesh.points,
                        faces=all_faces_as_quads(mesh.cells_dict),
                        name=name+"_mesh")

        fb = FloatingBody(mesh=cpt_mesh, name=name)
        return fb


    @staticmethod
    def from_file(filename: str, file_format=None, name=None) -> 'FloatingBody':
        """Create a FloatingBody from a mesh file using meshmagick."""
        from capytaine.io.mesh_loaders import load_mesh
        if name is None:
            name = filename
        mesh = load_mesh(filename, file_format, name=f"{name}_mesh")
        return FloatingBody(mesh, name=name)

    def __lt__(self, other: 'FloatingBody') -> bool:
        """Arbitrary order. The point is to sort together the problems involving the same body."""
        return self.name < other.name

    ##########
    #  Dofs  #
    ##########

    @property
    def nb_dofs(self) -> int:
        """Number of degrees of freedom."""
        return len(self.dofs)

    def add_translation_dof(self, direction=None, name=None, amplitude=1.0) -> None:
        """Add a new translation dof (in place).
        If no direction is given, the code tries to infer it from the name.

        Parameters
        ----------
        direction : array of shape (3,), optional
            the direction of the translation
        name : str, optional
            a name for the degree of freedom
        amplitude : float, optional
            amplitude of the dof (default: 1.0 m/s)
        """
        if direction is None:
            if name is not None and name.lower() in TRANSLATION_DOFS_DIRECTIONS:
                direction = TRANSLATION_DOFS_DIRECTIONS[name.lower()]
            else:
                raise ValueError("A direction needs to be specified for the dof.")

        if name is None:
            name = f"dof_{self.nb_dofs}_translation"

        direction = np.asarray(direction)
        assert direction.shape == (3,)

        motion = np.empty((self.mesh.nb_faces, 3))
        motion[:, :] = direction
        self.dofs[name] = amplitude * motion

    def add_rotation_dof(self, axis=None, name=None, amplitude=1.0) -> None:
        """Add a new rotation dof (in place).
        If no axis is given, the code tries to infer it from the name.

        Parameters
        ----------
        axis: Axis, optional
            the axis of the rotation
        name : str, optional
            a name for the degree of freedom
        amplitude : float, optional
            amplitude of the dof (default: 1.0)
        """
        if axis is None:
            if name is not None and name.lower() in ROTATION_DOFS_AXIS:
                axis_direction = ROTATION_DOFS_AXIS[name.lower()]
                for point_attr in ('rotation_center', 'center_of_mass', 'geometric_center'):
                    if hasattr(self, point_attr) and getattr(self, point_attr) is not None:
                        axis_point = getattr(self, point_attr)
                        LOG.info(f"The rotation dof {name} has been initialized around the point: "
                                 f"{self.name}.{point_attr} = {getattr(self, point_attr)}")
                        break
                else:
                    axis_point = np.array([0, 0, 0])
                    LOG.warning(f"The rotation dof {name} has been initialized "
                                f"around the origin of the domain (0, 0, 0).")
            else:
                raise ValueError("A direction needs to be specified for the dof.")
        else:
            axis_point = axis.point
            axis_direction = axis.vector

        if name is None:
            name = f"dof_{self.nb_dofs}_rotation"

        if self.mesh.nb_faces == 0:
            self.dofs[name] = np.empty((self.mesh.nb_faces, 3))
        else:
            motion = np.cross(axis_point - self.mesh.faces_centers, axis_direction)
            self.dofs[name] = amplitude * motion

    def add_all_rigid_body_dofs(self) -> None:
        """Add the six degrees of freedom of rigid bodies (in place)."""
        self.add_translation_dof(name="Surge")
        self.add_translation_dof(name="Sway")
        self.add_translation_dof(name="Heave")
        self.add_rotation_dof(name="Roll")
        self.add_rotation_dof(name="Pitch")
        self.add_rotation_dof(name="Yaw")

    @inplace_transformation
    def keep_only_dofs(self, dofs):
        for dof in list(self.dofs.keys()):
            if dof not in dofs:
                del self.dofs[dof]

        if hasattr(self, 'inertia_matrix'):
            self.inertia_matrix = self.inertia_matrix.sel(radiating_dof=dofs, influenced_dof=dofs)
        if hasattr(self, 'hydrostatic_stiffness'):
            self.hydrostatic_stiffness = self.hydrostatic_stiffness.sel(radiating_dof=dofs, influenced_dof=dofs)

        return self

    def add_dofs_labels_to_vector(self, vector):
        """Helper function turning a bare vector into a vector labelled by the name of the dofs of the body,
        to be used for instance for the computation of RAO."""
        return xr.DataArray(data=np.asarray(vector), dims=['influenced_dof'],
                            coords={'influenced_dof': list(self.dofs)},
                            )

    def add_dofs_labels_to_matrix(self, matrix):
        """Helper function turning a bare matrix into a matrix labelled by the name of the dofs of the body,
        to be used for instance for the computation of RAO."""
        return xr.DataArray(data=np.asarray(matrix), dims=['influenced_dof', 'radiating_dof'],
                            coords={'influenced_dof': list(self.dofs), 'radiating_dof': list(self.dofs)},
                            )

    ###################
    # Hydrostatics #
    ###################

    def surface_integral(self, data, **kwargs):
        """Returns integral of given data along wet surface area."""
        return np.sum(data * self.mesh.faces_areas, **kwargs)

    def waterplane_integral(self, data, **kwargs):
        """Returns integral of given data along water plane area."""
        return self.surface_integral(self.mesh.faces_normals[:,2] * data, **kwargs)

    @property
    def wet_surface_area(self):
        """Returns wet surface area."""
        return self.surface_integral(1)

    @property
    def volumes(self):
        """Returns volumes using x, y, z components of the FloatingBody."""
        norm_coord = self.mesh.faces_normals * self.mesh.faces_centers
        return self.surface_integral(norm_coord.T, axis=1)

    @property
    def volume(self):
        """Returns volume of the FloatingBody."""
        return np.mean(self.volumes)

    def disp_mass(self, *, rho=1000):
        return rho * self.volume

    @property
    def center_of_buoyancy(self):
        """Returns center of buoyancy of the FloatingBody."""
        volume = self.volume
        coords_sq_norm = self.mesh.faces_normals * self.mesh.faces_centers**2
        return self.surface_integral(coords_sq_norm.T, axis=1) / (2*volume)

    @property
    def waterplane_area(self):
        """Returns water plane area of the FloatingBody."""
        waterplane_area = -self.waterplane_integral(1)
        return waterplane_area

    @property
    def waterplane_center(self):
        """Returns water plane center of the FloatingBody.

        Note: Returns None if the FloatingBody is full submerged.
        """
        waterplane_area = self.waterplane_area
        if abs(waterplane_area) < 1e-10:
            return None
        else:
            waterplane_center = -self.waterplane_integral(
                self.mesh.faces_centers.T, axis=1) / waterplane_area
            return waterplane_center[:-1]

    @property
    def transversal_metacentric_radius(self):
        """Returns transversal metacentric radius of the body."""
        inertia_moment = -self.waterplane_integral(self.mesh.faces_centers[:,1]**2)
        return inertia_moment / self.volume

    @property
    def longitudinal_metacentric_radius(self):
        """Returns longitudinal metacentric radius of the body."""
        inertia_moment = -self.waterplane_integral(self.mesh.faces_centers[:,0]**2)
        return inertia_moment / self.volume

    @property
    def transversal_metacentric_height(self):
        """Returns transversal metacentric height of the body."""
        gb = self.center_of_mass - self.center_of_buoyancy
        return self.transversal_metacentric_radius - gb[2]

    @property
    def longitudinal_metacentric_height(self):
        """Returns longitudinal metacentric height of the body."""
        gb = self.center_of_mass - self.center_of_buoyancy
        return self.longitudinal_metacentric_radius - gb[2]

    def dof_normals(self, dof):
        """Returns dot product of the surface face normals and DOF"""
        return np.sum(self.mesh.faces_normals * dof, axis=1)

    def _infer_rotation_center(self):
        """Hacky way to infer the point around which the rotation dofs are defined.
        (Assuming all three rotation dofs are defined around the same point).
        In the future, should be replaced by something more robust.
        """
        if hasattr(self, "rotation_center"):
            return np.asarray(self.rotation_center)

        else:
            try:
                xc1 = self.dofs["Pitch"][:, 2] + self.mesh.faces_centers[:, 0]
                xc2 = -self.dofs["Yaw"][:, 1] + self.mesh.faces_centers[:, 0]
                yc1 = self.dofs["Yaw"][:, 0] + self.mesh.faces_centers[:, 1]
                yc2 = -self.dofs["Roll"][:, 2] + self.mesh.faces_centers[:, 1]
                zc1 = -self.dofs["Pitch"][:, 0] + self.mesh.faces_centers[:, 2]
                zc2 = self.dofs["Roll"][:, 1] + self.mesh.faces_centers[:, 2]

                # All items should be identical in a given vector
                assert np.isclose(xc1, xc1[0]).all()
                assert np.isclose(yc1, yc1[0]).all()
                assert np.isclose(zc1, zc1[0]).all()

                # Both vector should be identical
                assert np.allclose(xc1, xc2)
                assert np.allclose(yc1, yc2)
                assert np.allclose(zc1, zc2)

                return np.array([xc1[0], yc1[0], zc1[0]])

            except Exception as e:
                raise ValueError(
                        f"Failed to infer the rotation center of {self.name} to compute rigid body hydrostatics.\n"
                        f"Possible fix: add a `rotation_center` attibute to {self.name}.\n"
                        "Note that rigid body hydrostatic methods currently assume that the three rotation dofs have the same rotation center."
                        ) from e

    def each_hydrostatic_stiffness(self, influenced_dof_name, radiating_dof_name, *,
                                         influenced_dof_div=0.0, rho=1000.0, g=9.81):
        r"""
        Return the hydrostatic stiffness for a pair of DOFs.

        :math:`C_{ij} = \rho g\iint_S (\hat{n} \cdot V_j) (w_i + z D_i) dS`

        where :math:`\hat{n}` is surface normal,

        :math:`V_i = u_i \hat{n}_x + v_i \hat{n}_y + w_i \hat{n}_z` is DOF vector and

        :math:`D_i = \nabla \cdot V_i` is the divergence of the DOF.

        Parameters
        ----------
        influenced_dof_name : str
            Name of influenced DOF vector of the FloatingBody
        radiating_dof_name: str
            Name of radiating DOF vector of the FloatingBody
        influenced_dof_div: np.ndarray (Face_count), optional
            Influenced DOF divergence of the FloatingBody, by default 0.0.
        rho: float, optional
            water density, by default 1000.0
        g: float, optional
            Gravity acceleration, by default 9.81

        Returns
        -------
        hs_ij: xarray.variable
            hydrostatic_stiffness of ith DOF and jth DOF.

        Note
        ----
            This function computes the hydrostatic stiffness assuming :math:`D_{i} = 0`.
            If :math:`D_i \neq 0`, input the divergence interpolated to face centers.

            General integral equations are used for the rigid body modes and
            Neumann (1994) method is used for flexible modes.

        References
        ----------
            Newman, John Nicholas. "Wave effects on deformable bodies."Applied ocean
            research" 16.1 (1994): 47-59.
            http://resolver.tudelft.nl/uuid:0adff84c-43c7-43aa-8cd8-d4c44240bed8

        """
        # Newmann (1994) formula is not 'complete' as recovering the rigid body
        # terms is not possible. https://doi.org/10.1115/1.3058702.

        # Alternative is to use the general equation of hydrostatic and
        # restoring coefficient for rigid mdoes and use Neuman equation for elastic
        # modes.

        rigid_dof_names = ("Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw")
        dof_pair = (influenced_dof_name, radiating_dof_name)

        if set(dof_pair).issubset(set(rigid_dof_names)):
            if self.center_of_mass is None:
                raise ValueError(f"Trying to compute rigid-body hydrostatic stiffness for {self.name}, but no center of mass has been defined.\n"
                                 f"Suggested solution: define a `center_of_mass` attribute for the FloatingBody {self.name}.")
            mass = self.disp_mass(rho=rho) if self.mass is None else self.mass
            xc, yc, zc = self._infer_rotation_center()

            if dof_pair == ("Heave", "Heave"):
                norm_hs_stiff = self.waterplane_area
            elif dof_pair in [("Heave", "Roll"), ("Roll", "Heave")]:
                norm_hs_stiff = -self.waterplane_integral(self.mesh.faces_centers[:,1] - yc)
            elif dof_pair in [("Heave", "Pitch"), ("Pitch", "Heave")]:
                norm_hs_stiff = self.waterplane_integral(self.mesh.faces_centers[:,0] - xc)
            elif dof_pair == ("Roll", "Roll"):
                norm_hs_stiff = (
                        -self.waterplane_integral((self.mesh.faces_centers[:,1] - yc)**2)
                        + self.volume*(self.center_of_buoyancy[2] - zc) - mass/rho*(self.center_of_mass[2] - zc)
                )
            elif dof_pair in [("Roll", "Pitch"), ("Pitch", "Roll")]:
                norm_hs_stiff = self.waterplane_integral((self.mesh.faces_centers[:,0] - xc)
                                                          * (self.mesh.faces_centers[:,1] - yc))
            elif dof_pair == ("Roll", "Yaw"):
                norm_hs_stiff = - self.volume*(self.center_of_buoyancy[0] - xc) + mass/rho*(self.center_of_mass[0] - xc)
            elif dof_pair == ("Pitch", "Pitch"):
                norm_hs_stiff = (
                        -self.waterplane_integral((self.mesh.faces_centers[:,0] - xc)**2)
                        + self.volume*(self.center_of_buoyancy[2] - zc) - mass/rho*(self.center_of_mass[2] - zc)
                        )
            elif dof_pair == ("Pitch", "Yaw"):
                norm_hs_stiff = - self.volume*(self.center_of_buoyancy[1] - yc) + mass/rho*(self.center_of_mass[1] - yc)
            else:
                norm_hs_stiff = 0.0
        else:
            if self.mass is not None and np.isclose(self.mass, self.disp_mass(rho), rtol=1e-4):
                raise NotImplementedError(
                        f"Trying to compute the hydrostatic stiffness for dofs {radiating_dof_name} and {influenced_dof_name}"
                        f"of body {self.name}, which is not neutrally buoyant (mass={body.mass}, disp_mass={body.disp_mass(rho)}.\n"
                        f"This case has not been implemented in Capytaine. You need either a single rigid body or a neutrally buoyant body."
                        )

            # Neuman (1994) formula for flexible DOFs
            influenced_dof = np.array(self.dofs[influenced_dof_name])
            radiating_dof = np.array(self.dofs[radiating_dof_name])
            influenced_dof_div_array = np.array(influenced_dof_div)

            radiating_dof_normal = self.dof_normals(radiating_dof)
            z_influenced_dof_div = influenced_dof[:,2] + self.mesh.faces_centers[:,2] * influenced_dof_div_array
            norm_hs_stiff = self.surface_integral( -radiating_dof_normal * z_influenced_dof_div)

        hs_stiff = rho * g * norm_hs_stiff

        return xr.DataArray([[hs_stiff]],
                            dims=['influenced_dof', 'radiating_dof'],
                            coords={'influenced_dof': [influenced_dof_name],
                            'radiating_dof': [radiating_dof_name]},
                            name="hydrostatic_stiffness"
                            )

    def compute_hydrostatic_stiffness(self, *, divergence=None, rho=1000.0, g=9.81):
        r"""
        Compute hydrostatic stiffness matrix for all DOFs of the body.

        :math:`C_{ij} = \rho g\iint_S (\hat{n} \cdot V_j) (w_i + z D_i) dS`

        where :math:`\hat{n}` is surface normal,

        :math:`V_i = u_i \hat{n}_x + v_i \hat{n}_y + w_i \hat{n}_z` is DOF vector and

        :math:`D_i = \nabla \cdot V_i` is the divergence of the DOF.

        Parameters
        ----------
        divergence : dict mapping a dof name to an array of shape (nb_faces) or
                        xarray.DataArray of shape (nb_dofs × nb_faces), optional
            Divergence of the DOFs, by default None
        rho : float, optional
            Water density, by default 1000.0
        g: float, optional
            Gravity acceleration, by default 9.81

        Returns
        -------
        xr.DataArray
            Matrix of hydrostatic stiffness

        Note
        ----
            This function computes the hydrostatic stiffness assuming :math:`D_{i} = 0`.
            If :math:`D_i \neq 0`, input the divergence interpolated to face centers.

            General integral equations are used for the rigid body modes and
            Neumann (1994) method is used for flexible modes.

        References
        ----------
            Newman, John Nicholas. "Wave effects on deformable bodies."Applied ocean
            research" 16.1 (1994): 47-59.
            http://resolver.tudelft.nl/uuid:0adff84c-43c7-43aa-8cd8-d4c44240bed8

        """
        if len(self.dofs) == 0:
            raise AttributeError("Cannot compute hydrostatics stiffness on {} since no dof has been defined.".format(self.name))

        def divergence_dof(influenced_dof):
            if divergence is None:
                return 0.0
            elif isinstance(divergence, dict) and influenced_dof in divergence.keys():
                return divergence[influenced_dof]
            elif isinstance(divergence, xr.DataArray) and influenced_dof in divergence.coords["influenced_dof"]:
                return divergence.sel(influenced_dof=influenced_dof).values
            else:
                LOG.warning("Computing hydrostatic stiffness without the divergence of {}".format(influenced_dof))
                return 0.0

        hs_set =  xr.merge([
            self.each_hydrostatic_stiffness(
                influenced_dof_name, radiating_dof_name,
                influenced_dof_div = divergence_dof(influenced_dof_name),
                rho=rho, g=g
                )
            for radiating_dof_name in self.dofs
            for influenced_dof_name in self.dofs
            ])

        # Reorder dofs
        K = hs_set.hydrostatic_stiffness.sel(influenced_dof=list(self.dofs.keys()), radiating_dof=list(self.dofs.keys()))
        return K

    def compute_rigid_body_inertia(self, *, rho=1000, output_type="body_dofs"):
        """
        Inertia Mass matrix of the body for 6 rigid DOFs.

        Parameters
        ----------
        rho : float, optional
            Density of water, by default 1000.0
        output_type : {"body_dofs", "rigid_dofs", "all_dofs"}
            Type of DOFs for mass mat output, by default "body_dofs".

        Returns
        -------
        xarray.DataArray
            Inertia matrix

        Raises
        ------
        ValueError
            If output_type is not in {"body_dofs", "rigid_dofs", "all_dofs"}.
        """
        if self.center_of_mass is None:
            raise ValueError(f"Trying to compute rigid-body inertia matrix for {self.name}, but no center of mass has been defined.\n"
                             f"Suggested solution: define a `center_of_mass` attribute for the FloatingBody {self.name}.")

        rc = self._infer_rotation_center()
        fcs = (self.mesh.faces_centers - rc).T
        combinations = np.array([fcs[0]**2, fcs[1]**2, fcs[2]**2, fcs[0]*fcs[1],
                                 fcs[1]*fcs[2], fcs[2]*fcs[0]])
        integrals = np.array([
            [np.sum(normal_i * fcs[axis] * combination * self.mesh.faces_areas)
            for combination in combinations]
            for axis, normal_i in enumerate(self.mesh.faces_normals.T)])


        inertias = np.array([
            (integrals[0,1]   + integrals[0,2]   + integrals[1,1]/3
             + integrals[1,2]   + integrals[2,1] + integrals[2,2]/3)/3,
            (integrals[0,0]/3 + integrals[0,2]   + integrals[1,0]
             + integrals[1,2]   + integrals[2,0] + integrals[2,2]/3)/3,
            (integrals[0,0]/3 + integrals[0,1]   + integrals[1,0]
             + integrals[1,1]/3 + integrals[2,0] + integrals[2,1]  )/3,
            integrals[2,3],
            integrals[0,4],
            integrals[1,5]
        ])

        cog = self.center_of_mass - rc
        volume = self.volume
        volumic_inertia_matrix = np.array([
            [ volume        , 0              , 0               ,
              0             , volume*cog[2]  , -volume*cog[1]  ],
            [ 0             , volume         , 0               ,
             -volume*cog[2] , 0              , volume*cog[0]   ],
            [ 0             , 0              , volume          ,
              volume*cog[1] , -volume*cog[0] , 0 ]             ,
            [ 0             , -volume*cog[2] , volume*cog[1]   ,
              inertias[0]   , -inertias[3]   , -inertias[5]    ],
            [ volume*cog[2] , 0              , -volume*cog[0]  ,
             -inertias[3]   , inertias[1]    , -inertias[4]    ],
            [-volume*cog[1] , volume*cog[0]  , 0               ,
             -inertias[5]   , -inertias[4]   , inertias[2]     ],
        ])

        density = rho if self.mass is None else self.mass/volume
        inertia_matrix = density * volumic_inertia_matrix

        # Rigid DOFs
        rigid_dof_names = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]
        rigid_inertia_matrix_xr = xr.DataArray(data=np.asarray(inertia_matrix),
                            dims=['influenced_dof', 'radiating_dof'],
                            coords={'influenced_dof': rigid_dof_names,
                                    'radiating_dof': rigid_dof_names},
                            name="inertia_matrix")

        # Body DOFs (Default as np.nan)
        body_dof_names = list(self.dofs)
        body_dof_count = len(body_dof_names)
        other_dofs_inertia_matrix_xr = xr.DataArray(np.nan * np.zeros([body_dof_count, body_dof_count]),
                                    dims=['influenced_dof', 'radiating_dof'],
                                    coords={'influenced_dof': body_dof_names,
                                            'radiating_dof': body_dof_names},
                                    name="inertia_matrix")

        total_mass_xr = xr.merge([rigid_inertia_matrix_xr, other_dofs_inertia_matrix_xr], compat="override").inertia_matrix

        non_rigid_dofs = set(body_dof_names) - set(rigid_dof_names)

        if output_type == "body_dofs":
            if len(non_rigid_dofs) > 0:
                LOG.warning(f"Non-rigid dofs: {non_rigid_dofs} are detected and \
respective inertia coefficients are assigned as NaN.")

            inertia_matrix_xr = total_mass_xr.sel(influenced_dof=body_dof_names,
                                                  radiating_dof=body_dof_names)
        elif output_type == "rigid_dofs":
            inertia_matrix_xr = total_mass_xr.sel(influenced_dof=rigid_dof_names,
                                                  radiating_dof=rigid_dof_names)
        elif output_type == "all_dofs":
            if len(non_rigid_dofs) > 0:
                LOG.warning("Non-rigid dofs: {non_rigid_dofs} are detected and \
respective inertia coefficients are assigned as NaN.")

            inertia_matrix_xr = total_mass_xr
        else:
            raise ValueError(f"output_type should be either 'body_dofs', \
'all_dofs' or 'rigid_dofs'. Given output_type = '{output_type}'.")

        return inertia_matrix_xr


    def compute_hydrostatics(self, *, rho=1000.0, g=9.81, divergence=None):
        """Compute hydrostatics of the FloatingBody.

        Parameters
        ----------
        rho : float, optional
            Density of Water. The default is 1000.
        g: float, optional
            Gravity acceleration. The default is 9.81.
        divergence : np.ndarray, optional
            Divergence of the DOFs.

        Returns
        -------
        hydrostatics : dict
            All hydrostatics values of the FloatingBody.
        """
        if self.center_of_mass is None:
            raise ValueError(f"Trying to compute hydrostatics for {self.name}, but no center of mass has been defined.\n"
                             f"Suggested solution: define a `center_of_mass` attribute for the FloatingBody {self.name}.")

        self.keep_immersed_part()

        full_mesh_vertices = self.full_body.mesh.vertices
        coord_max = full_mesh_vertices.max(axis=0)
        coord_min = full_mesh_vertices.min(axis=0)
        full_length, full_breadth, depth = full_mesh_vertices.max(axis=0) - full_mesh_vertices.min(axis=0)

        vertices = self.mesh.vertices
        sub_length, sub_breadth, _ = vertices.max(axis=0) - vertices.min(axis=0)

        if abs(self.waterplane_area) > 1e-10:
            water_plane_idx = np.isclose(vertices[:,2], 0.0)
            water_plane = vertices[water_plane_idx][:,:-1]
            wl_length, wl_breadth = water_plane.max(axis=0) - water_plane.min(axis=0)
        else:
            wl_length, wl_breadth = 0.0, 0.0

        hydrostatics = {}
        hydrostatics["g"] = g
        hydrostatics["rho"] = rho
        hydrostatics["center_of_mass"] = self.center_of_mass

        hydrostatics["wet_surface_area"] = self.wet_surface_area
        hydrostatics["disp_volumes"] = self.volumes
        hydrostatics["disp_volume"] = self.volume
        hydrostatics["disp_mass"] = self.disp_mass(rho=rho)
        hydrostatics["center_of_buoyancy"] = self.center_of_buoyancy
        hydrostatics["waterplane_center"] = np.append(self.waterplane_center, 0.0)
        hydrostatics["waterplane_area"] = self.waterplane_area
        hydrostatics["transversal_metacentric_radius"] = self.transversal_metacentric_radius
        hydrostatics["longitudinal_metacentric_radius"] = self.longitudinal_metacentric_radius
        hydrostatics["transversal_metacentric_height"] = self.transversal_metacentric_height
        hydrostatics["longitudinal_metacentric_height"] = self.longitudinal_metacentric_height
        hydrostatics["hydrostatic_stiffness"] = self.compute_hydrostatic_stiffness(
            divergence=divergence, rho=rho, g=g)

        hydrostatics["length_overall"] = full_length
        hydrostatics["breadth_overall"] = full_breadth
        hydrostatics["depth"] = depth
        hydrostatics["draught"] = np.abs(coord_min[2])
        hydrostatics["length_at_waterline"] = wl_length
        hydrostatics["breadth_at_waterline"] = wl_breadth
        hydrostatics["length_overall_submerged"] = sub_length
        hydrostatics["breadth_overall_submerged"] = sub_breadth
        hydrostatics["inertia_matrix"] = self.compute_rigid_body_inertia(rho=rho)

        return hydrostatics


    ###################
    # Transformations #
    ###################

    def __add__(self, body_to_add: 'FloatingBody') -> 'FloatingBody':
        return self.join_bodies(body_to_add)

    def join_bodies(*bodies, name=None) -> 'FloatingBody':
        if name is None:
            name = "+".join(body.name for body in bodies)
        meshes = CollectionOfMeshes([body.mesh for body in bodies], name=f"{name}_mesh")
        dofs = FloatingBody.combine_dofs(bodies)

        if all(body.mass is not None for body in bodies):
            new_mass = sum(body.mass is not None for body in bodies)
        else:
            new_mass = None

        if (all(body.mass is not None for body in bodies)
                and all(body.center_of_mass is not None for body in bodies)):
            new_cog = sum(body.mass*body.center_of_mass is not None for body in bodies)/new_mass
        else:
            new_cog = None

        return FloatingBody(
            mesh=meshes, dofs=dofs, mass=new_mass, center_of_mass=new_cog, name=name
            )

    @staticmethod
    def combine_dofs(bodies) -> dict:
        """Combine the degrees of freedom of several bodies."""
        dofs = {}
        cum_nb_faces = accumulate(chain([0], (body.mesh.nb_faces for body in bodies)))
        total_nb_faces = sum(body.mesh.nb_faces for body in bodies)
        for body, nbf in zip(bodies, cum_nb_faces):
            # nbf is the cumulative number of faces of the previous subbodies,
            # that is the offset of the indices of the faces of the current body.
            for name, dof in body.dofs.items():
                new_dof = np.zeros((total_nb_faces, 3))
                new_dof[nbf:nbf+len(dof), :] = dof
                if '__' not in name:
                    new_dof_name = '__'.join([body.name, name])
                else:
                    # The body is probably a combination of bodies already.
                    # So for the associativity of the + operation,
                    # it is better to keep the same name.
                    new_dof_name = name
                dofs[new_dof_name] = new_dof
        return dofs

    def copy(self, name=None) -> 'FloatingBody':
        """Return a deep copy of the body.

        Parameters
        ----------
        name : str, optional
            a name for the new copy
        """
        new_body = copy.deepcopy(self)
        if name is None:
            new_body.name = f"copy_of_{self.name}"
            LOG.debug(f"Copy {self.name}.")
        else:
            new_body.name = name
            LOG.debug(f"Copy {self.name} under the name {name}.")
        return new_body

    def assemble_regular_array(self, distance, nb_bodies):
        """Create an regular array of identical bodies.

        Parameters
        ----------
        distance : float
            Center-to-center distance between objects in the array
        nb_bodies : couple of ints
            Number of objects in the x and y directions.

        Returns
        -------
        FloatingBody
        """
        bodies = (self.translated((i*distance, j*distance, 0), name=f"{i}_{j}") for j in range(nb_bodies[1]) for i in range(nb_bodies[0]))
        array = FloatingBody.join_bodies(*bodies)
        array.mesh = build_regular_array_of_meshes(self.mesh, distance, nb_bodies)
        array.name = f"array_of_{self.name}"
        return array

    def assemble_arbitrary_array(self, locations:np.ndarray):

        if not isinstance(locations, np.ndarray):
            raise TypeError('locations must be of type np.ndarray')
        assert locations.shape[1] == 2, 'locations must be of shape nx2, received {:}'.format(locations.shape)
        n = locations.shape[0]

        fb_list = []
        for idx, li in enumerate(locations):
            fb1 = self.copy()
            fb1.translate(np.append(li,0))
            fb1.name = 'arbitrary_array_body{:02d}'.format(idx)
            fb_list.append(fb1)

        arbitrary_array = fb_list[0].join_bodies(*fb_list[1:])

        return arbitrary_array

    def extract_faces(self, id_faces_to_extract, return_index=False):
        """Create a new FloatingBody by extracting some faces from the mesh.
        The dofs evolve accordingly.
        """
        if isinstance(self.mesh, CollectionOfMeshes):
            raise NotImplementedError  # TODO

        if return_index:
            new_mesh, id_v = Mesh.extract_faces(self.mesh, id_faces_to_extract, return_index)
        else:
            new_mesh = Mesh.extract_faces(self.mesh, id_faces_to_extract, return_index)
        new_body = FloatingBody(new_mesh)
        LOG.info(f"Extract floating body from {self.name}.")

        new_body.dofs = {}
        for name, dof in self.dofs.items():
            new_body.dofs[name] = dof[id_faces_to_extract, :]

        if return_index:
            return new_body, id_v
        else:
            return new_body

    def sliced_by_plane(self, plane):
        return FloatingBody(mesh=self.mesh.sliced_by_plane(plane), dofs=self.dofs, name=self.name)

    def minced(self, nb_slices=(8, 8, 4)):
        """Experimental method decomposing the mesh as a hierarchical structure.

        Parameters
        ----------
        nb_slices: Tuple[int, int, int]
            The number of slices in each of the x, y and z directions.
            Only powers of 2 are supported at the moment.

        Returns
        -------
        FloatingBody
        """
        minced_body = self.copy()

        # Extreme points of the mesh in each directions.
        x_min, x_max, y_min, y_max, z_min, z_max = self.mesh.axis_aligned_bbox
        sizes = [(x_min, x_max), (y_min, y_max), (z_min, z_max)]

        directions = [np.array(d) for d in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]]

        def _slice_positions_at_depth(i):
            """Helper function.

            Returns a list of floats as follows:
            i=1 -> [1/2]
            i=2 -> [1/4, 3/4]
            i=3 -> [1/8, 3/8, 5/8, 7/8]
                   ...
            """
            denominator = 2**i
            return [numerator/denominator for numerator in range(1, denominator, 2)]

        # GENERATE ALL THE PLANES THAT WILL BE USED TO MINCE THE MESH
        planes = []
        for direction, nb_slices_in_dir, (min_coord, max_coord) in zip(directions, nb_slices, sizes):
            planes_in_dir = []

            depth_of_treelike_structure = int(np.log2(nb_slices_in_dir))
            for i_depth in range(1, depth_of_treelike_structure+1):
                planes_in_dir_at_depth = []
                for relative_position in _slice_positions_at_depth(i_depth):
                    slice_position = (min_coord + relative_position*(max_coord-min_coord))*direction
                    plane = Plane(normal=direction, point=slice_position)
                    planes_in_dir_at_depth.append(plane)
                planes_in_dir.append(planes_in_dir_at_depth)
            planes.append(planes_in_dir)

        # SLICE THE MESH
        intermingled_x_y_z = chain.from_iterable(zip_longest(*planes))
        for planes in intermingled_x_y_z:
            if planes is not None:
                for plane in planes:
                    minced_body = minced_body.sliced_by_plane(plane)
        return minced_body

    @inplace_transformation
    def mirror(self, plane):
        self.mesh.mirror(plane)
        for dof in self.dofs:
            self.dofs[dof] -= 2 * np.outer(np.dot(self.dofs[dof], plane.normal), plane.normal)
        for point_attr in ('geometric_center', 'rotation_center', 'center_of_mass'):
            if point_attr in self.__dict__ and self.__dict__[point_attr] is not None:
                self.__dict__[point_attr] -= 2 * (np.dot(self.__dict__[point_attr], plane.normal) - plane.c) * plane.normal
        return self

    @inplace_transformation
    def translate(self, *args):
        self.mesh.translate(*args)
        for point_attr in ('geometric_center', 'rotation_center', 'center_of_mass'):
            if point_attr in self.__dict__ and self.__dict__[point_attr] is not None:
                self.__dict__[point_attr] += args[0]
        return self

    @inplace_transformation
    def rotate(self, axis, angle):
        self.mesh.rotate(axis, angle)
        for point_attr in ('geometric_center', 'rotation_center', 'center_of_mass'):
            if point_attr in self.__dict__ and self.__dict__[point_attr] is not None:
                self.__dict__[point_attr] = axis.rotate_points([self.__dict__[point_attr]], angle)
        for dof in self.dofs:
            self.dofs[dof] = axis.rotate_vectors(self.dofs[dof], angle)
        return self

    @inplace_transformation
    def clip(self, plane):
        # Keep of copy of the full mesh
        if self.full_body is None:
            self.full_body = self.copy()

        # Clip mesh
        LOG.info(f"Clipping {self.name} with respect to {plane}")
        self.mesh.clip(plane)

        # Clip dofs
        ids = self.mesh._clipping_data['faces_ids']
        for dof in self.dofs:
            if len(ids) > 0:
                self.dofs[dof] = np.array(self.dofs[dof])[ids]
            else:
                self.dofs[dof] = np.empty((0, 3))
        return self

    def clipped(self, plane, **kwargs):
        # Same API as for the other transformations
        return self.clip(plane, inplace=False, **kwargs)

    @inplace_transformation
    def keep_immersed_part(self, free_surface=0.0, sea_bottom=-np.infty):
        """Remove the parts of the mesh above the sea bottom and below the free surface."""
        self.clip(Plane(normal=(0, 0, 1), point=(0, 0, free_surface)))
        if sea_bottom > -np.infty:
            self.clip(Plane(normal=(0, 0, -1), point=(0, 0, sea_bottom)))
        return self

    #############
    #  Display  #
    #############

    def __str__(self):
        return self.name

    def __repr__(self):
        return (f"{self.__class__.__name__}(mesh={self.mesh.name}, "
                f"dofs={{{', '.join(self.dofs.keys())}}}, name={self.name})")

    def show(self, **kwargs):
        from capytaine.ui.vtk.body_viewer import FloatingBodyViewer
        viewer = FloatingBodyViewer()
        viewer.add_body(self, **kwargs)
        viewer.show()
        viewer.finalize()

    def show_matplotlib(self, *args, **kwargs):
        return self.mesh.show_matplotlib(*args, **kwargs)

    def animate(self, motion, *args, **kwargs):
        """Display a motion as a 3D animation.

        Parameters
        ==========
        motion: dict or pd.Series or str
            A dict or series mapping the name of the dofs to its amplitude.
            If a single string is passed, it is assumed to be the name of a dof
            and this dof with a unit amplitude will be displayed.
        """
        from capytaine.ui.vtk.animation import Animation
        if isinstance(motion, str):
            motion = {motion: 1.0}
        elif isinstance(motion, xr.DataArray):
            motion = {k: motion.sel(radiating_dof=k).data for k in motion.coords["radiating_dof"].data}

        if any(dof not in self.dofs for dof in motion):
            missing_dofs = set(motion.keys()) - set(self.dofs.keys())
            raise ValueError(f"Trying to animate the body {self.name} using dof(s) {missing_dofs}, but no dof of this name is defined for {self.name}.")

        animation = Animation(*args, **kwargs)
        animation._add_actor(self.mesh.merged(), faces_motion=sum(motion[dof_name] * dof for dof_name, dof in self.dofs.items() if dof_name in motion))
        return animation

