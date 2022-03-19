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
    name : str, optional
        a name for the body.
        If none is given, the one of the mesh is used.
    """

    def __init__(self, mesh=None, dofs=None, name=None):
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
                    if hasattr(self, point_attr):
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

        if hasattr(self, 'mass'):
            self.mass = self.mass.sel(radiating_dof=dofs, influenced_dof=dofs)
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

    def get_surface_integral(self, data, **kwargs):
        """Returns integral of given data along wet surface area."""
        return np.sum(data * self.mesh.faces_areas, **kwargs)

    def get_waterplane_integral(self, data, **kwargs):
        """Returns integral of given data along water plane area."""
        return self.get_surface_integral(self.mesh.faces_normals[:,2] * data, **kwargs)

    def get_wet_surface_area(self):
        """Returns wet surface area."""
        return self.get_surface_integral(1)

    def get_volumes(self):
        """Returns volumes using x, y, z components of the FloatingBody."""
        norm_coord = self.mesh.faces_normals * self.mesh.faces_centers
        return self.get_surface_integral(norm_coord.T, axis=1)

    def get_volume(self):
        """Returns volume of the FloatingBody."""
        volumes = self.get_volumes()
        return np.mean(volumes)

    def get_mass(self, density=1000):
        """Returns mass of the FloatingBody."""
        return density * self.get_volume()

    def get_buoyancy_center(self):
        """Returns center of buoyancy of the FloatingBody."""
        volume = self.get_volume()
        coords_sq_norm = self.mesh.faces_normals * self.mesh.faces_centers**2
        return self.get_surface_integral(coords_sq_norm.T, axis=1) / (2*volume)

    def get_waterplane_area(self, free_surface=0.0):
        """Returns water plane area of the FloatingBody."""
        # if self.mesh.vertices[:,2].max() < free_surface:
        #     waterplane_area = 0.0
        # else:
        waterplane_area = -self.get_waterplane_integral(1)
        return waterplane_area

    def get_waterplane_center(self, free_surface=0.0):
        """Returns water plane center of the FloatingBody.
        
        Note: Returns None if the FloatingBody is full submerged. 
        """
        # if self.mesh.vertices[:,2].max() < free_surface:
        #     waterplane_center = None
        # else:
        waterplane_area = self.get_waterplane_area()
        waterplane_center = -self.get_waterplane_integral(
            self.mesh.faces_centers.T, axis=1) / waterplane_area
        waterplane_center[2] = free_surface

        return waterplane_center

    def get_bmt(self):
        """Returns transversal metacentric radius of the body."""
        volume = self.get_volume()
        inertia_moment = -self.get_waterplane_integral(self.mesh.faces_centers[:,1]**2)
        bmt = inertia_moment / volume
        return bmt

    def get_bml(self):
        """Returns longitudinal metacentric radius of the body."""
        volume = self.get_volume()
        inertia_moment = -self.get_waterplane_integral(self.mesh.faces_centers[:,0]**2)
        bmt = inertia_moment / volume
        return bmt

    def get_gmt(self, cog=np.zeros(3)):
        """Returns transversal metacentric Height of the body."""
        gb = cog - self.get_buoyancy_center()
        return self.get_bmt() - gb[2]

    def get_gml(self, cog=np.zeros(3)):
        """Returns longitudinal metacentric Height of the body."""
        gb = cog - self.get_buoyancy_center()
        return self.get_bml() - gb[2]

    def get_dof_normals(self, dof):
        """Returns dot product of the surface face normals and DOF"""
        return np.sum(self.mesh.faces_normals * dof, axis=1)

    def get_hydrostatic_stiffnessij(self, dof_i_name, dof_j_name, divergence_i=0.0, 
                                      density=1000.0, gravity=9.80665, cog=[0,0,0]):
        r"""
        Return the hydrostatic stiffness for a set of DOFs.
        
        Parameters
        ----------
        dof_i_name : str
            Name of ith DOF vector of the FloatingBody
        dof_j_name: str
            Name of ith DOF vector of the FloatingBody
        divergence_i: np.ndarray (Face_count), optional
            ith DOF divergence of the FloatingBody, by default None. 
        density: float, optional
            water density, by default 1000. 
        gravity : float, optional
            Gravity, by default 9.80665
            
        Returns
        -------
        hydrostatic_stiffnessij: float
            hydrostatic_stiffness of ith DOF and jth DOF.
            
        Equation
        ------
        :math:`C_{ij} = \rho g\iint_S (\hat{n} \cdot V_j) (w_i + z D_i) dS`
        
        where :math:`\hat{n}` is surface normal, 
        
        :math:`V_i = u_i \hat{n}_x + v_i \hat{n}_y + w_i \hat{n}_z` is DOF vector and
        
        :math:`D_i = \nabla \cdot V_i` is the divergence of the DOF.
        
        NOTE
        -----
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
        
        dof_pair = (dof_i_name, dof_j_name)
        
        if set(dof_pair).issubset(set(rigid_dof_names)):
            if dof_pair == ("Heave", "Heave"):
                norm_hs_ij = self.get_waterplane_area()
            elif dof_pair in [("Heave", "Roll"), ("Roll", "Heave")]:
                norm_hs_ij = -self.get_waterplane_integral(self.mesh.faces_centers[:,1])
            elif dof_pair in [("Heave", "Pitch"), ("Pitch", "Heave")]:
                norm_hs_ij = self.get_waterplane_integral(self.mesh.faces_centers[:,0])
            elif dof_pair == ("Roll", "Roll"):
                norm_hs_ij = self.get_volume() * self.get_gmt()
            elif dof_pair in [("Roll", "Pitch"), ("Pitch", "Roll")]:
                norm_hs_ij = self.get_waterplane_integral(self.mesh.faces_centers[:,0]
                                                          * self.mesh.faces_centers[:,1])
            elif dof_pair == ("Roll", "Yaw"):
                norm_hs_ij = self.get_volume() * (-self.get_buoyancy_center()[0] + cog[0])
            elif dof_pair == ("Pitch", "Pitch"):
                norm_hs_ij = self.get_volume() * self.get_gml()
            elif dof_pair == ("Pitch", "Yaw"):
                norm_hs_ij = self.get_volume() * (-self.get_buoyancy_center()[1] + cog[1])
            else:
                norm_hs_ij = 0.0
        else:
            # Neuman (1994) formula for flexible DOFs
            dof_i_array = np.array(self.dofs[dof_i_name])
            dof_j_array = np.array(self.dofs[dof_j_name])
            divergence_i_array = np.array(divergence_i)
            
            dof_j_normal = self.get_dof_normals(dof_j_array)
            dof_z_div_i = dof_i_array[:,2] + self.mesh.faces_centers[:,2] * divergence_i_array
            dof_z_div_i_norm = -dof_j_normal * dof_z_div_i
            norm_hs_ij = self.get_surface_integral(dof_z_div_i_norm)
        
        hs_ij = density * gravity * norm_hs_ij
        return hs_ij

    def get_hydrostatic_stiffness(self, divergence=None, 
                                  density=1000.0, gravity=9.80665, cog=[0,0,0]):
        r"""
        Compute hydrostatic stiffness matrix for all DOFs of the body.
        
        Parameters
        ----------
        divergence : np.ndarray (DOF_count X Face_count), optional
            Divergence of the DOFs, by default None
        density : float, optional
            Water density, by default 1000.0
        gravity : float, optional
            Gravity, by default 9.80665
            
        Returns
        -------
        hydrostatic_stiffness: np.ndarray
            Matrix of hydrostatic stiffness
        
        Equation
        ------
        :math:`C_{ij} = \rho g \iint_S (\hat{n} \cdot V_j) (w_i + z D_i) dS`
        
        where :math:`\hat{n}` is surface normal, 
        
        :math:`V_i = u_i \hat{n}_x + v_i \hat{n}_y + w_i \hat{n}_z` is DOF vector and
        
        :math:`D_i = \nabla \cdot V_i` is the divergence of the DOF.
        
        NOTE: 
        -----
        This function computes the hydrostatic stiffness assuming :math:`D_{i} = 0`. 
        If :math:`D_i \neq 0`, input the divergence interpolated to face centers. 
        
        References
        ----------
        Newman, John Nicholas. "Wave effects on deformable bodies."Applied ocean 
        research" 16.1 (1994): 47-59.
        http://resolver.tudelft.nl/uuid:0adff84c-43c7-43aa-8cd8-d4c44240bed8

        """

        hydrostatic_stiffness =  np.array([
            [self.get_hydrostatic_stiffnessij(
                dof_i_name, dof_j_name, 
                divergence_i = 0.0 if divergence is None else divergence[i],
                density=density, gravity=gravity, cog=cog
                )
            for dof_j_name in self.dofs]
            for i, dof_i_name in enumerate(self.dofs)
            ])
        
        return hydrostatic_stiffness

    def get_rigid_dof_mass(self, cog=np.zeros(3), density=1000):
        """Interia Mass matrix of the body for 6 rigid DOFs."""
        fcs = (self.mesh.faces_centers).T
        combinations = np.array([fcs[0]**2, fcs[1]**2, fcs[2]**2, fcs[0]*fcs[1], 
                                 fcs[1]*fcs[2], fcs[2]*fcs[0]])
        integrals = np.array([
            [np.sum(normal_i * fcs[axis] * combination * self.mesh.faces_areas)
            for combination in combinations]
            for axis, normal_i in enumerate(self.mesh.faces_normals.T)])

        interias = density * np.array([
            (integrals[0,1]   + integrals[0,2]   + integrals[1,1]/3 \
             + integrals[1,2]   + integrals[2,1] + integrals[2,2]/3)/3,
            (integrals[0,0]/3 + integrals[0,2]   + integrals[1,0]   \
             + integrals[1,2]   + integrals[2,0] + integrals[2,2]/3)/3,
            (integrals[0,0]/3 + integrals[0,1]   + integrals[1,0]   \
             + integrals[1,1]/3 + integrals[2,0] + integrals[2,1]  )/3,
            integrals[2,3],
            integrals[0,4],
            integrals[1,5]
        ])

        mass = self.get_mass()
        mass_mat = np.array([
            [ mass       ,  0          ,  0           ,  
              0          ,  mass*cog[2], -mass*cog[1]],
            [ 0          ,  mass       ,  0           ,
             -mass*cog[2],  0          ,  mass*cog[0]],
            [ 0          ,  0          ,  mass        ,
              mass*cog[1], -mass*cog[0],  0          ],
            [ 0          , -mass*cog[2],  mass*cog[1] ,  
              interias[0], -interias[3], -interias[5]],
            [ mass*cog[2],  0          , -mass*cog[0] , 
             -interias[3],  interias[1], -interias[4]],
            [-mass*cog[1],  mass*cog[0],  0           ,
             -interias[5], -interias[4], interias[2]] ,
        ])

        return mass_mat


    def compute_hydrostatics(self, cog=np.zeros(3), density=1000, gravity=9.80665, 
                             free_surface=0.0, divergence=None):
        """Compute hydrostatics of the FloatingBody.
        
        Parameters
        ----------
        cog : np.ndarray, optional
            Center of gravity. The default is np.zeros(3).
        density : float, optional
            Density of Water. The default is 1000.
        gravity : float, optional
            Gravity. The default is 9.80665.
        free_surface : float, optional
            z coordinate of the free surface. The default is 0.0.
        divergence : np.ndarray, optional
            Divergence of the DOFs.
            
        Returns
        -------
        hydrostatics : dict
            All hydrostatics values of the FloatingBody.
        """
        vertices = self.mesh.vertices
        coord_max = vertices.max(axis=0)
        coord_min = vertices.min(axis=0)

        full_length, full_breadth, depth = vertices.max(axis=0) - vertices.min(axis=0)
        
        # Check whether the FloatingBody is below free surface. 
        
        if coord_max[2] >= free_surface:
            self.keep_immersed_part(free_surface=free_surface)
            
        self_coords = self.mesh.vertices
        water_plane_idx = np.isclose(self_coords[:,2], 0.0)
        water_plane = self_coords[water_plane_idx][:,:-1]
        wl_length, wl_breadth = water_plane.max(axis=0) - water_plane.min(axis=0)

        sub_length, sub_breadth, _ = self_coords.max(axis=0) \
                                    - self_coords.min(axis=0)

        hydrostatics = {}
        hydrostatics["grav"] = gravity
        hydrostatics["rho_water"] = density
        hydrostatics["cog"] = cog
        
        hydrostatics["wet_surface_area"] = self.get_wet_surface_area()
        hydrostatics["disp_volumes"] = self.get_volumes()
        hydrostatics["disp_volume"] = self.get_volume()
        hydrostatics["disp_mass"] = self.get_mass(density=density)
        hydrostatics["buoyancy_center"] = self.get_buoyancy_center()
        hydrostatics["waterplane_center"] = self.get_waterplane_center(
            free_surface=free_surface)
        hydrostatics["waterplane_area"] = self.get_waterplane_area(
            free_surface=free_surface)
        hydrostatics["transversal_metacentric_radius"] = self.get_bmt()
        hydrostatics["longitudinal_metacentric_radius"] = self.get_bml()
        hydrostatics["transversal_metacentric_height"] = self.get_gmt(cog=cog)
        hydrostatics["longitudinal_metacentric_height"] = self.get_gml(cog=cog)
        hydrostatics["stiffness_matrix"] = self.get_hydrostatic_stiffness(
            divergence=divergence, density=density, gravity=gravity, cog=cog)

        hydrostatics["length_overall"] = full_length
        hydrostatics["breadth_overall"] = full_breadth
        hydrostatics["depth"] = depth
        hydrostatics["draught"] = np.abs(coord_min[2])
        hydrostatics["length_at_waterline"] = wl_length
        hydrostatics["breadth_at_waterline"] = wl_breadth
        hydrostatics["length_overall_submerged"] = sub_length
        hydrostatics["breadth_overall_submerged"] = sub_breadth
        hydrostatics["inertia_matrix"] = self.get_rigid_dof_mass(
            cog=cog, density=density)

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
        return FloatingBody(mesh=meshes, dofs=dofs, name=name)

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
        array_mesh = build_regular_array_of_meshes(self.mesh, distance, nb_bodies)
        total_nb_faces = array_mesh.nb_faces
        array_dofs = {}
        for dof_name, dof in self.dofs.items():
            for i, j in product(range(nb_bodies[0]), range(nb_bodies[1])):
                shift_nb_faces = (j*nb_bodies[0] + i) * self.mesh.nb_faces
                new_dof = np.zeros((total_nb_faces, 3))
                new_dof[shift_nb_faces:shift_nb_faces+len(dof), :] = dof
                array_dofs[f'{i}_{j}__{dof_name}'] = new_dof
        return FloatingBody(mesh=array_mesh, dofs=array_dofs, name=f"array_of_{self.name}")

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
            if point_attr in self.__dict__:
                self.__dict__[point_attr] -= 2 * (np.dot(self.__dict__[point_attr], plane.normal) - plane.c) * plane.normal
        return self

    @inplace_transformation
    def translate(self, *args):
        self.mesh.translate(*args)
        for point_attr in ('geometric_center', 'rotation_center', 'center_of_mass'):
            if point_attr in self.__dict__:
                self.__dict__[point_attr] += args[0]
        return self

    @inplace_transformation
    def rotate(self, axis, angle):
        matrix = axis.rotation_matrix(angle)
        self.mesh.rotate(axis, angle)
        for point_attr in ('geometric_center', 'rotation_center', 'center_of_mass'):
            if point_attr in self.__dict__:
                self.__dict__[point_attr] = matrix @ self.__dict__[point_attr]
        for dof in self.dofs:
            self.dofs[dof] = (matrix @ self.dofs[dof].T).T
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
                self.dofs[dof] = self.dofs[dof][ids]
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
        animation = Animation(*args, **kwargs)
        animation._add_actor(self.mesh.merged(), faces_motion=sum(motion[dof_name] * dof for dof_name, dof in self.dofs.items() if dof_name in motion))
        return animation

