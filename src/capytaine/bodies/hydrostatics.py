# Copyright (C) 2025 Capytaine developers
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging

import numpy as np
import xarray as xr
from abc import ABC

from capytaine.bodies.dofs import TranslationDof, RotationDof

LOG = logging.getLogger(__name__)

class _HydrostaticsMixin(ABC):
    # This class is not meant to be instantiated but only to be inherited by other classes to give them more methods.

    @property
    def wet_surface_area(self):
        """Returns wet surface area."""
        return self.mesh.wet_surface_area

    @property
    def volumes(self):
        """Returns volumes using x, y, z components of the FloatingBody."""
        return self.mesh.volumes

    @property
    def volume(self):
        """Returns volume of the FloatingBody."""
        return self.mesh.volume

    @property
    def disp_volume(self):
        return self.mesh.disp_volume

    def disp_mass(self, *, rho=1000.0):
        return self.mesh.disp_mass(rho=rho)

    @property
    def center_of_buoyancy(self):
        """Returns center of buoyancy of the FloatingBody."""
        return self.mesh.center_of_buoyancy

    @property
    def waterplane_area(self):
        """Returns water plane area of the FloatingBody."""
        return self.mesh.waterplane_area

    @property
    def waterplane_center(self):
        """Returns water plane center of the FloatingBody.

        Note: Returns None if the FloatingBody is full submerged.
        """
        return self.mesh.waterplane_center

    @property
    def transversal_metacentric_radius(self):
        """Returns transversal metacentric radius of the mesh."""
        immersed_mesh = self.mesh.immersed_part()
        inertia_moment = -immersed_mesh.waterplane_integral(immersed_mesh.faces_centers[:,1]**2)
        return inertia_moment / self.disp_volume

    @property
    def longitudinal_metacentric_radius(self):
        """Returns longitudinal metacentric radius of the mesh."""
        immersed_mesh = self.mesh.immersed_part()
        inertia_moment = -immersed_mesh.waterplane_integral(immersed_mesh.faces_centers[:,0]**2)
        return inertia_moment / self.disp_volume

    @property
    def transversal_metacentric_height(self):
        """Returns transversal metacentric height of the mesh."""
        gb = self.center_of_mass - self.center_of_buoyancy
        return self.transversal_metacentric_radius - gb[2]

    @property
    def longitudinal_metacentric_height(self):
        """Returns longitudinal metacentric height of the mesh."""
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
                        f"Possible fix: add a `rotation_center` attribute to {self.name}.\n"
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
        # Newman (1994) formula is not 'complete' as recovering the rigid body
        # terms is not possible. https://doi.org/10.1115/1.3058702.

        # Alternative is to use the general equation of hydrostatic and
        # restoring coefficient for rigid modes and use Newman equation for elastic
        # modes.

        rigid_dof_names = ("Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw")
        dof_pair = (influenced_dof_name, radiating_dof_name)

        if set(dof_pair).issubset(set(rigid_dof_names)):
            if self.center_of_mass is None:
                raise ValueError(f"Trying to compute rigid-body hydrostatic stiffness for {self.name}, but no center of mass has been defined.\n"
                                 f"Suggested solution: define a `center_of_mass` attribute for the FloatingBody {self.name}.")
            mass = self.disp_mass(rho=rho) if self.mass is None else self.mass
            xc, yc, zc = self._infer_rotation_center()
            immersed_mesh = self.mesh.immersed_part()

            if dof_pair == ("Heave", "Heave"):
                norm_hs_stiff = immersed_mesh.waterplane_area
            elif dof_pair in [("Heave", "Roll"), ("Roll", "Heave")]:
                norm_hs_stiff = -immersed_mesh.waterplane_integral(immersed_mesh.faces_centers[:,1] - yc)
            elif dof_pair in [("Heave", "Pitch"), ("Pitch", "Heave")]:
                norm_hs_stiff = immersed_mesh.waterplane_integral(immersed_mesh.faces_centers[:,0] - xc)
            elif dof_pair == ("Roll", "Roll"):
                norm_hs_stiff = (
                        -immersed_mesh.waterplane_integral((immersed_mesh.faces_centers[:,1] - yc)**2)
                        + immersed_mesh.volume*(immersed_mesh.center_of_buoyancy[2] - zc) - mass/rho*(self.center_of_mass[2] - zc)
                )
            elif dof_pair in [("Roll", "Pitch"), ("Pitch", "Roll")]:
                norm_hs_stiff = immersed_mesh.waterplane_integral((immersed_mesh.faces_centers[:,0] - xc)
                                                          * (immersed_mesh.faces_centers[:,1] - yc))
            elif dof_pair == ("Roll", "Yaw"):
                norm_hs_stiff = - immersed_mesh.volume*(immersed_mesh.center_of_buoyancy[0] - xc) + mass/rho*(self.center_of_mass[0] - xc)
            elif dof_pair == ("Pitch", "Pitch"):
                norm_hs_stiff = (
                        -immersed_mesh.waterplane_integral((immersed_mesh.faces_centers[:,0] - xc)**2)
                        + immersed_mesh.volume*(immersed_mesh.center_of_buoyancy[2] - zc) - mass/rho*(self.center_of_mass[2] - zc)
                        )
            elif dof_pair == ("Pitch", "Yaw"):
                norm_hs_stiff = - immersed_mesh.volume*(immersed_mesh.center_of_buoyancy[1] - yc) + mass/rho*(self.center_of_mass[1] - yc)
            else:
                norm_hs_stiff = 0.0
        else:
            if self.mass is not None and not np.isclose(self.mass, self.disp_mass(rho=rho), rtol=1e-4):
                raise NotImplementedError(
                        f"Trying to compute the hydrostatic stiffness for dofs {radiating_dof_name} and {influenced_dof_name}"
                        f"of body {self.name}, which is not neutrally buoyant (mass={self.mass}, disp_mass={self.disp_mass(rho=rho)}).\n"
                        f"This case has not been implemented in Capytaine. You need either a single rigid body or a neutrally buoyant body."
                        )

            if np.any(self.mesh.faces_centers[:, 2] > 1e-2) and np.any(influenced_dof_div != 0.0):
                raise NotImplementedError(
                        "When computing hydrostatics of flexible dofs while providing the divergence of the dof, please make sure the mesh is clipped beforehand and provide the divergence only on the immersed faces of the clipped mesh."
                        )

            immersed_self = self.immersed_part()

            # Newman (1994) formula for flexible DOFs
            influenced_dof = np.array(immersed_self.dofs[influenced_dof_name])
            radiating_dof = np.array(immersed_self.dofs[radiating_dof_name])
            influenced_dof_div_array = np.array(influenced_dof_div)

            radiating_dof_normal = immersed_self.dof_normals(radiating_dof)
            z_influenced_dof_div = influenced_dof[:,2] + immersed_self.mesh.faces_centers[:,2] * influenced_dof_div_array
            norm_hs_stiff = immersed_self.mesh.surface_integral( -radiating_dof_normal * z_influenced_dof_div)

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
            if isinstance(influenced_dof, (TranslationDof, RotationDof)):
                return 0.0  # Dummy value that is not actually used afterwards.
            elif divergence is None:
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
            ], compat='no_conflicts', join="outer")

        # Reorder dofs
        K = hs_set.hydrostatic_stiffness.sel(influenced_dof=list(self.dofs.keys()), radiating_dof=list(self.dofs.keys()))
        return K

    def compute_rigid_body_inertia(self, *, rho=1000.0, output_type="body_dofs"):
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


        inertias = 1/self.volume * np.array([
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
        volumic_inertia_matrix = self.mesh.disp_volume * np.array([
            [1.0,     0,       0,       0,            cog[2],       -cog[1]     ],
            [0,       1.0,     0,       -cog[2],      0,            cog[0]      ],
            [0,       0,       1.0,     cog[1],       -cog[0],      0           ],
            [0,       -cog[2], cog[1],  inertias[0],  -inertias[3], -inertias[5]],
            [cog[2],  0,       -cog[0], -inertias[3], inertias[1],  -inertias[4]],
            [-cog[1], cog[0],  0,       -inertias[5], -inertias[4], inertias[2] ],
        ])

        density = rho if self.mass is None else self.mass/self.volume
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

        total_mass_xr = xr.merge([rigid_inertia_matrix_xr, other_dofs_inertia_matrix_xr], compat="override", join="outer").inertia_matrix

        non_rigid_dofs = set(body_dof_names) - set(rigid_dof_names)

        if output_type == "body_dofs":
            if len(non_rigid_dofs) > 0:
                LOG.warning(f"Non-rigid dofs {non_rigid_dofs} detected: their \
inertia coefficients are assigned as NaN.")

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

        full_mesh_vertices = self.mesh.vertices
        coord_max = full_mesh_vertices.max(axis=0)
        coord_min = full_mesh_vertices.min(axis=0)
        full_length, full_breadth, depth = full_mesh_vertices.max(axis=0) - full_mesh_vertices.min(axis=0)

        vertices = self.mesh.immersed_part().vertices
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
        hydrostatics["disp_volume"] = self.disp_volume
        hydrostatics["disp_mass"] = self.disp_mass(rho=rho)
        hydrostatics["center_of_buoyancy"] = self.center_of_buoyancy
        hydrostatics["waterplane_center"] = np.append(self.waterplane_center, 0.0)
        hydrostatics["waterplane_area"] = self.waterplane_area
        hydrostatics["transversal_metacentric_radius"] = self.transversal_metacentric_radius
        hydrostatics["longitudinal_metacentric_radius"] = self.longitudinal_metacentric_radius
        hydrostatics["transversal_metacentric_height"] = self.transversal_metacentric_height
        hydrostatics["longitudinal_metacentric_height"] = self.longitudinal_metacentric_height
        self.hydrostatic_stiffness = hydrostatics["hydrostatic_stiffness"] = self.compute_hydrostatic_stiffness(
            divergence=divergence, rho=rho, g=g)

        hydrostatics["length_overall"] = full_length
        hydrostatics["breadth_overall"] = full_breadth
        hydrostatics["depth"] = depth
        hydrostatics["draught"] = np.abs(coord_min[2])
        hydrostatics["length_at_waterline"] = wl_length
        hydrostatics["breadth_at_waterline"] = wl_breadth
        hydrostatics["length_overall_submerged"] = sub_length
        hydrostatics["breadth_overall_submerged"] = sub_breadth
        if any(dof.lower() in {"surge", "sway", "heave", "roll", "pitch", "yaw"}
               for dof in self.dofs) > 0: # If there is at least one rigid body dof:
            self.inertia_matrix = hydrostatics["inertia_matrix"] = self.compute_rigid_body_inertia(rho=rho)

        return hydrostatics
