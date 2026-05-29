from __future__ import annotations

import logging
from typing import Union, List, Optional, Literal
from functools import cached_property, lru_cache

import numpy as np
import xarray as xr

from capytaine.bodies.dofs import (
    AbstractDof,
    DofOnSubmesh,
)
from capytaine.bodies.abstract_bodies import AbstractBody
from capytaine.bodies.bodies import FloatingBody

LOG = logging.getLogger(__name__)


class Multibody(AbstractBody):
    def __init__(
        self,
        bodies: List[Union[FloatingBody, Multibody]],
        # own_dofs: Optional[Dict[str, np.array]] = None,
        *,
        name: Optional[str] = None
    ):

        # `bodies` might be FloatingBody or Multibody, we extract a list of FloatingBody.
        self.bodies: List[FloatingBody] = sum((b.bodies for b in bodies), [])

        if len(set(b.name for b in self.bodies)) < len(self.bodies):
            raise ValueError(
                "In Multibody, all component bodies must have a distinct name.\n"
                f"Got: {[b.name for b in self.bodies]}"
                )

        # if own_dofs is None:
        #     self.own_dofs = {}
        # else:
        #     self.own_dofs = own_dofs

        self._name = name

        # MESH MERGING
        # Need to merge the meshes of the component bodies while keeping track of the origin of each panel.
        self.mesh, self.body_masks = self.bodies[0].mesh.join_meshes(
            *[b.mesh for b in self.bodies[1:]],
            return_masks=True
        )

        # DOF MERGING
        for body in self.bodies:
            body._check_dofs_shape_consistency()
        componenents_dofs = {}
        total_nb_faces = self.mesh.nb_faces
        for body, body_mask in zip(self.bodies, self.body_masks):
            for name, dof in body.dofs.items():
                if isinstance(dof, AbstractDof):
                    new_dof = DofOnSubmesh(dof, body_mask)
                else:
                    new_dof = np.zeros((total_nb_faces, 3))
                    new_dof[body_mask, :] = dof

                if '__' not in name:
                    new_dof_name = '__'.join([body.name, name])
                else:
                    # The body is probably a combination of bodies already.
                    # So for the associativity of the + operation,
                    # it is better to keep the same name.
                    new_dof_name = name
                componenents_dofs[new_dof_name] = new_dof

        self.dofs = componenents_dofs

        # LIDS MERGING
        if all(body.lid_mesh is None for body in self.bodies):
            self.lid_mesh = None
            self.mesh_including_lid = self.mesh
            self.hull_mask = np.full((self.mesh.nb_faces,), True)
        else:
            lid_meshes = [body.lid_mesh.copy() for body in self.bodies if body.lid_mesh is not None]
            joined_lid = lid_meshes[0].join_meshes(*lid_meshes[1:], name=f"{self.name}_lid_mesh")
            self.lid_mesh = joined_lid
            self.mesh_including_lid, masks = self.mesh.join_meshes(self.lid_mesh, return_masks=True)
            self.hull_mask = masks[0]

        # HYDROSTATICS MERGING
        # Keep legacy behavior of former mesh joining
        for matrix_name in ["inertia_matrix", "hydrostatic_stiffness"]:
            if all(hasattr(body, matrix_name) for body in bodies):
                from scipy.linalg import block_diag
                setattr(self, matrix_name, self.add_dofs_labels_to_matrix(
                        block_diag(*[getattr(body, matrix_name) for body in bodies])
                        ))

        LOG.debug(f"New multibody: {self.__str__()}.")

    @lru_cache
    def as_FloatingBody(self):
        from capytaine.bodies.bodies import FloatingBody
        if all(body.mass is not None for body in self.bodies):
            total_mass = sum(body.mass for body in self.bodies)
        else:
            total_mass = None

        if (all(body.mass is not None for body in self.bodies)
                and all(body.center_of_mass is not None for body in self.bodies)):
            new_cog = sum(body.mass*np.asarray(body.center_of_mass) for body in self.bodies)/total_mass
        else:
            new_cog = None

        return FloatingBody(
                mesh=self.mesh,
                dofs=self.dofs,
                lid_mesh=self.lid_mesh,
                mass=total_mass,
                center_of_mass=new_cog,
                name=self.name,
                )

    @property
    def name(self):
        if self._name is None:
            return '+'.join(b.name for b in self.bodies)
        else:
            return self._name

    def __str__(self):
        short_bodies = ', '.join(b.__short_str__() for b in self.bodies)
        if self._name is not None:
            name_str = f", name={self._name}"
        else:
            name_str = ""
        return f"Multibody([{short_bodies}]{name_str})"

    def __short_str__(self):
        if self._name is not None:
            return f"Multibody(..., name={self._name})"
        else:
            return str(self)

    def _check_dofs_shape_consistency(self):
        # TODO
        ...

    @cached_property
    def minimal_computable_wavelength(self):
        return min(b.minimal_computable_wavelength for b in self.bodies)

    def first_irregular_frequency_estimate(self, *args, **kwargs):
        return min(b.first_irregular_frequency_estimate(*args, **kwargs) for b in self.bodies)

    @property
    def nb_dofs(self):
        return len(self.dofs)

    def immersed_part(self, *args, **kwargs):
        new_multibody = Multibody(
                [b.immersed_part() for b in self.bodies],
                # own_dofs=None,  # TODO
                )
        if hasattr(self, 'inertia_matrix'):
            new_multibody.inertia_matrix = self.inertia_matrix
        if hasattr(self, 'hydrostatic_stiffness'):
            new_multibody.hydrostatic_stiffness = self.hydrostatic_stiffness
        return new_multibody

    def integrate_pressure(self, pressure):
        return self.as_FloatingBody().integrate_pressure(pressure)

    @cached_property
    def center_of_buoyancy(self):
        return {b.name: b.center_of_buoyancy for b in self.bodies}

    @cached_property
    def center_of_mass(self):
        return {b.name: b.center_of_mass for b in self.bodies}

    @cached_property
    def volume(self):
        return {b.name: b.volume for b in self.bodies}

    @cached_property
    def mass(self):
        return {b.name: b.mass for b in self.bodies}

    def _combine_component_matrices(self, matrices):
        for m, b in zip(matrices, self.bodies):
            m.coords['radiating_dof'] = np.array([b.name + '__' + k for k in m.coords['radiating_dof'].values])
            m.coords['influenced_dof'] = np.array([b.name + '__' + k for k in m.coords['influenced_dof'].values])

        return xr.concat(
            matrices,
            dim="radiating_dof",
            fill_value=0.0,
            join="outer"
        ).sel(
            radiating_dof=list(self.dofs.keys()),
            influenced_dof=list(self.dofs.keys())
        )

    def compute_hydrostatic_stiffness(self, *, rho=1000.0, g=9.81):
        return self._combine_component_matrices([b.compute_hydrostatic_stiffness(rho=rho, g=g) for b in self.bodies])

    def compute_rigid_body_inertia(self, rho=1000.0):
        return self._combine_component_matrices([b.compute_rigid_body_inertia(rho=rho) for b in self.bodies])

    # --- Geometric transforms ---

    def copy(self, name=None) -> Multibody:
        return Multibody(
            [b.copy(name=b.name) for b in self.bodies],
            name=name,
        )

    def rename(self, name: 'str') -> Multibody:
        return self.copy(name=name)

    def translated(self, shift, *, name=None) -> Multibody:
        return Multibody(
            [b.translated(shift, name=b.name) for b in self.bodies],
            name=name,
        )

    def rotated_with_matrix(self, R, *, name=None) -> Multibody:
        return Multibody(
            [b.rotated_with_matrix(R, name=b.name) for b in self.bodies],
            name=name,
        )

    def mirrored(self, plane: Literal['xOz', 'yOz']) -> Multibody:
        mirrored_bodies = []
        for b in self.bodies:
            mb = b.mirrored(plane)
            mb.name = b.name
            mirrored_bodies.append(mb)
        return Multibody(mirrored_bodies)

    def clipped(self, *, origin, normal, name=None) -> Multibody:
        return Multibody(
            [b.clipped(origin=origin, normal=normal, name=b.name) for b in self.bodies],
            name=name,
        )
