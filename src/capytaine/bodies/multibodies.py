from __future__ import annotations

import logging
from itertools import chain, accumulate
from typing import Union, Dict, List, Optional
from functools import cached_property, lru_cache

import numpy as np
import xarray as xr


LOG = logging.getLogger(__name__)


class Multibody:
    def __init__(
        self,
        bodies: List[Union[FloatingBody, Multibody]],
        own_dofs: Optional[Dict[str, np.array]] = None,
        *,
        name: Optional[str] = None
    ):
        self.bodies = bodies

        if len(set(b.name for b in self.bodies)) < len(self.bodies):
            raise ValueError(
                "In Multibody, all component bodies must have a distinct name.\n"
                f"Got: {[b.name for b in self.bodies]}"
                )

        if own_dofs is None:
            self.own_dofs = {}
        else:
            self.own_dofs = own_dofs

        if name is None:
            self.name = '+'.join(b.name for b in self.bodies)
        else:
            self.name = name

        # for matrix_name in ["inertia_matrix", "hydrostatic_stiffness"]:
        #     if all(hasattr(body, matrix_name) for body in bodies):
        #         from scipy.linalg import block_diag
        #         setattr(self, matrix_name, self.add_dofs_labels_to_matrix(
        #                 block_diag(*[getattr(body, matrix_name) for body in bodies])
        #                 ))

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
                mass=self.mass,
                center_of_mass=new_cog,
                name=self.name,
                )

    def __str__(self):
        short_bodies = ', '.join(b.__short_str__() for b in self.bodies)
        short_dofs = '{' + ', '.join('"{}": ...'.format(d) for d in self.own_dofs) + '}'
        return f"Multibody({short_bodies}, own_dofs={short_dofs})"

    def __short_str__(self):
        return str(self)

    def __repr__(self):
        return str(self)

    def _check_dofs_shape_consistency(self):
        # TODO
        ...

    @cached_property
    def minimal_computable_wavelength(self):
        return min(b.minimal_computable_wavelength for b in self.bodies)

    def first_irregular_frequency_estimate(self, *args, **kwargs):
        return min(b.first_irregular_frequency_estimate(*args, **kwargs) for b in self.bodies)

    @cached_property
    def mesh(self):
        return self.bodies[0].mesh.join_meshes(*[b.mesh for b in self.bodies[1:]])

    @cached_property
    def lid_mesh(self):
        if all(body.lid_mesh is None for body in self.bodies):
            return None
        else:
            lid_meshes = [body.lid_mesh.copy() for body in self.bodies if body.lid_mesh is not None]
            joined_lid = lid_meshes[0].join_meshes(*lid_meshes[1:], name=f"{self.name}_lid_mesh")
            return joined_lid

    @cached_property
    def mesh_including_lid(self):
        return self.bodies[0].mesh_including_lid.join_meshes(*[b.mesh_including_lid for b in self.bodies[1:]])

    @cached_property
    def hull_mask(self):
        return np.concatenate([b.hull_mask for b in self.bodies])

    @property
    def nb_dofs(self):
        return sum(b.nb_dofs for b in self.bodies) + len(self.own_dofs)

    @property
    def dofs(self):
        for body in self.bodies:
            body._check_dofs_shape_consistency()

        componenents_dofs = {}
        cum_nb_faces = accumulate(chain([0], (body.mesh.nb_faces for body in self.bodies)))
        total_nb_faces = sum(body.mesh.nb_faces for body in self.bodies)
        for body, nbf in zip(self.bodies, cum_nb_faces):
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
                componenents_dofs[new_dof_name] = new_dof

        return {**componenents_dofs, **self.own_dofs}

    def immersed_part(self, *args, **kwargs):
        return Multibody(
                [b.immersed_part() for b in self.bodies],
                own_dofs=None,  # TODO
                )

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
            fill_value=0.0
        ).sel(
            radiating_dof=list(self.dofs.keys()),
            influenced_dof=list(self.dofs.keys())
        )

    def compute_hydrostatic_stiffness(self, *, rho=1000.0, g=9.81):
        return self._combine_component_matrices([b.compute_hydrostatic_stiffness(rho=rho, g=g) for b in self.bodies])

    def compute_rigid_body_inertia(self, rho=1000.0):
        return self._combine_component_matrices([b.compute_rigid_body_inertia(rho=rho) for b in self.bodies])
