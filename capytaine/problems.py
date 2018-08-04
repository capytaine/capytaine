#!/usr/bin/env python
# coding: utf-8
"""Definition of the problems to solve with the BEM solver."""
# This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import logging

from attr import attrs, attrib, astuple
import numpy as np

from capytaine._Wavenumber import invert_xtanhx
from capytaine.results import LinearPotentialFlowResult, DiffractionResult, RadiationResult
from capytaine.tools.Airy_wave import Airy_wave_velocity


LOG = logging.getLogger(__name__)


@attrs(slots=True, cmp=False)
class LinearPotentialFlowProblem:
    """General class of a potential flow problem.

    Stores:
    * the environmental variables (gravity and fluid density),
    * the shape of the domain (position of the free surface and of the sea bottom),
    * the frequency of interest,
    * the meshed floating body,
    * the Neumann boundary conditions on the body.
    """

    body = attrib(default=None)
    free_surface = attrib(default=0.0)
    sea_bottom = attrib(default=-np.infty)
    omega = attrib(default=1.0)
    g = attrib(default=9.81)
    rho = attrib(default=1000.0)

    boundary_condition = attrib(default=None, repr=False)

    @free_surface.validator
    def _check_free_surface(self, _, free_surface):
        if free_surface not in [0, np.infty]:
            raise NotImplementedError(
                "Only z=0 and z=∞ are accepted values for the free surface position at the moment.")
        elif free_surface == np.infty and self.sea_bottom != -np.infty:
            raise NotImplementedError(
                "The case without free surface but with a sea bottom has not been implemented yet.")

    @sea_bottom.validator
    def _check_depth(self, _, sea_bottom):
        if self.free_surface < sea_bottom:
            raise ValueError("Sea bottom is above the free surface.")

    # @body.validator
    # def _check_body_position(self, _, body):
    #     if body is not None:
    #         if (any(body.mesh.vertices[:, 2] > self.free_surface + 1e-3)
    #                 or any(body.mesh.vertices[:, 2] < self.sea_bottom - 1e-3)):
    #             LOG.warning(f"""The mesh of the body {body.name} is not inside the domain.\n
    #                             Check the values of free_surface and sea_bottom\n
    #                             or use body.get_immersed_part() to clip the mesh.""")

    @boundary_condition.validator
    def _check_size_of_boundary_condition(self, _, bc):
        if self.body is None:
            if bc is not None:
                LOG.warning(f"""The problem {self} has no body but has a boundary condition.""")
        else:
            if bc is not None and len(bc) != self.body.mesh.nb_faces:
                LOG.warning(f"""The size of the boundary condition in {self} does not match the
                            number of faces in the body.""")

    def __lt__(self, other):
        if isinstance(other, LinearPotentialFlowProblem):
            return astuple(self)[:6] < astuple(other)[:6]
        else:
            return NotImplemented

    @property
    def depth(self):
        return self.free_surface - self.sea_bottom

    @property
    def wavenumber(self):
        # TODO: Store the value?
        if self.depth == np.infty or self.omega**2*self.depth/self.g > 20:
            return self.omega**2/self.g
        else:
            return invert_xtanhx(self.omega**2*self.depth/self.g)/self.depth

    @property
    def wavelength(self):
        return 2*np.pi/self.wavenumber

    @property
    def period(self):
        return 2*np.pi/self.omega

    @property
    def dimensionless_omega(self):
        return self.omega**2*self.depth/self.g

    @property
    def dimensionless_wavenumber(self):
        return self.wavenumber*self.depth

    @property
    def influenced_dofs(self):
        # TODO: let the user choose the influenced dofs
        return self.body.dofs

    def make_results_container(self):
        return LinearPotentialFlowResult(self)


@attrs(slots=True, cmp=False)
class DiffractionProblem(LinearPotentialFlowProblem):
    """Particular LinearPotentialFlowProblem whose boundary conditions have
    been computed from an incoming Airy wave."""
    # body = attrib(default=None)
    angle = attrib(default=0.0)  # Angle of the incoming wave.
    boundary_condition = attrib(default=None, init=False, repr=False)
    convention = attrib(default="Nemoh", repr=False)

    # @body.validator
    # def _check_dofs(self, attribute, body):
    #     if body is not None and len(self.body.dofs) == 0:
    #         LOG.warning(f"In {self}: the body has no degrees of freedom defined.\n"
    #                     f"The problem will be solved but the Froude-Krylov forces won't be computed.")

    def __str__(self):
        parameters = [f"body={self.body.name}, omega={self.omega:.3f}, depth={self.depth}, angle={self.angle:.3f}, "]
        if not self.free_surface == 0.0:
            parameters.append(f"free_surface={self.free_surface}, ")
        if not self.g == 9.81:
            parameters.append(f"g={self.g}, ")
        if not self.rho == 1000:
            parameters.append(f"rho={self.rho}, ")
        return "DiffractionProblem(" + ''.join(parameters)[:-2] + ")"

    def __attrs_post_init__(self):
        if self.body is not None:
            self.boundary_condition = -(
                Airy_wave_velocity(self.body.mesh.faces_centers, self, convention=self.convention)
                * self.body.mesh.faces_normals
            ).sum(axis=1)

    def make_results_container(self):
        return DiffractionResult(self)


@attrs(slots=True, cmp=False)
class RadiationProblem(LinearPotentialFlowProblem):
    """Particular LinearPotentialFlowProblem whose boundary conditions have
    been computed from the degree of freedom of the body."""
    # body = attrib(default=None)
    radiating_dof = attrib(default=None)
    boundary_condition = attrib(default=None, init=False, repr=False)

    # @body.validator
    # def _check_dofs(self, attribute, body):
    #     if len(body.dofs) == 0:
    #         LOG.error(f"In {self}: the body has no degrees of freedom defined.")
    #         raise ValueError("The body in a radiation problem needs to have degrees of freedom")

    def __str__(self):
        parameters = [f"body={self.body.name}, omega={self.omega:.3f}, depth={self.depth}, radiating_dof={self.radiating_dof}, "]
        if not self.free_surface == 0.0:
            parameters.append(f"free_surface={self.free_surface}, ")
        if not self.g == 9.81:
            parameters.append(f"g={self.g}, ")
        if not self.rho == 1000:
            parameters.append(f"rho={self.rho}, ")
        return "RadiationProblem(" + ''.join(parameters)[:-2] + ")"

    def __attrs_post_init__(self):
        """Set the boundary condition"""
        if self.radiating_dof is None:
            self.radiating_dof = next(iter(self.body.dofs))
            dof = self.body.dofs[self.radiating_dof]
            self.boundary_condition = np.sum(dof * self.body.mesh.faces_normals, axis=1)
        elif self.radiating_dof in self.body.dofs:
            dof = self.body.dofs[self.radiating_dof]
            self.boundary_condition = np.sum(dof * self.body.mesh.faces_normals, axis=1)
        else:
            LOG.error(f"In {self}: the radiating degree of freedom {self.radiating_dof} is not one of"
                      f"the degrees of freedom of the body.\n"
                      f"The dofs of the body are {list(self.body.dofs.keys())}")
            raise ValueError("Unrecognized degree of freedom name.")

    def make_results_container(self):
        return RadiationResult(self)

