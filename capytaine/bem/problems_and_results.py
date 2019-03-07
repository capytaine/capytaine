#!/usr/bin/env python
# coding: utf-8
"""Definition of the problems to solve with the BEM solver."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

from attr import attrs, attrib, astuple, Factory, asdict

import numpy as np
from scipy.optimize import newton


LOG = logging.getLogger(__name__)


@attrs(cmp=False)
class LinearPotentialFlowProblem:
    """General class of a potential flow problem.

    Stores:

    * the environmental variables (gravity and fluid density),
    * the shape of the domain (position of the free surface and of the sea bottom),
    * the frequency of interest,
    * the meshed floating body,
    * the Neumann boundary conditions on the body.
    """
    default_parameters = {'rho': 1000.0, 'g': 9.81, 'free_surface': 0.0, 'sea_bottom': -np.infty, 'omega': 1.0}

    body = attrib(default=None)
    free_surface = attrib(default=default_parameters['free_surface'])
    sea_bottom = attrib(default=default_parameters['sea_bottom'])
    omega = attrib(default=default_parameters['omega'])
    g = attrib(default=default_parameters['g'])
    rho = attrib(default=default_parameters['rho'])

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

    @body.validator
    def _check_body_position(self, _, body):
        if body is not None:
            if (any(body.mesh.faces_centers[:, 2] > self.free_surface)
                    or any(body.mesh.faces_centers[:, 2] < self.sea_bottom)):
                LOG.warning(f"""The mesh of the body {body.name} is not inside the domain.\n
                                Check the values of free_surface and sea_bottom\n
                                or use body.keep_immersed_part() to clip the mesh.""")

    @omega.validator
    def _check_frequency(self, _, omega):
        if omega in {0, np.infty} and self.depth != np.infty:
            raise NotImplementedError(f"omega={omega} is only implemented for infinite depth.")

    def __str__(self):
        """Do not display default values in str(problem)."""
        parameters = [f"body={self.body.name if self.body is not None else 'None'}",
                      f"omega={self.omega:.3f}",
                      f"depth={self.depth}"]
        try:
            parameters.extend(self._str_other_attributes())
        except AttributeError:
            pass

        if not self.free_surface == self.default_parameters['free_surface']:
            parameters.append(f"free_surface={self.free_surface}")
        if not self.g == self.default_parameters['g']:
            parameters.append(f"g={self.g}")
        if not self.rho == self.default_parameters['rho']:
            parameters.append(f"rho={self.rho}")

        return self.__class__.__name__ + "(" + ', '.join(parameters) + ")"

    def __eq__(self, other):
        if isinstance(other, LinearPotentialFlowProblem):
            return astuple(self)[:6] == astuple(other)[:6]
        else:
            return NotImplemented

    def __lt__(self, other):
        # Arbitrary order. Used for ordering of problems: problems with similar body are grouped together.
        if isinstance(other, LinearPotentialFlowProblem):
            return astuple(self)[:6] < astuple(other)[:6]
        else:
            return NotImplemented

    @property
    def depth(self):
        return self.free_surface - self.sea_bottom

    @property
    def wavenumber(self):
        if self.depth == np.infty or self.omega**2*self.depth/self.g > 20:
            return self.omega**2/self.g
        else:
            return newton(lambda x: x*np.tanh(x) - self.omega**2*self.depth/self.g, x0=1.0)/self.depth

    @property
    def wavelength(self):
        if self.wavenumber == 0.0:
            return np.infty
        else:
            return 2*np.pi/self.wavenumber

    @property
    def period(self):
        if self.omega == 0.0:
            return np.infty
        else:
            return 2*np.pi/self.omega

    @property
    def dimensionless_omega(self):
        if self.depth != np.infty:
            return self.omega**2*self.depth/self.g
        else:
            raise AttributeError("Dimensionless omega is defined only for finite depth problems.")

    @property
    def dimensionless_wavenumber(self):
        if self.depth != np.infty:
            return self.wavenumber*self.depth
        else:
            raise AttributeError("Dimensionless wavenumber is defined only for finite depth problems.")

    @property
    def influenced_dofs(self):
        # TODO: let the user choose the influenced dofs
        return self.body.dofs

    def make_results_container(self):
        return LinearPotentialFlowResult(self)


@attrs(cmp=False)
class DiffractionProblem(LinearPotentialFlowProblem):
    """Particular LinearPotentialFlowProblem whose boundary conditions have
    been computed from an incoming Airy wave."""

    wave_direction = attrib(default=0.0)  # Angle of the incoming wave.
    convention = attrib(default="Nemoh", repr=False)

    def __attrs_post_init__(self):
        from capytaine.bem.airy_waves import airy_waves_velocity
        if self.body is not None:
            self.boundary_condition = -(
                    airy_waves_velocity(self.body.mesh.faces_centers, self, convention=self.convention)
                    * self.body.mesh.faces_normals
            ).sum(axis=1)

            if len(self.body.dofs) == 0:
                LOG.warning(f"The body {self.body.name} used in diffraction problem has no dofs!")

    def _str_other_attributes(self):
        return [f"wave_direction={self.wave_direction:.3f}"]

    def make_results_container(self):
        return DiffractionResult(self)


@attrs(cmp=False)
class RadiationProblem(LinearPotentialFlowProblem):
    """Particular LinearPotentialFlowProblem whose boundary conditions have
    been computed from the degree of freedom of the body."""

    radiating_dof = attrib(default=None)

    def __str__(self):
        """Do not display default values in str(problem)."""
        parameters = [f"body={self.body.name if self.body is not None else 'None'}, "
                      f"omega={self.omega:.3f}, depth={self.depth}, radiating_dof={self.radiating_dof}, "]
        if not self.free_surface == 0.0:
            parameters.append(f"free_surface={self.free_surface}, ")
        if not self.g == 9.81:
            parameters.append(f"g={self.g}, ")
        if not self.rho == 1000:
            parameters.append(f"rho={self.rho}, ")
        return "RadiationProblem(" + ''.join(parameters)[:-2] + ")"

    def __attrs_post_init__(self):
        """Set the boundary condition"""
        if self.body is None:
            self.boundary_condition = None
            return

        if len(self.body.dofs) == 0:
            raise ValueError(f"Body {self.body.name} does not have any degrees of freedom.")

        if self.radiating_dof is None:
            self.radiating_dof = next(iter(self.body.dofs))

        if self.radiating_dof not in self.body.dofs:
            LOG.error(f"In {self}: the radiating degree of freedom {self.radiating_dof} is not one of"
                      f"the degrees of freedom of the body.\n"
                      f"The dofs of the body are {list(self.body.dofs.keys())}")
            raise ValueError("Unrecognized degree of freedom name.")

        dof = self.body.dofs[self.radiating_dof]
        self.boundary_condition = np.sum(dof * self.body.mesh.faces_normals, axis=1)

    def _str_other_attributes(self):
        return [f"radiating_dof={self.radiating_dof}"]

    def make_results_container(self):
        return RadiationResult(self)


@attrs
class LinearPotentialFlowResult:
    problem = attrib()

    sources = attrib(default=None, init=False, repr=False)
    potential = attrib(default=None, init=False, repr=False)

    fs_elevation = attrib(default=Factory(dict), init=False, repr=False)

    __str__ = LinearPotentialFlowProblem.__str__

    def __getattr__(self, name):
        """Direct access to the attributes of the included problem."""
        try:
            return getattr(self.problem, name)
        except AttributeError:
            raise AttributeError(f"{self.__class__} does not have a attribute named {name}.")

    @property
    def settings_dict(self):
        settings = asdict(self.problem)
        # Keep only the name of the body, not the full object.
        settings['body_name'] = self.body.name
        del settings['body']
        # Keep only water_depth  # TODO: Remove.
        settings['water_depth'] = self.free_surface - self.sea_bottom
        del settings['free_surface']
        del settings['sea_bottom']
        return settings


@attrs
class DiffractionResult(LinearPotentialFlowResult):
    forces = attrib(default=Factory(dict), init=False, repr=False)

    def store_force(self, dof, force):
        self.forces[dof] = 1j*self.omega*force

    @property
    def records(self):
        from capytaine.bem.airy_waves import froude_krylov_force
        FK = froude_krylov_force(self.problem)
        return [dict(self.settings_dict, influenced_dof=dof,
                     diffraction_force=self.forces[dof], Froude_Krylov_force=FK[dof])
                for dof in self.influenced_dofs]


@attrs
class RadiationResult(LinearPotentialFlowResult):
    added_masses = attrib(default=Factory(dict), init=False, repr=False)
    radiation_dampings = attrib(default=Factory(dict), init=False, repr=False)

    def store_force(self, dof, force):
        self.added_masses[dof] = force.real
        if self.problem.omega == np.infty:
            self.radiation_dampings[dof] = 0
        else:
            self.radiation_dampings[dof] = self.problem.omega * force.imag

    @property
    def records(self):
        return [dict(self.settings_dict, influenced_dof=dof,
                     added_mass=self.added_masses[dof], radiation_damping=self.radiation_dampings[dof])
                for dof in self.influenced_dofs]
