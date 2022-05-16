#!/usr/bin/env python
# coding: utf-8
"""Definition of the problems to solve with the BEM solver, and the results of this resolution."""
# Copyright (C) 2017-2020 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np
import pandas as pd
from scipy.optimize import newton

from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.bem.airy_waves import airy_waves_velocity, froude_krylov_force

LOG = logging.getLogger(__name__)

_default_parameters = {'rho': 1000.0, 'g': 9.81, 'omega': 1.0,
                      'free_surface': 0.0, 'water_depth': np.infty,
                      'convention': "nemoh"}


class LinearPotentialFlowProblem:
    """General class of a potential flow problem.

    Parameters
    ----------
    body: FloatingBody, optional
        The body interacting with the waves
    free_surface: float, optional
        The position of the free surface (accepted values: 0 and np.infty)
    sea_bottom: float, optional
        The position of the sea bottom
    omega: float, optional
        The frequency of the waves in rad/s
    rho: float, optional
        The density of water in kg/m3 (default: 1000.0)
    g: float, optional
        The acceleration of gravity in m/s2 (default: 9.81)
    boundary_condition: np.ndarray of shape (body.mesh.nb_faces,)
        The Neumann boundary condition on the floating body

    TODO: more consistent use of free_surface and sea_bottom vs. water_depth
    """

    def __init__(self, *,
                 body=None,
                 free_surface=_default_parameters['free_surface'],
                 sea_bottom=-_default_parameters['water_depth'],
                 omega=_default_parameters['omega'],
                 rho=_default_parameters['rho'],
                 g=_default_parameters['g'],
                 boundary_condition=None):

        self.body = body
        self.free_surface = free_surface
        self.sea_bottom = sea_bottom
        self.omega = omega
        self.rho = rho
        self.g = g
        self.boundary_condition = boundary_condition

        self._check_data()

    def _check_data(self):
        """Sanity checks on the data."""

        if self.free_surface not in [0.0, np.infty]:
            raise NotImplementedError(
                f"Free surface is {self.free_surface}. "
                "Only z=0 and z=∞ are accepted values for the free surface position."
            )

        if self.free_surface == np.infty and self.sea_bottom != -np.infty:
            raise NotImplementedError(
                "Problems with a sea bottom but no free surface have not been implemented."
            )

        if self.free_surface < self.sea_bottom:
            raise ValueError("Sea bottom is above the free surface.")

        if self.omega in {0, np.infty} and self.depth != np.infty:
            raise NotImplementedError(
                f"omega={self.omega} is only implemented for infinite depth."
            )

        if self.body is not None:
            if ((isinstance(self.body.mesh, CollectionOfMeshes) and len(self.body.mesh) == 0)
                    or len(self.body.mesh.faces) == 0):
                raise ValueError(f"The mesh of the body {self.body.name} is empty.")

            if (any(self.body.mesh.faces_centers[:, 2] >= self.free_surface)
                    or any(self.body.mesh.faces_centers[:, 2] <= self.sea_bottom)):

                LOG.warning(
                    f"The mesh of the body {self.body.name} is not inside the domain.\n"
                    "Check the position of the free_surface and the sea_bottom\n"
                    "or use body.keep_immersed_part() to clip the mesh."
                )

        if self.boundary_condition is not None:
            if len(self.boundary_condition.shape) != 1:
                raise ValueError("Expected a 1-dimensional array as boundary_condition")

            if self.boundary_condition.shape[0] != self.body.mesh.nb_faces:
                raise ValueError(
                    f"The shape of the boundary condition ({self.boundary_condition.shape})"
                    f"does not match the number of faces of the mesh ({self.body.mesh.nb_faces})."
                )

    @property
    def body_name(self):
        return self.body.name if self.body is not None else 'None'

    def _asdict(self):
        return {"body_name": self.body_name,
                "water_depth": self.water_depth,
                "omega": self.omega,
                "rho": self.rho,
                "g": self.g}

    @staticmethod
    def _group_for_parallel_resolution(problems):
        """Given a list of problems, returns a list of groups of problems, such
        that each group should be executed in the same process to benefit from
        caching.
        """
        problems_params = pd.DataFrame([pb._asdict() for pb in problems])
        groups_of_indices = problems_params.groupby(["body_name", "water_depth", "omega", "rho", "g"]).groups.values()
        groups_of_problems = [[problems[i] for i in grp] for grp in groups_of_indices]
        return groups_of_problems

    def __str__(self):
        """Do not display default values in str(problem)."""
        parameters = [f"body={self.body_name}",
                      f"omega={self.omega:.3f}",
                      f"depth={self.depth}"]
        try:
            parameters.extend(self._str_other_attributes())
        except AttributeError:
            pass

        if not self.free_surface == _default_parameters['free_surface']:
            parameters.append(f"free_surface={self.free_surface}")
        if not self.g == _default_parameters['g']:
            parameters.append(f"g={self.g}")
        if not self.rho == _default_parameters['rho']:
            parameters.append(f"rho={self.rho}")

        return self.__class__.__name__ + "(" + ', '.join(parameters) + ")"

    def _astuple(self):
        return (self.body, self.free_surface, self.sea_bottom, self.omega, self.rho, self.g)

    def __eq__(self, other):
        if isinstance(other, LinearPotentialFlowProblem):
            return self._astuple() == other._astuple()
        else:
            return NotImplemented

    def __lt__(self, other):
        # Arbitrary order. Used for ordering of problems: problems with same body are grouped together.
        if isinstance(other, LinearPotentialFlowProblem):
            return self._astuple()[:6] < other._astuple()[:6]
            # Not the whole tuple, because when using inheriting classes,
            # "radiating_dof" cannot be compared with "wave_direction"
        else:
            return NotImplemented

    @property
    def water_depth(self):
        return self.free_surface - self.sea_bottom

    @property
    def depth(self):
        return self.water_depth

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
        return self.body.dofs if self.body is not None else set()

    def make_results_container(self):
        return LinearPotentialFlowResult(self)


class DiffractionProblem(LinearPotentialFlowProblem):
    """Particular LinearPotentialFlowProblem with boundary conditions
    computed from an incoming Airy wave."""

    def __init__(self, *,
                 body=None,
                 free_surface=_default_parameters['free_surface'],
                 sea_bottom=-_default_parameters['water_depth'],
                 omega=_default_parameters['omega'],
                 rho=_default_parameters['rho'],
                 g=_default_parameters['g'],
                 wave_direction=0.0,
                 convention=_default_parameters['convention']):

        self.wave_direction = wave_direction
        self.convention = convention

        super().__init__(body=body, free_surface=free_surface, sea_bottom=sea_bottom,
                         omega=omega, rho=rho, g=g)

        if not (-2*np.pi-1e-3 <= self.wave_direction <= 2*np.pi+1e-3):
            LOG.warning(f"The value {self.wave_direction} has been provided for the wave direction, and it does not look like an angle in radians. "
                         "The wave direction in Capytaine is defined in radians and not in degrees, so the result might not be what you expect. "
                         "If you were actually giving an angle in radians, use the modulo operator to give a value between -2π and 2π to disable this warning.")

        if self.body is not None:

            self.boundary_condition = -(
                    airy_waves_velocity(self.body.mesh.faces_centers, self, convention=self.convention)
                    * self.body.mesh.faces_normals
            ).sum(axis=1)

            if len(self.body.dofs) == 0:
                LOG.warning(f"The body {self.body.name} used in diffraction problem has no dofs!")

    def _astuple(self):
        return super()._astuple() + (self.wave_direction,)

    def _asdict(self):
        d = super()._asdict()
        d["wave_direction"] = self.wave_direction
        d["convention"] = self.convention
        return d

    def _str_other_attributes(self):
        return [f"wave_direction={self.wave_direction:.3f}"]

    def make_results_container(self):
        return DiffractionResult(self)


class RadiationProblem(LinearPotentialFlowProblem):
    """Particular LinearPotentialFlowProblem whose boundary conditions have
    been computed from the degree of freedom of the body."""

    def __init__(self, *, body=None,
                 free_surface=_default_parameters['free_surface'],
                 sea_bottom=-_default_parameters['water_depth'],
                 omega=_default_parameters['omega'],
                 rho=_default_parameters['rho'],
                 g=_default_parameters['g'],
                 radiating_dof=None):

        self.radiating_dof = radiating_dof

        super().__init__(body=body, free_surface=free_surface, sea_bottom=sea_bottom,
                         omega=omega, rho=rho, g=g)

        if self.body is not None:

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

    def _astuple(self):
        return super()._astuple() + (self.radiating_dof,)

    def _asdict(self):
        d = super()._asdict()
        d["radiating_dof"] = self.radiating_dof
        return d

    def _str_other_attributes(self):
        return [f"radiating_dof={self.radiating_dof}"]

    def make_results_container(self):
        return RadiationResult(self)


class LinearPotentialFlowResult:

    def __init__(self, problem):
        self.problem = problem

        self.sources = None
        self.potential = None
        self.fs_elevation = {}

        # Copy data from problem
        self.body               = self.problem.body
        self.free_surface       = self.problem.free_surface
        self.sea_bottom         = self.problem.sea_bottom
        self.omega              = self.problem.omega
        self.rho                = self.problem.rho
        self.g                  = self.problem.g
        self.boundary_condition = self.problem.boundary_condition
        self.water_depth        = self.problem.water_depth
        self.depth              = self.problem.depth
        self.wavenumber         = self.problem.wavenumber
        self.wavelength         = self.problem.wavelength
        self.period             = self.problem.period
        self.body_name          = self.problem.body_name
        self.influenced_dofs    = self.problem.influenced_dofs

    __str__ = LinearPotentialFlowProblem.__str__


class DiffractionResult(LinearPotentialFlowResult):

    def __init__(self, problem):
        super().__init__(problem)
        self.wave_direction = self.problem.wave_direction
        self.forces = {}

    def store_force(self, dof, force):
        self.forces[dof] = 1j*self.omega*force

    @property
    def records(self):
        params = self.problem._asdict()
        FK = froude_krylov_force(self.problem)
        return [dict(**params,
                     influenced_dof=dof,
                     diffraction_force=self.forces[dof],
                     Froude_Krylov_force=FK[dof])
                for dof in self.influenced_dofs]


class RadiationResult(LinearPotentialFlowResult):

    def __init__(self, problem):
        super().__init__(problem)
        self.radiating_dof = self.problem.radiating_dof
        self.added_masses = {}
        self.radiation_dampings = {}

    def store_force(self, dof, force):
        self.added_masses[dof] = force.real
        if self.problem.omega == np.infty:
            self.radiation_dampings[dof] = 0
        else:
            self.radiation_dampings[dof] = self.problem.omega * force.imag

    @property
    def records(self):
        params = self.problem._asdict()
        return [dict(params,
                     influenced_dof=dof,
                     added_mass=self.added_masses[dof],
                     radiation_damping=self.radiation_dampings[dof])
                for dof in self.influenced_dofs]
