"""Definition of the problems to solve with the BEM solver, and the results of this resolution."""
# Copyright (C) 2017-2023 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging

import numpy as np
import pandas as pd
from scipy.optimize import newton

from capytaine.tools.deprecation_handling import _get_water_depth
from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.bem.airy_waves import airy_waves_velocity, froude_krylov_force
from capytaine.tools.symbolic_multiplication import SymbolicMultiplication

LOG = logging.getLogger(__name__)

_default_parameters = {'rho': 1000.0, 'g': 9.81, 'omega': 1.0,
                      'free_surface': 0.0, 'water_depth': np.inf,
                       'wave_direction': 0.0, 'forward_speed': 0.0}



class LinearPotentialFlowProblem:
    """General class of a potential flow problem.

    At most one of the following parameter must be provided: omega, period, wavenumber or wavelength.
    Internally only omega is stored, hence setting another parameter can lead to small rounding errors.

    Parameters
    ----------
    body: FloatingBody, optional
        The body interacting with the waves
    free_surface: float, optional
        The position of the free surface (accepted values: 0 and np.inf)
    water_depth: float, optional
        The depth of water in m (default: np.inf)
    sea_bottom: float, optional
        The position of the sea bottom (deprecated: please prefer setting water_depth)
    omega: float, optional
        The angular frequency of the waves in rad/s
    period: float, optional
        The period of the waves in s
    wavenumber: float, optional
        The angular wave number of the waves in rad/m
    wavelength: float, optional
        The wave length of the waves in m
    forward_speed: float, optional
        The speed of the body (in m/s, in the x direction, default: 0.0)
    rho: float, optional
        The density of water in kg/m3 (default: 1000.0)
    g: float, optional
        The acceleration of gravity in m/s2 (default: 9.81)
    boundary_condition: np.ndarray of shape (body.mesh.nb_faces,), optional
        The Neumann boundary condition on the floating body
    """

    def __init__(self, *,
                 body=None,
                 free_surface=_default_parameters['free_surface'],
                 water_depth=None, sea_bottom=None,
                 omega=None, period=None, wavenumber=None, wavelength=None,
                 forward_speed=_default_parameters['forward_speed'],
                 rho=_default_parameters['rho'],
                 g=_default_parameters['g'],
                 wave_direction=_default_parameters['wave_direction'],
                 boundary_condition=None):

        self.body = body
        self.free_surface = float(free_surface)
        self.rho = float(rho)
        self.g = float(g)
        self.forward_speed = float(forward_speed)
        self.wave_direction = float(wave_direction)  # Required for (diffraction problem) and (radiation problems with forward speed).

        self.boundary_condition = boundary_condition

        self.water_depth = _get_water_depth(free_surface, water_depth, sea_bottom, default_water_depth=_default_parameters["water_depth"])
        self.omega, self.period, self.wavenumber, self.wavelength, self.provided_freq_type = \
                self._get_frequencies(omega=omega, period=period, wavenumber=wavenumber, wavelength=wavelength)

        self._check_data()

        if forward_speed != 0.0:
            dopplered_omega = self.omega - self.wavenumber*self.forward_speed*np.cos(self.wave_direction)
            self.encounter_omega, self.encounter_period, self.encounter_wavenumber, self.encounter_wavelength, _ = \
                    self._get_frequencies(omega=abs(dopplered_omega))

            if dopplered_omega >= 0.0:
                self.encounter_wave_direction = self.wave_direction
            else:
                self.encounter_wave_direction = self.wave_direction + np.pi
        else:
            self.encounter_omega = self.omega
            self.encounter_period = self.period
            self.encounter_wavenumber = self.wavenumber
            self.encounter_wavelength = self.wavelength
            self.encounter_wave_direction = self.wave_direction


    def _get_frequencies(self, *, omega=None, period=None, wavenumber=None, wavelength=None):
        frequency_data = dict(omega=omega, period=period, wavenumber=wavenumber, wavelength=wavelength)
        nb_provided_frequency_data = 4 - list(frequency_data.values()).count(None)

        if nb_provided_frequency_data > 1:
            raise ValueError("Settings a problem requires at most one of the following: omega (angular frequency) OR period OR wavenumber OR wavelength.\n"
                             "Received {} of them: {}".format(nb_provided_frequency_data, {k: v for k, v in frequency_data.items() if v is not None}))

        if nb_provided_frequency_data == 0:
            provided_freq_type = 'omega'
            frequency_data = {'omega': _default_parameters['omega']}
        else:
            provided_freq_type = [k for (k, v) in frequency_data.items() if v is not None][0]

        if ((float(frequency_data[provided_freq_type]) == 0.0 and provided_freq_type in {'omega', 'wavenumber'})
            or (float(frequency_data[provided_freq_type]) == np.inf and provided_freq_type in {'period', 'wavelength'})):
                omega = SymbolicMultiplication("0")
                wavenumber = SymbolicMultiplication("0")
                period = SymbolicMultiplication("∞")
                wavelength = SymbolicMultiplication("∞")
        elif ((float(frequency_data[provided_freq_type]) == 0.0 and provided_freq_type in {'period', 'wavelength'})
            or (float(frequency_data[provided_freq_type]) == np.inf and provided_freq_type in {'omega', 'wavenumber'})):
                omega = SymbolicMultiplication("∞")
                wavenumber = SymbolicMultiplication("∞")
                period = SymbolicMultiplication("0")
                wavelength = SymbolicMultiplication("0")
        else:

            if provided_freq_type in {'omega', 'period'}:
                if provided_freq_type == 'omega':
                    omega = frequency_data['omega']
                    period = 2*np.pi/omega
                else:  # provided_freq_type is 'period'
                    period = frequency_data['period']
                    omega = 2*np.pi/period

                if self.water_depth == np.inf:
                    wavenumber = omega**2/self.g
                else:
                    wavenumber = newton(lambda k: k*np.tanh(k*self.water_depth) - omega**2/self.g, x0=1.0)
                wavelength = 2*np.pi/wavenumber

            else:  # provided_freq_type is 'wavelength' or 'wavenumber'
                if provided_freq_type == 'wavelength':
                    wavelength = frequency_data['wavelength']
                    wavenumber = 2*np.pi/wavelength
                else:  # provided_freq_type is 'wavenumber'
                    wavenumber = frequency_data['wavenumber']
                    wavelength = 2*np.pi/wavenumber

                omega = np.sqrt(self.g*wavenumber*np.tanh(wavenumber*self.water_depth))
                period = 2*np.pi/omega

        return omega, period, wavenumber, wavelength, provided_freq_type

    def _check_data(self):
        """Sanity checks on the data."""

        if self.free_surface not in {0.0, np.inf}:
            raise NotImplementedError(
                f"Free surface is {self.free_surface}. "
                "Only z=0 and z=∞ are accepted values for the free surface position."
            )

        if not (-2*np.pi-1e-3 <= self.wave_direction <= 2*np.pi+1e-3):
            LOG.warning(f"The value {self.wave_direction} has been provided for the wave direction, and it does not look like an angle in radians. "
                         "The wave direction in Capytaine is defined in radians and not in degrees, so the result might not be what you expect. "
                         "If you were actually giving an angle in radians, use the modulo operator to give a value between -2π and 2π to disable this warning.")

        if self.free_surface == np.inf and self.water_depth != np.inf:

            raise NotImplementedError(
                "Problems with a sea bottom but no free surface have not been implemented."
            )

        if self.water_depth < 0.0:
            raise ValueError("`water_depth` should be strictly positive (provided water depth: {self.water_depth}).")

        if float(self.omega) in {0, np.inf}:
            if self.water_depth != np.inf:
                LOG.warning(
                        f"Default Green function allows for {self.provided_freq_type}={float(self.__getattribute__(self.provided_freq_type))} only for infinite depth (provided water depth: {self.water_depth})."
                        )

            if self.forward_speed != 0.0:
                raise NotImplementedError(
                        f"omega={float(self.omega)} is only implemented without forward speed (provided forward speed: {self.forward_speed})."
                )


        if self.body is not None:
            if ((isinstance(self.body.mesh, CollectionOfMeshes) and len(self.body.mesh) == 0)
                    or len(self.body.mesh.faces) == 0):
                raise ValueError(f"The mesh of the body {self.body.__short_str__()} is empty.")

            panels_above_fs = self.body.mesh.faces_centers[:, 2] >= self.free_surface + 1e-8
            panels_below_sb = self.body.mesh.faces_centers[:, 2] <= -self.water_depth
            if (any(panels_above_fs) or any(panels_below_sb)):

                if not any(panels_below_sb):
                    issue = f"{np.count_nonzero(panels_above_fs)} panels above the free surface"
                elif not any(panels_above_fs):
                    issue = f"{np.count_nonzero(panels_below_sb)} panels below the sea bottom"
                else:
                    issue = (f"{np.count_nonzero(panels_above_fs)} panels above the free surface " +
                             f"and {np.count_nonzero(panels_below_sb)} panels below the sea bottom")

                LOG.warning(
                        f"The mesh of the body {self.body.__short_str__()} has {issue}.\n" +
                        "It has been clipped to fit inside the domain.\n" +
                        "To remove this warning, clip the mesh manually with the `immersed_part()` method."
                )

                self.body = self.body.immersed_part(free_surface=self.free_surface,
                                                    water_depth=self.water_depth)

        if self.boundary_condition is not None:
            if len(self.boundary_condition.shape) != 1:
                raise ValueError(f"Expected a 1-dimensional array as boundary_condition. Provided boundary condition's shape: {self.boundary_condition.shape}.")

            if self.boundary_condition.shape[0] != self.body.mesh_including_lid.nb_faces:
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
                "omega": float(self.omega),
                "encounter_omega": float(self.encounter_omega),
                "period": float(self.period),
                "wavelength": float(self.wavelength),
                "wavenumber": float(self.wavenumber),
                "forward_speed": self.forward_speed,
                "wave_direction": self.wave_direction,
                "encounter_wave_direction": self.encounter_wave_direction,
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
        parameters = [f"body={self.body.__short_str__() if self.body is not None else None}",
                      f"{self.provided_freq_type}={float(self.__getattribute__(self.provided_freq_type)):.3f}",
                      f"water_depth={self.water_depth}"]

        if not self.forward_speed == _default_parameters['forward_speed']:
            parameters.append(f"forward_speed={self.forward_speed:.3f}")

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

    def __repr__(self):
        return self.__str__()

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    def __rich_repr__(self):
        yield "body", self.body, None
        yield self.provided_freq_type, self.__getattribute__(self.provided_freq_type)
        yield "water_depth", self.water_depth, _default_parameters["water_depth"]
        try:
            yield from self._specific_rich_repr()
        except:
            pass
        yield "g", self.g, _default_parameters["g"]
        yield "rho", self.rho, _default_parameters["rho"]

    def _astuple(self):
        return (self.body, self.free_surface, self.water_depth,
                float(self.omega), float(self.period), float(self.wavenumber), float(self.wavelength),
                self.forward_speed, self.rho, self.g)

    def __eq__(self, other):
        if isinstance(other, LinearPotentialFlowProblem):
            return self._astuple() == other._astuple()
        else:
            return NotImplemented

    def __lt__(self, other):
        # Arbitrary order. Used for ordering of problems: problems with same body are grouped together.
        if isinstance(other, LinearPotentialFlowProblem):
            return self._astuple()[:9] < other._astuple()[:9]
            # Not the whole tuple, because when using inheriting classes,
            # "radiating_dof" cannot be compared with "wave_direction"
        else:
            return NotImplemented

    @property
    def depth(self):
        return self.water_depth

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
                 water_depth=None, sea_bottom=None,
                 omega=None, period=None, wavenumber=None, wavelength=None,
                 forward_speed=_default_parameters['forward_speed'],
                 rho=_default_parameters['rho'],
                 g=_default_parameters['g'],
                 wave_direction=_default_parameters['wave_direction']):

        super().__init__(body=body, free_surface=free_surface, water_depth=water_depth, sea_bottom=sea_bottom,
                         omega=omega, period=period, wavenumber=wavenumber, wavelength=wavelength, wave_direction=wave_direction,
                         forward_speed=forward_speed, rho=rho, g=g)

        if float(self.omega) in {0.0, np.inf}:
            raise NotImplementedError("DiffractionProblem does not support zero or infinite frequency.")

        if self.body is not None:

            self.boundary_condition = -(
                    airy_waves_velocity(self.body.mesh.faces_centers, self)
                    * self.body.mesh.faces_normals
            ).sum(axis=1)
            # Note that even with forward speed, this is computed based on the
            # frequency and not the encounter frequency.

            if self.body.lid_mesh is not None:
                self.boundary_condition = np.concatenate([self.boundary_condition, np.zeros(self.body.lid_mesh.nb_faces)])

            if len(self.body.dofs) == 0:
                LOG.warning(f"The body {self.body.name} used in diffraction problem has no dofs!")

    def _str_other_attributes(self):
        return [f"wave_direction={self.wave_direction:.3f}"]

    def _specific_rich_repr(self):
        yield "wave_direction", self.wave_direction, _default_parameters["wave_direction"]

    def make_results_container(self, *args, **kwargs):
        return DiffractionResult(self, *args, **kwargs)


class RadiationProblem(LinearPotentialFlowProblem):
    """Particular LinearPotentialFlowProblem whose boundary conditions have
    been computed from the degree of freedom of the body."""

    def __init__(self, *, body=None,
                 free_surface=_default_parameters['free_surface'],
                 water_depth=None, sea_bottom=None,
                 omega=None, period=None, wavenumber=None, wavelength=None,
                 forward_speed=_default_parameters['forward_speed'],
                 wave_direction=_default_parameters['wave_direction'],
                 rho=_default_parameters['rho'],
                 g=_default_parameters['g'],
                 radiating_dof=None):

        self.radiating_dof = radiating_dof

        super().__init__(body=body, free_surface=free_surface, water_depth=water_depth, sea_bottom=sea_bottom,
                         omega=omega, period=period, wavenumber=wavenumber, wavelength=wavelength,
                         wave_direction=wave_direction, forward_speed=forward_speed, rho=rho, g=g)

        if self.body is not None:

            if len(self.body.dofs) == 0:
                raise ValueError(f"Body {self.body.name} does not have any degrees of freedom.")

            if self.radiating_dof is None:
                self.radiating_dof = next(iter(self.body.dofs))

            if self.radiating_dof not in self.body.dofs:
                raise ValueError(f"In {self}:\n"
                                 f"the radiating dof {repr(self.radiating_dof)} is not one of the degrees of freedom of the body.\n"
                                 f"The dofs of the body are {list(self.body.dofs.keys())}")

            dof = self.body.dofs[self.radiating_dof]

            self.boundary_condition = -1j * self.encounter_omega * np.sum(dof * self.body.mesh.faces_normals, axis=1)

            if self.forward_speed != 0.0:
                if self.radiating_dof.lower() == "pitch":
                    ddofdx_dot_n = np.array([nz for (nx, ny, nz) in self.body.mesh.faces_normals])
                elif self.radiating_dof.lower() == "yaw":
                    ddofdx_dot_n = np.array([-ny for (nx, ny, nz) in self.body.mesh.faces_normals])
                elif self.radiating_dof.lower() in {"surge", "sway", "heave", "roll"}:
                    ddofdx_dot_n = 0.0
                else:
                    raise NotImplementedError(
                            "Radiation problem with forward speed is currently only implemented for a single rigid body.\n"
                            "Only radiating dofs with name in {'Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'} are supported.\n"
                            f"Got instead `radiating_dof={self.radiating_dof}`"
                            )
                self.boundary_condition += self.forward_speed * ddofdx_dot_n

            if self.body.lid_mesh is not None:
                self.boundary_condition = np.concatenate([self.boundary_condition, np.zeros(self.body.lid_mesh.nb_faces)])


    def _astuple(self):
        return super()._astuple() + (self.radiating_dof,)

    def _asdict(self):
        d = super()._asdict()
        d["radiating_dof"] = self.radiating_dof
        return d

    def _str_other_attributes(self):
        if self.forward_speed != 0.0:
            return [f"wave_direction={self.wave_direction:.3f}, radiating_dof=\'{self.radiating_dof}\'"]
        else:
            return [f"radiating_dof=\'{self.radiating_dof}\'"]

    def _specific_rich_repr(self):
        yield "radiating_dof", self.radiating_dof

    def make_results_container(self, *args, **kwargs):
        return RadiationResult(self, *args, **kwargs)


class LinearPotentialFlowResult:

    def __init__(self, problem, forces=None, sources=None, potential=None, pressure=None):
        self.problem = problem

        self.forces = forces if forces is not None else {}
        self.sources = sources
        self.potential = potential
        self.pressure = pressure

        self.fs_elevation = {}  # Only used in legacy `get_free_surface_elevation`. To be removed?

        # Copy data from problem
        self.body               = self.problem.body
        self.free_surface       = self.problem.free_surface
        self.omega              = self.problem.omega
        self.period             = self.problem.period
        self.wavenumber         = self.problem.wavenumber
        self.wavelength         = self.problem.wavelength
        self.forward_speed      = self.problem.forward_speed
        self.wave_direction     = self.problem.wave_direction
        self.encounter_omega    = self.problem.encounter_omega
        self.encounter_period   = self.problem.encounter_period
        self.encounter_wavenumber = self.problem.encounter_wavenumber
        self.encounter_wavelength = self.problem.encounter_wavelength
        self.encounter_wave_direction = self.problem.encounter_wave_direction
        self.rho                = self.problem.rho
        self.g                  = self.problem.g
        self.boundary_condition = self.problem.boundary_condition
        self.water_depth        = self.problem.water_depth
        self.depth              = self.problem.water_depth
        self.provided_freq_type = self.problem.provided_freq_type
        self.body_name          = self.problem.body_name
        self.influenced_dofs    = self.problem.influenced_dofs

    @property
    def force(self):
        # Just an alias
        return self.forces

    __str__ = LinearPotentialFlowProblem.__str__
    __repr__ = LinearPotentialFlowProblem.__repr__
    _repr_pretty_ = LinearPotentialFlowProblem._repr_pretty_
    __rich_repr__ = LinearPotentialFlowProblem.__rich_repr__


class DiffractionResult(LinearPotentialFlowResult):

    def __init__(self, problem, *args, **kwargs):
        super().__init__(problem, *args, **kwargs)

    _str_other_attributes = DiffractionProblem._str_other_attributes
    _specific_rich_repr = DiffractionProblem._specific_rich_repr

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

    def __init__(self, problem, *args, **kwargs):
        super().__init__(problem, *args, **kwargs)
        self.radiating_dof = self.problem.radiating_dof

    _str_other_attributes = RadiationProblem._str_other_attributes
    _specific_rich_repr = RadiationProblem._specific_rich_repr

    @property
    def added_mass(self):
        return {dof: float(np.real(force)/(self.encounter_omega*self.encounter_omega)) for (dof, force) in self.forces.items()}

    @property
    def radiation_damping(self):
        if float(self.encounter_omega) in {0.0, np.inf} and self.forward_speed == 0.0:
            return {dof: 0.0 for dof in self.forces.keys()}
        else:
            return {dof: float(np.imag(force)/self.encounter_omega) for (dof, force) in self.forces.items()}

    # Aliases for backward compatibility
    added_masses = added_mass
    radiation_dampings = radiation_damping

    @property
    def records(self):
        params = self.problem._asdict()
        return [dict(params,
                     influenced_dof=dof,
                     added_mass=self.added_mass[dof],
                     radiation_damping=self.radiation_damping[dof])
                for dof in self.influenced_dofs]
