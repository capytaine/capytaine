# Copyright (C) 2017-2024 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>
"""Solver for the BEM problem.

.. code-block:: python

    problem = RadiationProblem(...)
    result = BEMSolver(green_functions=..., engine=...).solve(problem)

"""

import logging

import numpy as np

from datetime import datetime

from rich.progress import track

from capytaine.bem.problems_and_results import LinearPotentialFlowProblem
from capytaine.green_functions.delhommeau import Delhommeau
from capytaine.bem.engines import BasicMatrixEngine
from capytaine.io.xarray import problems_from_dataset, assemble_dataset, kochin_data_array
from capytaine.tools.optional_imports import silently_import_optional_dependency
from capytaine.tools.lists_of_points import _normalize_points, _normalize_free_surface_points
from capytaine.tools.symbolic_multiplication import supporting_symbolic_multiplication

LOG = logging.getLogger(__name__)

class BEMSolver:
    """
    Solver for linear potential flow problems.

    Parameters
    ----------
    green_function: AbstractGreenFunction, optional
        Object handling the computation of the Green function.
        (default: :class:`~capytaine.green_function.delhommeau.Delhommeau`)
    engine: MatrixEngine, optional
        Object handling the building of matrices and the resolution of linear systems with these matrices.
        (default: :class:`~capytaine.bem.engines.BasicMatrixEngine`)

    Attributes
    ----------
    exportable_settings : dict
        Settings of the solver that can be saved to reinit the same solver later.
    """

    def __init__(self, *, green_function=None, engine=None):
        self.green_function = Delhommeau() if green_function is None else green_function
        self.engine = BasicMatrixEngine() if engine is None else engine

        try:
            self.exportable_settings = {
                **self.green_function.exportable_settings,
                **self.engine.exportable_settings
            }
        except AttributeError:
            pass

    def __str__(self):
        return f"BEMSolver(engine={self.engine}, green_function={self.green_function})"

    def __repr__(self):
        return self.__str__()

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    @classmethod
    def from_exported_settings(settings):
        raise NotImplementedError

    def solve(self, problem, method='indirect', keep_details=True, _check_wavelength=True):
        """Solve the linear potential flow problem.

        Parameters
        ----------
        problem: LinearPotentialFlowProblem
            the problem to be solved
        method: string, optional
            select boundary integral approach indirect (i.e.Nemoh)/direct (i.e.WAMIT) (default: indirect)
        keep_details: bool, optional
            if True, store the sources and the potential on the floating body in the output object
            (default: True)
        _check_wavelength: bool, optional (default: True)
            If True, the frequencies are compared to the mesh resolution and
            the estimated first irregular frequency to warn the user.

        Returns
        -------
        LinearPotentialFlowResult
            an object storing the problem data and its results
        """
        LOG.info("Solve %s.", problem)

        if _check_wavelength:
            self._check_wavelength_and_mesh_resolution([problem])
            self._check_wavelength_and_irregular_frequencies([problem])

        if problem.forward_speed != 0.0:
            omega, wavenumber = problem.encounter_omega, problem.encounter_wavenumber
        else:
            omega, wavenumber = problem.omega, problem.wavenumber

        linear_solver = supporting_symbolic_multiplication(self.engine.linear_solver)
        if (method == 'direct'):
            if problem.forward_speed != 0.0:
                raise NotImplementedError("Direct solver is not able to solve problems with forward speed.")

            S, D = self.engine.build_matrices(
                    problem.body.mesh_including_lid, problem.body.mesh_including_lid,
                    problem.free_surface, problem.water_depth, wavenumber,
                    self.green_function, adjoint_double_layer=False
                    )

            potential = linear_solver(D, S @ problem.boundary_condition)
            pressure = 1j * omega * problem.rho * potential
            sources = None
        else:
            S, K = self.engine.build_matrices(
                    problem.body.mesh_including_lid, problem.body.mesh_including_lid,
                    problem.free_surface, problem.water_depth, wavenumber,
                    self.green_function, adjoint_double_layer=True
                    )

            sources = linear_solver(K, problem.boundary_condition)
            potential = S @ sources
            pressure = 1j * omega * problem.rho * potential
            if problem.forward_speed != 0.0:
                result = problem.make_results_container(sources=sources)
                # Temporary result object to compute the ∇Φ term
                nabla_phi = self._compute_potential_gradient(problem.body.mesh_including_lid, result)
                pressure += problem.rho * problem.forward_speed * nabla_phi[:, 0]

        pressure_on_hull = pressure[:problem.body.mesh.nb_faces]  # Discards pressure on lid if any
        forces = problem.body.integrate_pressure(pressure_on_hull)

        if not keep_details:
            result = problem.make_results_container(forces)
        else:
            result = problem.make_results_container(forces, sources, potential, pressure)

        LOG.debug("Done!")

        return result

    def solve_all(self, problems, *, method='indirect', n_jobs=1, progress_bar=True, _check_wavelength=True, **kwargs):
        """Solve several problems.
        Optional keyword arguments are passed to `BEMSolver.solve`.

        Parameters
        ----------
        problems: list of LinearPotentialFlowProblem
            several problems to be solved
        method: string, optional
            select boundary integral approach indirect (i.e.Nemoh)/direct (i.e.WAMIT) (default: indirect)
        n_jobs: int, optional (default: 1)
            the number of jobs to run in parallel using the optional dependency `joblib`
            By defaults: do not use joblib and solve sequentially.
        progress_bar: bool, optional (default: True)
            Display a progress bar while solving
        _check_wavelength: bool, optional (default: True)
            If True, the frequencies are compared to the mesh resolution and
            the estimated first irregular frequency to warn the user.

        Returns
        -------
        list of LinearPotentialFlowResult
            the solved problems
        """
        if _check_wavelength:
            self._check_wavelength_and_mesh_resolution(problems)
            self._check_wavelength_and_irregular_frequencies(problems)

        if n_jobs == 1:  # force sequential resolution
            problems = sorted(problems)
            if progress_bar:
                problems = track(problems, total=len(problems), description="Solving BEM problems")
            return [self.solve(pb, method=method, _check_wavelength=False, **kwargs) for pb in problems]
        else:
            joblib = silently_import_optional_dependency("joblib")
            if joblib is None:
                raise ImportError(f"Setting the `n_jobs` argument to {n_jobs} requires the missing optional dependency 'joblib'.")
            groups_of_problems = LinearPotentialFlowProblem._group_for_parallel_resolution(problems)
            parallel = joblib.Parallel(return_as="generator", n_jobs=n_jobs)
            groups_of_results = parallel(joblib.delayed(self.solve_all)(grp, method=method, n_jobs=1, progress_bar=False, _check_wavelength=False, **kwargs) for grp in groups_of_problems)
            if progress_bar:
                groups_of_results = track(groups_of_results,
                                          total=len(groups_of_problems),
                                          description=f"Solving BEM problems with {n_jobs} threads:")
            results = [res for grp in groups_of_results for res in grp]  # flatten the nested list
            return results

    @staticmethod
    def _check_wavelength_and_mesh_resolution(problems):
        """Display a warning if some of the problems have a mesh resolution
        that might not be sufficient for the given wavelength."""
        LOG.debug("Check wavelength with mesh resolution.")
        risky_problems = [pb for pb in problems
                          if 0.0 < pb.wavelength < pb.body.minimal_computable_wavelength]
        nb_risky_problems = len(risky_problems)
        if nb_risky_problems == 1:
            pb = risky_problems[0]
            freq_type = risky_problems[0].provided_freq_type
            freq = pb.__getattribute__(freq_type)
            LOG.warning(f"Mesh resolution for {pb}:\n"
                        f"The resolution of the mesh of the body {pb.body.__short_str__()} might "
                        f"be insufficient for {freq_type}={freq}.\n"
                        "This warning appears because the largest panel of this mesh "
                        f"has radius {pb.body.mesh.faces_radiuses.max():.3f} > wavelength/8."
                        )
        elif nb_risky_problems > 1:
            freq_type = risky_problems[0].provided_freq_type
            freqs = np.array([float(pb.__getattribute__(freq_type)) for pb in risky_problems])
            LOG.warning(f"Mesh resolution for {nb_risky_problems} problems:\n"
                        "The resolution of the mesh might be insufficient "
                        f"for {freq_type} ranging from {freqs.min():.3f} to {freqs.max():.3f}.\n"
                        "This warning appears when the largest panel of this mesh "
                        "has radius > wavelength/8."
                        )

    @staticmethod
    def _check_wavelength_and_irregular_frequencies(problems):
        """Display a warning if some of the problems might encounter irregular frequencies."""
        LOG.debug("Check wavelength with estimated irregular frequency.")
        risky_problems = [pb for pb in problems
                          if pb.body.first_irregular_frequency_estimate(g=pb.g) < pb.omega < np.inf]
        nb_risky_problems = len(risky_problems)
        if nb_risky_problems >= 1:
            if any(pb.body.lid_mesh is None for pb in problems):
                recommendation = "Setting a lid for the floating body is recommended."
            else:
                recommendation = "The lid might need to be closer to the free surface."
            if nb_risky_problems == 1:
                pb = risky_problems[0]
                freq_type = risky_problems[0].provided_freq_type
                freq = pb.__getattribute__(freq_type)
                LOG.warning(f"Irregular frequencies for {pb}:\n"
                            f"The body {pb.body.__short_str__()} might display irregular frequencies "
                            f"for {freq_type}={freq}.\n"
                            + recommendation
                            )
            elif nb_risky_problems > 1:
                freq_type = risky_problems[0].provided_freq_type
                freqs = np.array([float(pb.__getattribute__(freq_type)) for pb in risky_problems])
                LOG.warning(f"Irregular frequencies for {nb_risky_problems} problems:\n"
                            "Irregular frequencies might be encountered "
                            f"for {freq_type} ranging from {freqs.min():.3f} to {freqs.max():.3f}.\n"
                            + recommendation
                            )

    def fill_dataset(self, dataset, bodies, *, method='indirect', n_jobs=1, _check_wavelength=True, **kwargs):
        """Solve a set of problems defined by the coordinates of an xarray dataset.

        Parameters
        ----------
        dataset : xarray Dataset
            dataset containing the problems parameters: frequency, radiating_dof, water_depth, ...
        bodies : FloatingBody or list of FloatingBody
            The body or bodies involved in the problems
            They should all have different names.
        method: string, optional
            select boundary integral approach indirect (i.e.Nemoh)/direct (i.e.WAMIT) (default: indirect)
        n_jobs: int, optional (default: 1)
            the number of jobs to run in parallel using the optional dependency `joblib`
            By defaults: do not use joblib and solve sequentially.
        progress_bar: bool, optional (default: True)
            Display a progress bar while solving
        _check_wavelength: bool, optional (default: True)
            If True, the frequencies are compared to the mesh resolution and
            the estimated first irregular frequency to warn the user.

        Returns
        -------
        xarray Dataset
        """
        attrs = {'start_of_computation': datetime.now().isoformat(),
                 **self.exportable_settings}
        problems = problems_from_dataset(dataset, bodies)
        if 'theta' in dataset.coords:
            results = self.solve_all(problems, keep_details=True, method=method, n_jobs=n_jobs, _check_wavelength=_check_wavelength)
            kochin = kochin_data_array(results, dataset.coords['theta'])
            dataset = assemble_dataset(results, attrs=attrs, **kwargs)
            dataset.update(kochin)
        else:
            results = self.solve_all(problems, keep_details=False, method=method, n_jobs=n_jobs, _check_wavelength=_check_wavelength)
            dataset = assemble_dataset(results, attrs=attrs, **kwargs)
        return dataset


    def compute_potential(self, points, result):
        """Compute the value of the potential at given points for a previously solved potential flow problem.

        Parameters
        ----------
        points: array of shape (3,) or (N, 3), or 3-ple of arrays returned by meshgrid, or cpt.Mesh or cpt.CollectionOfMeshes object
            Coordinates of the point(s) at which the potential should be computed
        result: LinearPotentialFlowResult
            The return of the BEM solver

        Returns
        -------
        complex-valued array of shape (1,) or (N,) or (nx, ny, nz) or (mesh.nb_faces,) depending of the kind of input
            The value of the potential at the points

        Raises
        ------
        Exception: if the :code:`LinearPotentialFlowResult` object given as input does not contain the source distribution.
        """
        points, output_shape = _normalize_points(points, keep_mesh=True)
        if result.sources is None:
            raise Exception(f"""The values of the sources of {result} cannot been found.
            They probably have not been stored by the solver because the option keep_details=True have not been set or the direct method has been used.
            Please re-run the resolution with the indirect method and keep_details=True.""")

        S, _ = self.green_function.evaluate(points, result.body.mesh_including_lid, result.free_surface, result.water_depth, result.encounter_wavenumber)
        potential = S @ result.sources  # Sum the contributions of all panels in the mesh
        return potential.reshape(output_shape)

    def _compute_potential_gradient(self, points, result):
        points, output_shape = _normalize_points(points, keep_mesh=True)

        if result.sources is None:
            raise Exception(f"""The values of the sources of {result} cannot been found.
            They probably have not been stored by the solver because the option keep_details=True have not been set.
            Please re-run the resolution with this option.""")

        _, gradG = self.green_function.evaluate(points, result.body.mesh_including_lid, result.free_surface, result.water_depth, result.encounter_wavenumber,
                                                early_dot_product=False)
        velocities = np.einsum('ijk,j->ik', gradG, result.sources)  # Sum the contributions of all panels in the mesh
        return velocities.reshape((*output_shape, 3))

    def compute_velocity(self, points, result):
        """Compute the value of the velocity vector at given points for a previously solved potential flow problem.

        Parameters
        ----------
        points: array of shape (3,) or (N, 3), or 3-ple of arrays returned by meshgrid, or cpt.Mesh or cpt.CollectionOfMeshes object
            Coordinates of the point(s) at which the velocity should be computed
        result: LinearPotentialFlowResult
            The return of the BEM solver

        Returns
        -------
        complex-valued array of shape (3,) or (N,, 3) or (nx, ny, nz, 3) or (mesh.nb_faces, 3) depending of the kind of input
            The value of the velocity at the points

        Raises
        ------
        Exception: if the :code:`LinearPotentialFlowResult` object given as input does not contain the source distribution.
        """
        nabla_phi = self._compute_potential_gradient(points, result)
        if result.forward_speed != 0.0:
            nabla_phi[..., 0] -= result.forward_speed
        return nabla_phi

    def compute_pressure(self, points, result):
        """Compute the value of the pressure at given points for a previously solved potential flow problem.

        Parameters
        ----------
        points: array of shape (3,) or (N, 3), or 3-ple of arrays returned by meshgrid, or cpt.Mesh or cpt.CollectionOfMeshes object
            Coordinates of the point(s) at which the pressure should be computed
        result: LinearPotentialFlowResult
            The return of the BEM solver

        Returns
        -------
        complex-valued array of shape (1,) or (N,) or (nx, ny, nz) or (mesh.nb_faces,) depending of the kind of input
            The value of the pressure at the points

        Raises
        ------
        Exception: if the :code:`LinearPotentialFlowResult` object given as input does not contain the source distribution.
        """
        if result.forward_speed != 0:
            pressure = 1j * result.encounter_omega * result.rho * self.compute_potential(points, result)
            nabla_phi = self._compute_potential_gradient(points, result)
            pressure += result.rho * result.forward_speed * nabla_phi[..., 0]
        else:
            pressure = 1j * result.omega * result.rho * self.compute_potential(points, result)
        return pressure


    def compute_free_surface_elevation(self, points, result):
        """Compute the value of the free surface elevation at given points for a previously solved potential flow problem.

        Parameters
        ----------
        points: array of shape (2,) or (N, 2), or 2-ple of arrays returned by meshgrid, or cpt.Mesh or cpt.CollectionOfMeshes object
            Coordinates of the point(s) at which the free surface elevation should be computed
        result: LinearPotentialFlowResult
            The return of the BEM solver

        Returns
        -------
        complex-valued array of shape (1,) or (N,) or (nx, ny, nz) or (mesh.nb_faces,) depending of the kind of input
            The value of the free surface elevation at the points

        Raises
        ------
        Exception: if the :code:`LinearPotentialFlowResult` object given as input does not contain the source distribution.
        """
        points, output_shape = _normalize_free_surface_points(points, keep_mesh=True)

        if result.forward_speed != 0:
            fs_elevation = -1/result.g * (-1j*result.encounter_omega) * self.compute_potential(points, result)
            nabla_phi = self._compute_potential_gradient(points, result)
            fs_elevation += -1/result.g * result.forward_speed * nabla_phi[..., 0]
        else:
            fs_elevation = -1/result.g * (-1j*result.omega) * self.compute_potential(points, result)

        return fs_elevation.reshape(output_shape)


    ## Legacy

    def get_potential_on_mesh(self, result, mesh, chunk_size=50):
        """Compute the potential on a mesh for the potential field of a previously solved problem.
        Since the interaction matrix does not need to be computed in full to compute the matrix-vector product,
        only a few lines are evaluated at a time to reduce the memory cost of the operation.

        The newer method :code:`compute_potential` should be preferred in the future.

        Parameters
        ----------
        result : LinearPotentialFlowResult
            the return of the BEM solver
        mesh : Mesh or CollectionOfMeshes
            a mesh
        chunk_size: int, optional
            Number of lines to compute in the matrix.
            (legacy, should be passed as an engine setting instead).

        Returns
        -------
        array of shape (mesh.nb_faces,)
            potential on the faces of the mesh

        Raises
        ------
        Exception: if the :code:`Result` object given as input does not contain the source distribution.
        """
        LOG.info(f"Compute potential on {mesh.name} for {result}.")

        if result.sources is None:
            raise Exception(f"""The values of the sources of {result} cannot been found.
            They probably have not been stored by the solver because the option keep_details=True have not been set or the direct method has been used.
            Please re-run the resolution with the indirect method and keep_details=True.""")

        if chunk_size > mesh.nb_faces:
            S = self.engine.build_S_matrix(
                mesh,
                result.body.mesh_including_lid,
                result.free_surface, result.water_depth, result.wavenumber,
                self.green_function
            )
            phi = S @ result.sources

        else:
            phi = np.empty((mesh.nb_faces,), dtype=np.complex128)
            for i in range(0, mesh.nb_faces, chunk_size):
                faces_to_extract = list(range(i, min(i+chunk_size, mesh.nb_faces)))
                S = self.engine.build_S_matrix(
                    mesh.extract_faces(faces_to_extract),
                    result.body.mesh_including_lid,
                    result.free_surface, result.water_depth, result.wavenumber,
                    self.green_function
                )
                phi[i:i+chunk_size] = S @ result.sources

        LOG.debug(f"Done computing potential on {mesh.name} for {result}.")

        return phi

    def get_free_surface_elevation(self, result, free_surface, keep_details=False):
        """Compute the elevation of the free surface on a mesh for a previously solved problem.

        The newer method :code:`compute_free_surface_elevation` should be preferred in the future.

        Parameters
        ----------
        result : LinearPotentialFlowResult
            the return of the solver
        free_surface : FreeSurface
            a meshed free surface
        keep_details : bool, optional
            if True, keep the free surface elevation in the LinearPotentialFlowResult (default:False)

        Returns
        -------
        array of shape (free_surface.nb_faces,)
            the free surface elevation on each faces of the meshed free surface

        Raises
        ------
        Exception: if the :code:`Result` object given as input does not contain the source distribution.
        """
        if result.forward_speed != 0.0:
            raise NotImplementedError("For free surface elevation with forward speed, please use the `compute_free_surface_elevation` method.")

        fs_elevation = 1j*result.omega/result.g * self.get_potential_on_mesh(result, free_surface.mesh)
        if keep_details:
            result.fs_elevation[free_surface] = fs_elevation
        return fs_elevation
