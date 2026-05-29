# Copyright (C) 2017-2026 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>
"""Solver for the BEM problem.

.. code-block:: python

    problem = RadiationProblem(...)
    result = BEMSolver(engine=..., method=...).solve(problem)

"""

import os
import shutil
import textwrap
import logging
from datetime import datetime

import numpy as np
from rich.progress import track

from capytaine.bem.problems_and_results import LinearPotentialFlowProblem, DiffractionProblem
from capytaine.bem.engines import BasicMatrixEngine
from capytaine.bem.problems_checks import (
    _check_wavelength_and_mesh_resolution,
    _check_wavelength_and_water_depth,
    _check_wavelength_and_irregular_frequencies,
    _check_ram
)
from capytaine.io.xarray import problems_from_dataset, assemble_dataset, kochin_data_array
from capytaine.tools.memory_monitor import MemoryMonitor
from capytaine.tools.optional_imports import import_optional_dependency
from capytaine.tools.lists_of_points import _normalize_points, _normalize_free_surface_points
from capytaine.tools.symbolic_multiplication import supporting_symbolic_multiplication
from capytaine.tools.timer import Timer
from capytaine.ui.env_vars import look_for_boolean_var
from capytaine.ui.error_messages import display_grouped_errors

LOG = logging.getLogger(__name__)

# Mapping between a dtype and its complex version
COMPLEX_DTYPE = {
                    np.float32: np.complex64,
                    np.float64: np.complex128,
                    np.complex64 : np.complex64,
                    np.complex128 : np.complex128
                }

class BEMSolver:
    """
    Solver for linear potential flow problems.

    Parameters
    ----------
    engine: MatrixEngine, optional
        Object handling the building of matrices and the resolution of linear systems with these matrices.
        (default: :class:`~capytaine.bem.engines.BasicMatrixEngine`)
    method: string, optional
        select boundary integral equation used to solve the problems.
        Accepted values: "indirect" (as in e.g. Nemoh), "direct" (as in e.g. WAMIT)
        Default value: "indirect"
    green_function: AbstractGreenFunction, optional
        For convenience and backward compatibility, the Green function can be
        set here if the engine is the default one.
        This argument is just passed to the default engine at initialization.

    Attributes
    ----------
    timer: Timer
        Storing the time spent on each subtasks of the resolution
    exportable_settings : dict
        Settings of the solver that can be saved to reinit the same solver later.
    """

    def __init__(self, *, green_function=None, engine=None, method="indirect"):

        if engine is None:
            self.engine = BasicMatrixEngine(green_function=green_function)
        else:
            if green_function is not None:
                raise ValueError("If you are not using the default engine, set the Green function in the engine.\n"
                                 "Setting the Green function in the solver is only a shortcut to set up "
                                 "the Green function of the default engine since Capytaine version 3.0")
            self.engine = engine

        if method.lower() not in {"direct", "indirect"}:
            raise ValueError(f"Unrecognized method when initializing solver: {repr(method)}. Expected \"direct\" or \"indirect\".")
        self.method = method.lower()

        self.reset_timer()

        self.exportable_settings = {
            **self.engine.exportable_settings,
            "method": self.method,
        }

    def __str__(self):
        return f"BEMSolver(engine={self.engine}, method={self.method})"

    def __repr__(self):
        return self.__str__()

    def reset_timer(self):
        self.timer = Timer(default_tags={"process": 0})

    def timer_summary(self):
        df = self.timer.as_dataframe()
        df["step"] = df["step"].where(
                df["step"].str.startswith("Total"), "  " + df["step"]
                )
        total = (
                df.groupby(["step", "process"])
                ["timing"].sum()
                .unstack()
                )
        return total

    def displayed_total_summary(self, width=None):
        total = self.timer_summary()
        if width is None:
            width = shutil.get_terminal_size(fallback=(80, 20)).columns - 25
        total_str = total.to_string(
                float_format="{:.2f}".format,
                line_width=width,
                )
        return textwrap.indent(total_str, "  ")

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    @classmethod
    def from_exported_settings(settings):
        raise NotImplementedError

    def _solve(self, problem, method=None, keep_details=True, _check_wavelength=True):
        """Called by BEMSolver.solve. See the documentation therein."""
        LOG.info("Solve %s.", problem)

        if _check_wavelength:
            _check_wavelength_and_mesh_resolution([problem])
            _check_wavelength_and_water_depth([problem])
            _check_wavelength_and_irregular_frequencies([problem])

        if isinstance(problem, DiffractionProblem) and float(problem.encounter_omega) in {0.0, np.inf}:
            raise ValueError("Diffraction problems at zero or infinite frequency are not defined")
            # This error used to be raised when initializing the problem.
            # It is now raised here, in order to be catchable by
            # _solve_and_catch_errors, such that batch resolution
            # can include this kind of problems without the full batch
            # failing.
            # Note that if this error was not raised here, the resolution
            # would still fail with a less explicit error message.

        if problem.forward_speed != 0.0:
            omega, wavenumber = problem.encounter_omega, problem.encounter_wavenumber
        else:
            omega, wavenumber = problem.omega, problem.wavenumber
        gf_params = dict(free_surface=problem.free_surface, water_depth=problem.water_depth, wavenumber=wavenumber)

        linear_solver = supporting_symbolic_multiplication(self.engine.linear_solver)
        method = method if method is not None else self.method
        if (method == 'direct'):
            if problem.forward_speed != 0.0:
                raise NotImplementedError("Direct solver is not able to solve problems with forward speed.")

            with self.timer(step="Green function"):
                S, D = self.engine.build_matrices(
                        problem.body.mesh_including_lid, problem.body.mesh_including_lid,
                        **gf_params, adjoint_double_layer=False, diagonal_term_in_double_layer=True,
                        )
            with self.timer(step="Matrix-vector product"):
                rhs = S @ problem.boundary_condition
            with self.timer(step="Linear solver"):
                rhs = rhs.astype(COMPLEX_DTYPE[D.dtype.type])
                potential = linear_solver(D, rhs)
            pressure = 1j * omega * problem.rho * potential
            sources = None
        else:
            with self.timer(step="Green function"):
                S, K = self.engine.build_matrices(
                        problem.body.mesh_including_lid, problem.body.mesh_including_lid,
                        **gf_params, adjoint_double_layer=True, diagonal_term_in_double_layer=True,
                        )
            with self.timer(step="Linear solver"):
                rhs = problem.boundary_condition.astype(COMPLEX_DTYPE[K.dtype.type])
                sources = linear_solver(K, rhs)
            with self.timer(step="Matrix-vector product"):
                potential = S @ sources
            pressure = 1j * omega * problem.rho * potential
            if problem.forward_speed != 0.0:
                result = problem.make_results_container(sources=sources)
                # Temporary result object to compute the ∇Φ term
                nabla_phi = self._compute_potential_gradient(problem.body.mesh_including_lid, result)
                pressure += problem.rho * problem.forward_speed * nabla_phi[:, 0]

        pressure_on_hull = pressure[problem.body.hull_mask]  # Discards pressure on lid if any
        forces = problem.body.integrate_pressure(pressure_on_hull)
        if not keep_details:
            result = problem.make_results_container(forces)
        else:
            result = problem.make_results_container(forces, sources, potential, pressure)

        LOG.debug("Done!")

        return result

    def solve(self, problem, method=None, keep_details=True, n_threads=None, _check_wavelength=True):
        """Solve the linear potential flow problem.

        Parameters
        ----------
        problem: LinearPotentialFlowProblem
            the problem to be solved
        keep_details: bool, optional
            if True, store the sources and the potential on the floating body in the output object
            (default: True)
        n_threads: int, optional
            the number of threads to use for the resolution.
            Requires the optional package `threadpoolctl` to set.
        _check_wavelength: bool, optional (default: True)
            If True, the frequencies are compared to the mesh resolution and
            the estimated first irregular frequency to warn the user.
        method: str, optional
            select boundary integral equation used to solve the problems.
            It is recommended to set the method more globally when initializing the solver.
            If provided here, the value in argument of `solve` overrides the global one.

        Returns
        -------
        LinearPotentialFlowResult
            an object storing the problem data and its results
        """
        # Thin wrapper around _solve adding the timer and the threading control
        with self.timer(step="Total solve function"):
            if n_threads is None:
                return self._solve(problem, method=method, keep_details=keep_details, _check_wavelength=_check_wavelength)
            else:
                threadpoolctl = import_optional_dependency(
                        "threadpoolctl",
                        error_message=f"Setting the `n_threads` argument to {n_threads} with `n_jobs=1` requires the missing optional dependency 'threadpoolctl'."
                        )
                with threadpoolctl.threadpool_limits(limits=n_threads):
                    return self._solve(problem, method=method, keep_details=keep_details, _check_wavelength=_check_wavelength)


    def _solve_and_catch_errors(self, problem, *args, _display_errors, **kwargs):
        """Same as BEMSolver.solve() but returns a
        FailedLinearPotentialFlowResult when the resolution failed."""
        try:
            res = self.solve(problem, *args, **kwargs)
        except Exception as e:
            res = problem.make_failed_results_container(e)
            if _display_errors:
                display_grouped_errors([res])
        return res


    def solve_all(self, problems, *, method=None, n_jobs=1, n_threads=None, progress_bar=None, _check_wavelength=True, _display_errors=True, **kwargs):
        """Solve several problems.
        Optional keyword arguments are passed to `BEMSolver.solve`.

        Parameters
        ----------
        problems: list of LinearPotentialFlowProblem
            several problems to be solved
        method: string, optional
            select boundary integral equation used to solve the problems.
            It is recommended to set the method more globally when initializing the solver.
            If provided here, the value in argument of `solve_all` overrides the global one.
        n_jobs: int, optional (default: 1)
            the number of jobs to run in parallel using the optional dependency ``joblib``
            By defaults: do not use joblib and solve sequentially.
        n_threads: int, optional
            the number of threads used to solve each problem.
            The total number of used CPU will be n_jobs×n_threads.
            By default: use as much as possible.
            Requires the optional dependency ``threadpoolctl`` if ``n_jobs==1``.
            Also controlled by the environment variables ``OMP_NUM_THREADS`` and ``MKL_NUM_THREADS``.
        progress_bar: bool, optional
            Display a progress bar while solving.
            If no value is provided to this method directly,
            check whether the environment variable `CAPYTAINE_PROGRESS_BAR` is defined
            and otherwise default to True.
        _check_wavelength: bool, optional (default: True)
            If True, the frequencies are compared to the mesh resolution and
            the estimated first irregular frequency to warn the user.

        Returns
        -------
        list of LinearPotentialFlowResult
            the solved problems
        """
        if _check_wavelength:
            _check_wavelength_and_mesh_resolution(problems)
            _check_wavelength_and_water_depth(problems)
            _check_wavelength_and_irregular_frequencies(problems)

        _check_ram(problems, self.engine, n_jobs)

        if progress_bar is None:
            progress_bar = look_for_boolean_var("CAPYTAINE_PROGRESS_BAR", default=True)

        groups_of_problems = LinearPotentialFlowProblem._group_for_parallel_resolution(problems)  # or not parallel resolution

        with MemoryMonitor():
            if n_jobs == 1:  # force sequential resolution

                def solve_single_pb(pb):
                    return self._solve_and_catch_errors(pb, method=method, _display_errors=True,
                                                        _check_wavelength=False, n_threads=n_threads, **kwargs)

                if progress_bar:
                    groups_of_problems = track(groups_of_problems, total=len(groups_of_problems), description="Solving BEM problems")

                results = [solve_single_pb(pb) for group in groups_of_problems for pb in group]

            else:
                joblib = import_optional_dependency(
                    "joblib",
                    error_message=f"Setting the `n_jobs` argument to {n_jobs} requires the missing optional dependency 'joblib'."
                )

                def solver_group_of_problems_in_other_process(group):
                    # Meant to be called in another process by joblib's backend.

                    self.reset_timer()
                    # Timer data will be concatenated to the Timer of the main process.
                    # We reset the timer at each call to avoid concatenating
                    # the same data twice in the main process.

                    group_results = [self._solve_and_catch_errors(
                        pb,
                        method=method,
                        _display_errors=False,  # Will be displayed by the main process, see below
                        _check_wavelength=False,  # Already checked, see above
                        n_threads=None,  # Should be controlled by joblib outside of this function
                        **kwargs
                    ) for pb in group]

                    return group_results, self.timer, os.getpid()

                with joblib.parallel_config(backend='loky', inner_max_num_threads=n_threads):
                    parallel = joblib.Parallel(return_as="generator", n_jobs=n_jobs)
                    groups_of_results = parallel(joblib.delayed(solver_group_of_problems_in_other_process)(group) for group in groups_of_problems)

                if progress_bar:
                    # The progress bar is on the results iterator, because the inputs are consumed immediately by joblib
                    groups_of_results = track(groups_of_results,
                                              total=len(groups_of_problems),
                                              description=f"Solving BEM problems with {n_jobs} processes:")

                results = []
                process_id_mapping = {}
                for grp_results, other_timer, process_id in groups_of_results:
                    results.extend(grp_results)
                    display_grouped_errors(grp_results)
                    if process_id not in process_id_mapping:
                        process_id_mapping[process_id] = len(process_id_mapping) + 1
                    self.timer.add_data_from_other_timer(other_timer, process=process_id_mapping[process_id])

        LOG.info("Solver timer summary (in seconds):\n%s", self.displayed_total_summary())
        return results


    def fill_dataset(self, dataset, bodies, *, method=None, n_jobs=1, n_threads=None, _check_wavelength=True, progress_bar=None, **kwargs):
        """Solve a set of problems defined by the coordinates of an xarray dataset.

        Parameters
        ----------
        dataset : xarray Dataset
            dataset containing the problems parameters: frequency, radiating_dof, water_depth, ...
        bodies : FloatingBody or Multibody or list of FloatingBody or list of Multibody
            The body or bodies involved in the problems
            They should all have different names.
        method: string, optional
            select boundary integral equation used to solve the problems.
            It is recommended to set the method more globally when initializing the solver.
            If provided here, the value in argument of `fill_dataset` overrides the global one.
        n_jobs: int, optional (default: 1)
            the number of jobs to run in parallel using the optional dependency ``joblib``.
            By defaults: do not use joblib and solve sequentially.
        n_threads: int, optional
            the number of threads used to solve each problem.
            The total number of used CPU will be n_jobs×n_threads.
            By default: use as much as possible.
            Requires the optional dependency ``threadpoolctl`` if ``n_jobs==1``.
            Also controlled by the environment variables ``OMP_NUM_THREADS`` and ``MKL_NUM_THREADS``.
        progress_bar: bool, optional
            Display a progress bar while solving.
            If no value is provided to this method directly,
            check whether the environment variable `CAPYTAINE_PROGRESS_BAR` is defined
            and otherwise default to True.
        _check_wavelength: bool, optional (default: True)
            If True, the frequencies are compared to the mesh resolution and
            the estimated first irregular frequency to warn the user.

        Returns
        -------
        xarray Dataset
        """
        attrs = {'start_of_computation': datetime.now().isoformat(),
                 **self.exportable_settings}
        if method is not None:  # Overrides the method in self.exportable_settings
            attrs["method"] = method
        problems = problems_from_dataset(dataset, bodies)
        if 'theta' in dataset.coords:
            results = self.solve_all(problems, keep_details=True, method=method, n_jobs=n_jobs, n_threads=n_threads, _check_wavelength=_check_wavelength, progress_bar=progress_bar)
            kochin = kochin_data_array(results, dataset.coords['theta'])
            dataset = assemble_dataset(results, attrs=attrs, **kwargs)
            dataset.update(kochin)
        else:
            results = self.solve_all(problems, keep_details=False, method=method, n_jobs=n_jobs, n_threads=n_threads, _check_wavelength=_check_wavelength, progress_bar=progress_bar)
            dataset = assemble_dataset(results, attrs=attrs, **kwargs)
        return dataset


    def compute_potential(self, points, result):
        """Compute the value of the potential at given points for a previously solved potential flow problem.

        Parameters
        ----------
        points: array of shape (3,) or (N, 3), or 3-ple of arrays returned by meshgrid, or MeshLike object
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
        gf_params = dict(free_surface=result.free_surface, water_depth=result.water_depth, wavenumber=result.encounter_wavenumber)

        points, output_shape = _normalize_points(points)
        if result.sources is None:
            raise Exception(f"""The values of the sources of {result} cannot been found.
            They probably have not been stored by the solver because the option keep_details=True have not been set or the direct method has been used.
            Please re-run the resolution with the indirect method and keep_details=True.""")

        with self.timer(step="Post-processing potential"):
            S = self.engine.build_S_matrix(points, result.body.mesh_including_lid, **gf_params)
            potential = S @ result.sources  # Sum the contributions of all panels in the mesh
        return potential.reshape(output_shape)

    def _compute_potential_gradient(self, points, result):
        points, output_shape = _normalize_points(points, keep_mesh=True)
        # keep_mesh, because we need the normal vectors associated with each collocation points to compute the fullK matrix

        if result.sources is None:
            raise Exception(f"""The values of the sources of {result} cannot been found.
            They probably have not been stored by the solver because the option keep_details=True have not been set.
            Please re-run the resolution with this option.""")

        gf_params = dict(free_surface=result.free_surface, water_depth=result.water_depth, wavenumber=result.encounter_wavenumber)
        with self.timer(step="Post-processing velocity"):
            gradG = self.engine.build_fullK_matrix(points, result.body.mesh_including_lid, **gf_params)
            # gradG is either:
            # - an array of shape (3, nb_points, mesh_including_lid.nb_faces)
            # - a 3-ple of matrices of the shape (nb_points, mesh_including_lid), that could be stored as LazyMatrix.
            vx = gradG[0] @ result.sources
            vy = gradG[1] @ result.sources
            vz = gradG[2] @ result.sources
            # The matrix-vector product here computes the integral over the mesh of the contributions of each panel in the mesh
            velocities = np.stack([vx, vy, vz], axis=-1)
            # velocities.shape = (nb_points, 3)
        return velocities.reshape((*output_shape, 3))

    def compute_velocity(self, points, result):
        """Compute the value of the velocity vector at given points for a previously solved potential flow problem.

        Parameters
        ----------
        points: array of shape (3,) or (N, 3), or 3-ple of arrays returned by meshgrid, or MeshLike object
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
        points: array of shape (3,) or (N, 3), or 3-ple of arrays returned by meshgrid, or MeshLike object
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
        points: array of shape (2,) or (N, 2), or 2-ple of arrays returned by meshgrid, or MeshLike object
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
        points, output_shape = _normalize_free_surface_points(points)

        if result.forward_speed != 0:
            fs_elevation = -1/result.g * (-1j*result.encounter_omega) * self.compute_potential(points, result)
            nabla_phi = self._compute_potential_gradient(points, result)
            fs_elevation += -1/result.g * result.forward_speed * nabla_phi[..., 0]
        else:
            fs_elevation = -1/result.g * (-1j*result.omega) * self.compute_potential(points, result)

        return fs_elevation.reshape(output_shape)
