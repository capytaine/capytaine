#!/usr/bin/env python
# coding: utf-8
"""Solver for the BEM problem based on Nemoh's Green function.

Example
-------

::

    problem = RadiationProblem(...)
    result = Nemoh().solve(problem)

"""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from datetime import datetime
from functools import lru_cache

import numpy as np
import xarray as xr

from capytaine.matrices import linear_solvers
from capytaine.matrices.builders import identity_like
from capytaine.bem.hierarchical_toeplitz_matrices import hierarchical_toeplitz_matrices
from capytaine.bem.prony_decomposition import find_best_exponential_decomposition
import capytaine.bem.NemohCore as NemohCore
from capytaine.io.xarray import problems_from_dataset, assemble_dataset, kochin_data_array


LOG = logging.getLogger(__name__)

# Caching the tabulation of the integrals for Delhommeau's Green function.
tabulated_integrals = lru_cache(maxsize=1)(NemohCore.initialize_green_wave.initialize_tabulated_integrals)


class Nemoh:
    """Solver for the BEM problem based on Nemoh's Green function.

    Parameters
    ----------
    tabulation_nb_integration_points: int, optional
        Number of points for the evaluation of the tabulated elementary integrals w.r.t. :math:`theta`
        used for the computation of the Green function (default: 251)
    finite_depth_prony_decomposition_method: string, optional
        The implementation of the Prony decomposition used to compute the finite depth Green function.
    use_symmetries: bool, optional
        if True, use the symmetries of the meshes when computing matrices and solving linear system
    ACA_distance: float, optional
        Above this distance, the ACA is used to approximate the matrix with a low-rank block.
    ACA_tol: float, optional
        The tolerance of the ACA when building a low-rank matrix.
    cache_rankine_matrices: bool, optional
        If True, cache intermediate matrices when computing the influence matrices,
        if False, cache only the final result.
    matrix_cache_size: int, optional
        number of matrices to keep in cache
    linear_solver: str or function, optional
        Setting of the numerical solver for linear problems Ax = b.
        It can be set with the name of a preexisting solver (available: "direct" [default], "gmres", "store_lu")
        or by passing directly a solver function.

    Attributes
    ----------
    tabulated_integrals: 3-ple of arrays
        Tabulated integrals for the computation of the Green function.
    """

    defaults_settings = dict(
        tabulation_nb_integration_points=251,
        finite_depth_prony_decomposition_method='fortran',
        linear_solver='gmres',
        hierarchical_matrices=True,
        ACA_distance=np.infty,
        ACA_tol=1e-2,
        matrix_cache_size=1,
        cache_rankine_matrices=False,
    )

    available_linear_solvers = {'direct': linear_solvers.solve_directly,
                                'gmres': linear_solvers.solve_gmres}

    def __init__(self, **settings):

        # Check that all the given settings are relevant.
        for key in settings:
            if key not in self.defaults_settings:
                LOG.warning(f"Unrecognized solver option: {key}")

        settings = {**self.defaults_settings, **settings}  # Complete settings with default settings

        self.settings = settings  # Keep a copy for saving in output dataset

        LOG.info("Initialize Nemoh's Green function.")
        self.tabulated_integrals = tabulated_integrals(328, 46, settings['tabulation_nb_integration_points'])

        if settings['linear_solver'] in Nemoh.available_linear_solvers:
            self.linear_solver = Nemoh.available_linear_solvers[settings['linear_solver']]
        else:
            self.linear_solver = settings['linear_solver']

        if settings['matrix_cache_size'] > 0 and settings['cache_rankine_matrices']:
            if settings['hierarchical_matrices']:
                # If the rankine matrix is cached, the recursive decomposition of the matrix
                # has to be done before the caching.
                # Otherwise, only the latest blocks of the matrices would be cached,
                # which is not the expected behavior.
                # It follows that the decomposition has to be done twice:
                # once for the Rankine part and once for the wave part.
                self.build_matrices_rankine = hierarchical_toeplitz_matrices(
                    self.build_matrices_rankine,
                    ACA_tol=settings['ACA_tol'],
                    ACA_distance=settings['ACA_distance'],
                    dtype=np.float64
                )
                self.build_matrices_wave = hierarchical_toeplitz_matrices(
                    self.build_matrices_wave,
                    ACA_tol=settings['ACA_tol'],
                    ACA_distance=settings['ACA_distance'],
                    dtype=np.complex128
                )

            self.build_matrices_rankine = lru_cache(maxsize=settings['matrix_cache_size'])(self.build_matrices_rankine)
            self.build_matrices = lru_cache(maxsize=settings['matrix_cache_size'])(self.build_matrices)

        else:
            if settings['hierarchical_matrices']:
                # If the rankine matrix is not cached, the recursive decomposition of the matrix
                # can be done at the top level.
                self.build_matrices = hierarchical_toeplitz_matrices(
                    self.build_matrices,
                    ACA_tol=settings['ACA_tol'],
                    ACA_distance=settings['ACA_distance'],
                    dtype=np.complex128
                )
            if settings['matrix_cache_size'] > 0:
                self.build_matrices = lru_cache(maxsize=settings['matrix_cache_size'])(self.build_matrices)

    def exportable_settings(self):
        settings = self.settings.copy()
        if not settings['hierarchical_matrices']:
            del settings['ACA_distance']
            del settings['ACA_tol']
        if settings['matrix_cache_size'] == 0:
            del settings['cache_rankine_matrices']
        settings['linear_solver'] = str(settings['linear_solver'])
        return settings

    def solve(self, problem, keep_details=True):
        """Solve the BEM problem using Nemoh.

        Parameters
        ----------
        problem: LinearPotentialFlowProblem
            the problem to be solved
        keep_details: bool, optional
            if True, store the sources and the potential on the floating body in the output object
            (default: True)

        Returns
        -------
        LinearPotentialFlowResult
            an object storing the problem data and its results
        """

        LOG.info("Solve %s.", problem)

        if problem.wavelength < 8*problem.body.mesh.faces_radiuses.max():
            LOG.warning(f"Resolution of the mesh (8×max_radius={8*problem.body.mesh.faces_radiuses.max():.2e}) "
                        f"might be insufficient for this wavelength (wavelength={problem.wavelength:.2e})!")

        S, K = self.build_matrices(
            problem.body.mesh, problem.body.mesh,
            free_surface=problem.free_surface, sea_bottom=problem.sea_bottom, wavenumber=problem.wavenumber
        )

        sources = self.linear_solver(K, problem.boundary_condition)
        potential = S @ sources

        result = problem.make_results_container()
        if keep_details:
            result.sources = sources
            result.potential = potential

        for influenced_dof_name, influenced_dof_vectors in problem.influenced_dofs.items():
            # Scalar product on each face:
            influenced_dof_normal = np.sum(influenced_dof_vectors * problem.body.mesh.faces_normals, axis=1)
            # Sum over all faces:
            integrated_potential = - problem.rho * np.sum(potential * influenced_dof_normal * problem.body.mesh.faces_areas)
            # Store result:
            result.store_force(influenced_dof_name, integrated_potential)
            # Depending of the type of problem, the force will be kept as a complex-valued Froude-Krylov force
            # or stored as a couple of added mass and radiation damping coefficients.

        LOG.debug("Done!")

        return result

    def solve_all(self, problems, **kwargs):
        """Solve several problems.
        Optional keyword arguments are passed to `Nemoh.solve`.

        Parameters
        ----------
        problems: list of LinearPotentialFlowProblem
            several problems to be solved

        Returns
        -------
        list of LinearPotentialFlowResult
            the solved problems
        """
        return [self.solve(problem, **kwargs) for problem in sorted(problems)]

    def fill_dataset(self, dataset, bodies, **kwargs):
        """Solve a set of problems defined by the coordinates of an xarray dataset.

        Parameters
        ----------
        dataset : xarray Dataset
            dataset containing the problems parameters: frequency, radiating_dof, water_depth, ...
        bodies : list of FloatingBody
            the bodies involved in the problems

        Returns
        -------
        xarray Dataset
        """
        attrs = {'start_of_computation': datetime.now().isoformat(),
                 **self.exportable_settings()}
        problems = problems_from_dataset(dataset, bodies)
        if 'theta' in dataset.coords:
            results = self.solve_all(problems, keep_details=True)
            kochin = kochin_data_array(results, dataset.coords['theta'])
            dataset = assemble_dataset(results, attrs=attrs, **kwargs)
            dataset['kochin'] = kochin
        else:
            results = self.solve_all(problems, keep_details=False)
            dataset = assemble_dataset(results, attrs=attrs, **kwargs)
        return dataset

    #######################
    #  Building matrices  #
    #######################

    def build_matrices(self, mesh1, mesh2, free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0):
        r"""Build the S and K influence matrices between mesh1 and mesh2.
        In brief, it calls `build_matrices_rankine` and `build_matrices_wave` and sum their outputs.
        It also adds :math:`\mathbb{I}/2` matrix to :math:`V` to get :math:`K`.

        Parameters
        ----------
        mesh1: Mesh or CollectionOfMeshes
            mesh of the receiving body (where the potential is measured)
        mesh2: Mesh or CollectionOfMeshes
            mesh of the source body (over which the source distribution is integrated)
        free_surface: float, optional
            position of the free surface (default: :math:`z = 0`)
        sea_bottom: float, optional
            position of the sea bottom (default: :math:`z = -\infty`)
        wavenumber: float, optional
            wavenumber (default: 1.0)

        Returns
        -------
        couple of matrix-like objects (either 2D arrays or BlockMatrix objects)
            couple of influence matrices
        """

        LOG.debug(f"\tEvaluating matrix of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name} "
                  f"for depth={free_surface-sea_bottom} and wavenumber={wavenumber}.")

        Srankine, Vrankine = self.build_matrices_rankine(mesh1, mesh2, free_surface, sea_bottom, wavenumber)

        if (free_surface == np.infty or
                (free_surface - sea_bottom == np.infty and wavenumber in (0, np.infty))):
            # No more terms in the Green function
            return Srankine, Vrankine + identity_like(Vrankine)/2

        Swave, Vwave = self.build_matrices_wave(mesh1, mesh2, free_surface, sea_bottom, wavenumber)

        # The real valued matrices Srankine and Vrankine are automatically recasted as complex in the sum.
        Swave += Srankine
        if mesh1 is mesh2:
            Vwave += Vrankine + identity_like(Vrankine)/2
        else:
            Vwave += Vrankine
        return Swave, Vwave

    def build_matrices_rankine(self, mesh1, mesh2, free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0):
        r"""Build the Rankine part of the S and V influence matrices between mesh1 and mesh2.

        Parameters
        ----------
        mesh1: Mesh or CollectionOfMeshes
            mesh of the receiving body (where the potential is measured)
        mesh2: Mesh or CollectionOfMeshes
            mesh of the source body (over which the source distribution is integrated)
        free_surface: float, optional
            position of the free surface (default: :math:`z = 0`)
        sea_bottom: float, optional
            position of the sea bottom (default: :math:`z = -\infty`)
        wavenumber: float, optional
            wavenumber (default: 1.0)

        Returns
        -------
        couple of real-valued matrix-like objects (either 2D arrays or BlockMatrix objects)
            couple of influence matrices
        """
        # RANKINE TERM

        S, V = NemohCore.green_rankine.build_matrices_rankine_source(
            mesh1.faces_centers, mesh1.faces_normals,
            mesh2.vertices,      mesh2.faces + 1,
            mesh2.faces_centers, mesh2.faces_normals,
            mesh2.faces_areas,   mesh2.faces_radiuses,
                                 )

        if free_surface == np.infty:
            # No free surface, no more terms in the Green function
            return S, V

        # REFLECTION TERM

        def reflect_vector(x):
            y = x.copy()
            y[:, 2] *= -1
            return y

        if free_surface - sea_bottom == np.infty:
            # INFINITE DEPTH
            def reflect_point(x):
                y = x.copy()
                # y[:, 2] = 2*free_surface - x[:, 2]
                y[:, 2] *= -1
                y[:, 2] += 2*free_surface
                return y
        else:
            # FINITE DEPTH
            def reflect_point(x):
                y = x.copy()
                # y[:, 2] = 2*sea_bottom - x[:, 2]
                y[:, 2] *= -1
                y[:, 2] += 2*sea_bottom
                return y

        Srefl, Vrefl = NemohCore.green_rankine.build_matrices_rankine_source(
            reflect_point(mesh1.faces_centers), reflect_vector(mesh1.faces_normals),
            mesh2.vertices,      mesh2.faces + 1,
            mesh2.faces_centers, mesh2.faces_normals,
            mesh2.faces_areas,   mesh2.faces_radiuses,
                                 )

        if free_surface - sea_bottom < np.infty or wavenumber == 0.0:
            S += Srefl
            V += Vrefl
        else:
            S -= Srefl
            V -= Vrefl

        return S, V

    def build_matrices_wave(self, mesh1, mesh2, free_surface, sea_bottom, wavenumber):
        r"""Build the wave part of the influence matrices between mesh1 and mesh2.

        Parameters
        ----------
        mesh1: Mesh or CollectionOfMeshes
            mesh of the receiving body (where the potential is measured)
        mesh2: Mesh or CollectionOfMeshes
            mesh of the source body (over which the source distribution is integrated)
        free_surface: float, optional
            position of the free surface (default: :math:`z = 0`)
        sea_bottom: float, optional
            position of the sea bottom (default: :math:`z = -\infty`)
        wavenumber: float, optional
            wavenumber (default: 1.0)

        Returns
        -------
        couple of complex-valued matrix-like objects (either 2D arrays or BlockMatrix objects)
            couple of influence matrices
        """
        depth = free_surface - sea_bottom
        if depth == np.infty:
            return NemohCore.green_wave.build_matrices_wave_source(
                mesh1.faces_centers, mesh1.faces_normals,
                mesh2.faces_centers, mesh2.faces_areas,
                wavenumber, 0.0,
                *self.tabulated_integrals,
                np.empty(1), np.empty(1),  # Dummy arrays that won't actually be used by the fortran code.
                mesh1 is mesh2
            )
        else:
            a_exp, lamda_exp = find_best_exponential_decomposition(
                wavenumber*depth*np.tanh(wavenumber*depth),
                wavenumber*depth,
                method=self.settings['finite_depth_prony_decomposition_method'],
            )

            return NemohCore.green_wave.build_matrices_wave_source(
                mesh1.faces_centers, mesh1.faces_normals,
                mesh2.faces_centers, mesh2.faces_areas,
                wavenumber, depth,
                *self.tabulated_integrals,
                lamda_exp, a_exp,
                mesh1 is mesh2
                )

    #######################
    #  Compute potential  #
    #######################

    def get_potential_on_mesh(self, result, mesh, chunk_size=50):
        """Compute the potential on a mesh for the potential field of a previously solved problem.
        Since the interaction matrix does not need to be computed in full to compute the matrix-vector product,
        only a few lines are evaluated at a time to reduce the memory cost of the operation.

        Parameters
        ----------
        result : LinearPotentialFlowResult
            the return of Nemoh's solver
        mesh : Mesh or CollectionOfMeshes
            a mesh
        chunk_size: int, optional
            Number of lines to compute in the matrix.

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
            They probably have not been stored by the solver because the option keep_details=True have not been set.
            Please re-run the resolution with this option.""")

        if chunk_size > mesh.nb_faces:
            S, _ = self.build_matrices(
                mesh,
                result.body.mesh,
                free_surface=result.free_surface,
                sea_bottom=result.sea_bottom,
                wavenumber=result.wavenumber
            )
            phi = S @ result.sources

        else:
            phi = np.empty((mesh.nb_faces,), dtype=np.complex128)

            for i in range(0, mesh.nb_faces, chunk_size):
                S, _ = self.build_matrices(
                    mesh.extract_faces(list(range(i, i+chunk_size))),
                    result.body.mesh,
                    free_surface=result.free_surface,
                    sea_bottom=result.sea_bottom,
                    wavenumber=result.wavenumber
                )
                phi[i:i+chunk_size] = S @ result.sources

        LOG.debug(f"Done computing potential on {mesh.name} for {result}.")

        return phi

    def get_free_surface_elevation(self, result, free_surface, keep_details=False):
        """Compute the elevation of the free surface on a mesh for a previously solved problem.

        Parameters
        ----------
        result : LinearPotentialFlowResult
            the return of Nemoh's solver
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
        fs_elevation = 1j*result.omega/result.g * self.get_potential_on_mesh(result, free_surface.mesh)
        if keep_details:
            result.fs_elevation[free_surface] = fs_elevation
        return fs_elevation

