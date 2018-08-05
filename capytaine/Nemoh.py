# coding: utf-8
# This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""Solver for the BEM problem based on Nemoh's Green function.

Example
-------

::

    problem = RadiationProblem(...)
    result = Nemoh().solve(problem)

"""

import logging
from functools import lru_cache
from copy import deepcopy

import numpy as np

from capytaine.tools.exponential_decomposition import find_best_exponential_decomposition
from capytaine.Toeplitz_matrices import identity_matrix_of_same_shape_as, solve
from capytaine.symmetries import use_symmetries
import capytaine.NemohCore as NemohCore


LOG = logging.getLogger(__name__)

FLOAT_PRECISION = np.float64


class Nemoh:
    """Solver for the BEM problem based on Nemoh's Green function.

    Parameters
    ----------
    store_matrices_in_cache: bool, optional
        If True, store the last computed influence matrices in __cache__ for later reuse (default: True)
    npinte: int, optional
        Number of points for the evaluation of the integral w.r.t. :math:`theta` in the Green function (default: 251)
    max_stored_exponential_decompositions: int, optional
        Number of stored exponential decomposition (default: 50)

    Attributes
    ----------
    XR: array of shape (328)
    XZ: array of shape (46)
    APD: array of shape (328, 46, 2, 2)
        Tabulated integrals for the Green functions
    __cache__: dict of dict of arrays
        Store last computations of influence matrices
    """
    def __init__(self, npinte=251):
        LOG.info("Initialize Nemoh's Green function.")
        self.XR, self.XZ, self.APD = NemohCore.initialize_green_wave.initialize_tabulated_integrals(328, 46, npinte)

    def solve(self, problem, keep_details=False):
        """Solve the BEM problem using Nemoh.

        Parameters
        ----------
        problem: LinearPotentialFlowProblem
            the problem to be solved
        keep_details: bool, optional
            if True, store the sources and the potential on the floating body in the output object (default: False)

        Returns
        -------
        LinearPotentialFlowResult
            an object storing the problem data and its results
        """

        LOG.info("Solve %s.", problem)

        if problem.wavelength < 8*problem.body.mesh.faces_radiuses.max():
            LOG.warning(f"Resolution of the mesh (8×max_radius={8*problem.body.mesh.faces_radiuses.max():.2e}) "
                        f"might be insufficient for this wavelength (wavelength={problem.wavelength:.2e})!")

        S, V = self.build_matrices(
            problem.body.mesh, problem.body.mesh,
            free_surface=problem.free_surface, sea_bottom=problem.sea_bottom, wavenumber=problem.wavenumber
        )

        identity = identity_matrix_of_same_shape_as(V)
        sources = solve(V + identity/2, problem.boundary_condition)
        potential = S @ sources

        result = problem.make_results_container()
        if keep_details:
            result.sources = sources
            result.potential = potential

        for influenced_dof_name, influenced_dof in problem.body.dofs.items():
            influenced_dof = np.sum(influenced_dof * problem.body.mesh.faces_normals, axis=1)
            integrated_potential = - problem.rho * potential @ (influenced_dof * problem.body.mesh.faces_areas)
            result.store_force(influenced_dof_name, integrated_potential)
            # Depending of the type of problem, the force will be kept as a complex-valued Froude-Krylov force
            # or stored as a couple of added mass and damping radiation coefficients.

        LOG.debug("Done!")

        return result

    def solve_all(self, problems, processes=1):
        """Solve several problems in parallel.

        Running::

            solver.solve_all(problems)

        is more or less equivalent to::

             [solver.solve(problem) for problem in problems]

        but in parallel with some optimizations for faster resolution.

        Parameters
        ----------
        problems: list of LinearPotentialFlowProblem
            several problems to be solved
        processes: int, optional
            number of parallel processes (default: 1)

        Return
        ------
        list of LinearPotentialFlowResult
            the solved problems
        """
        from multiprocessing import Pool
        with Pool(processes=processes) as pool:
            results = pool.map(self.solve, sorted(problems))
        return results

    #######################
    #  Building matrices  #
    #######################

    def build_matrices(self, mesh1, mesh2, free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0):
        """
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
        matrix-like
            couple of influence matrices
        """

        LOG.debug(f"\tEvaluating matrix of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name}"
                  f"for depth={free_surface-sea_bottom} and wavenumber={wavenumber}.")

        # S, V = self._build_matrices_0(mesh1, mesh2)
        S0, V0 = self._build_matrices_0(mesh1, mesh2)
        S = deepcopy(S0)
        V = deepcopy(V0)
        # Do a copy to avoid interfering with the cache. 
        # TODO: fix that when both matrices are in a single array.

        if free_surface < np.infty:

            S1, V1 = self._build_matrices_1(mesh1, mesh2, free_surface, sea_bottom)
            S += S1
            V += V1

            S2, V2 = self._build_matrices_2(mesh1, mesh2, free_surface, sea_bottom, wavenumber)
            S += S2
            V += V2

        return S, V

    @lru_cache(maxsize=1)
    @use_symmetries
    def _build_matrices_0(self, mesh1, mesh2):
        """Compute the first part of the influence matrices of mesh1 on mesh2."""
        S, V = NemohCore.green_rankine.build_matrices_rankine_source(
            mesh1.faces_centers, mesh1.faces_normals,
            mesh2.vertices,      mesh2.faces + 1,
            mesh2.faces_centers, mesh2.faces_normals,
            mesh2.faces_areas,   mesh2.faces_radiuses,
            )
        return S.astype(np.complex128), V.astype(np.complex128)
        # The real valued matrices are converted to complex to ease the combination with the other matrices later.
        # TODO: Move the conversion outside of the scope of the cache.

    @lru_cache(maxsize=1)
    @use_symmetries
    def _build_matrices_1(self, mesh1, mesh2, free_surface, sea_bottom):
        """Compute the second part of the influence matrices of mesh1 on mesh2."""
        depth = free_surface - sea_bottom

        def reflect_vector(x):
            y = x.copy()
            y[:, 2] = -x[:, 2]
            return y

        if depth == np.infty:
            def reflect_point(x):
                y = x.copy()
                y[:, 2] = 2*free_surface - x[:, 2]
                return y
        else:
            def reflect_point(x):
                y = x.copy()
                y[:, 2] = 2*sea_bottom - x[:, 2]
                return y

        S1, V1 = NemohCore.green_rankine.build_matrices_rankine_source(
            reflect_point(mesh1.faces_centers), reflect_vector(mesh1.faces_normals),
            mesh2.vertices,      mesh2.faces + 1,
            mesh2.faces_centers, mesh2.faces_normals,
            mesh2.faces_areas,   mesh2.faces_radiuses,
            )

        if depth == np.infty:
            return -S1.astype(np.complex128), -V1.astype(np.complex128)
        else:
            return S1.astype(np.complex128), V1.astype(np.complex128)

    @lru_cache(maxsize=1)
    @use_symmetries
    def _build_matrices_2(self, mesh1, mesh2, free_surface, sea_bottom, wavenumber):
        """Compute the third part (wave part) of the influence matrices of mesh1 on mesh2."""
        depth = free_surface - sea_bottom
        if depth == np.infty:
            S2, V2 = NemohCore.green_wave.build_matrices_wave_source(
                mesh1.faces_centers, mesh1.faces_normals,
                mesh2.faces_centers, mesh2.faces_areas,
                wavenumber, 0.0,
                self.XR, self.XZ, self.APD,
                np.empty(1), np.empty(1),  # Dummy arrays that won't actually be used by the fortran code.
                mesh1 is mesh2
                )
        else:
            a_exp, lamda_exp = find_best_exponential_decomposition(wavenumber*depth*np.tanh(wavenumber*depth),
                                                                   wavenumber*depth)

            S2, V2 = NemohCore.green_wave.build_matrices_wave_source(
                mesh1.faces_centers, mesh1.faces_normals,
                mesh2.faces_centers, mesh2.faces_areas,
                wavenumber, depth,
                self.XR, self.XZ, self.APD,
                lamda_exp, a_exp,
                mesh1 is mesh2
                )

        return S2, V2

    #######################
    #  Compute potential  #
    #######################

    def get_potential_on_mesh(self, result, mesh):
        """Compute the potential on a mesh for the potential field of a previously solved problem.

        Parameters
        ----------
        result : LinearPotentialFlowResult
            the return of Nemoh's solver
        mesh : Mesh or CollectionOfMeshes
            a mesh

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

        S, _ = self.build_matrices(
            mesh,
            result.body.mesh,
            free_surface=result.free_surface,
            sea_bottom=result.sea_bottom,
            wavenumber=result.wavenumber
        )

        phi = S @ result.sources

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

