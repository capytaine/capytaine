################################################################################
#                                Green function                                #
################################################################################
import logging
from functools import lru_cache

import numpy as np

from capytaine.bem.prony_decomposition import find_best_exponential_decomposition
import capytaine.bem.NemohCore as NemohCore

tabulated_integrals = lru_cache(maxsize=1)(NemohCore.initialize_green_wave.initialize_tabulated_integrals)
LOG = logging.getLogger(__name__)

class Delhommeau:
    """
    Parameters
    ----------
    tabulation_nb_integration_points: int, optional
        Number of points for the evaluation of the tabulated elementary integrals w.r.t. :math:`theta`
        used for the computation of the Green function (default: 251)
    finite_depth_prony_decomposition_method: string, optional
        The implementation of the Prony decomposition used to compute the finite depth Green function.

    Attributes
    ----------
    tabulated_integrals: 3-ple of arrays
        Tabulated integrals for the computation of the Green function.
    """
    def __init__(self,
                 tabulation_nb_integration_points=251,
                 finite_depth_prony_decomposition_method='fortran',
                 ):
        self.tabulated_integrals = tabulated_integrals(328, 46, tabulation_nb_integration_points)

        self.finite_depth_prony_decomposition_method = finite_depth_prony_decomposition_method

        self.exportable_settings = {
            'green_function': 'Delhommeau',
            'tabulation_nb_integration_points': tabulation_nb_integration_points,
            'finite_depth_prony_decomposition_method': finite_depth_prony_decomposition_method,
        }

    def evaluate(self, mesh1, mesh2, free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0):
        Srankine, Vrankine = self.evaluate_rankine(mesh1, mesh2, free_surface, sea_bottom, wavenumber)

        if (free_surface == np.infty or
                (free_surface - sea_bottom == np.infty and wavenumber in (0, np.infty))):
            # No more terms in the Green function
            return Srankine, Vrankine

        Swave, Vwave = self.evaluate_wave(mesh1, mesh2, free_surface, sea_bottom, wavenumber)

        # The real valued matrices Srankine and Vrankine are automatically recasted as complex in the sum.
        Swave += Srankine
        Vwave += Vrankine
        return Swave, Vwave

    def evaluate_rankine(self, mesh1, mesh2, free_surface=0.0, sea_bottom=-np.infty, wavenumber=1.0):
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

    def evaluate_wave(self, mesh1, mesh2, free_surface, sea_bottom, wavenumber):
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
                method=self.finite_depth_prony_decomposition_method,
            )

            return NemohCore.green_wave.build_matrices_wave_source(
                mesh1.faces_centers, mesh1.faces_normals,
                mesh2.faces_centers, mesh2.faces_areas,
                wavenumber, depth,
                *self.tabulated_integrals,
                lamda_exp, a_exp,
                mesh1 is mesh2
                )

################################################################################
#                                   Engines                                    #
################################################################################
import logging

import numpy as np

from capytaine.matrices import linear_solvers
from capytaine.matrices.builders import identity_like

LOG = logging.getLogger(__name__)

class BasicEngine:
    """
    Parameters
    ----------
    matrix_cache_size: int, optional
        number of matrices to keep in cache
    linear_solver: str or function, optional
        Setting of the numerical solver for linear problems Ax = b.
        It can be set with the name of a preexisting solver
        (available: "direct" and "gmres", the latter is the default choice)
        or by passing directly a solver function.
    """
    available_linear_solvers = {'direct': linear_solvers.solve_directly,
                                'gmres': linear_solvers.solve_gmres}

    def __init__(self,
                 linear_solver=linear_solvers.solve_gmres,
                 matrix_cache_size=1,
                 ):

        if linear_solver in self.available_linear_solvers:
            self.linear_solver = self.available_linear_solvers[linear_solver]
        else:
            self.linear_solver = linear_solver

        if matrix_cache_size > 0:
            self.build_matrices = lru_cache(maxsize=matrix_cache_size)(self.build_matrices)

        self.exportable_settings = {
            'engine': 'BasicEngine',
            'matrix_cache_size': matrix_cache_size,
            'linear_solver': str(linear_solver),
        }

    def build_matrices(self, problem, green_function):
        S, V = green_function.evaluate(
            problem.body.mesh, problem.body.mesh,
            free_surface=problem.free_surface, sea_bottom=problem.sea_bottom, wavenumber=problem.wavenumber
        )

        return S, V + identity_like(V)/2

    def build_S_matrix_for_reconstruction(self, problem, mesh, green_function):
        if chunk_size > mesh.nb_faces:
            S, _ = green_function.evaluate(
                mesh,
                result.body.mesh,
                free_surface=result.free_surface,
                sea_bottom=result.sea_bottom,
                wavenumber=result.wavenumber
            )
            return S

        else:
            raise NotImplementedError

            for i in range(0, mesh.nb_faces, chunk_size):
                S, _ = green_function.evaluate(
                        mesh.extract_faces(list(range(i, i+chunk_size))),
                        result.body.mesh,
                        free_surface=result.free_surface,
                        sea_bottom=result.sea_bottom,
                        wavenumber=result.wavenumber
                    )
            return S



from capytaine.bem.hierarchical_toeplitz_matrices import hierarchical_toeplitz_matrices

class HierarchicalToeplitzMatrices:
    """

    Parameters
    ----------
    ACA_distance: float, optional
        Above this distance, the ACA is used to approximate the matrix with a low-rank block.
    ACA_tol: float, optional
        The tolerance of the ACA when building a low-rank matrix.
    matrix_cache_size: int, optional
        number of matrices to keep in cache
    """
    def __init__(self,
                 ACA_distance=np.infty,
                 ACA_tol=1e-2,
                 matrix_cache_size=1,
                 ):
        self.build_matrices = lru_cache(maxsize=matrix_cache_size)(self.build_matrices)
        self.ACA_distance = ACA_distance
        self.ACA_tol = ACA_tol

        self.linear_solver = solve_gmres
        self.build_matrices = hierarchical_toeplitz_matrices(
            BasicEngine.build_matrices,
            ACA_tol=ACA_tol,
            ACA_distance=ACA_distance,
            dtype=np.complex128
        )

        self.exportable_settings = {
            'engine': 'HierarchicalToeplitzMatrices',
            'ACA_distance': ACA_distance,
            'ACA_tol': ACA_tol,
            'matrix_cache_size': matrix_cache_size,
        }


################################################################################
#                                    Solver                                    #
################################################################################
import logging

import numpy as np

from datetime import datetime
from capytaine.io.xarray import problems_from_dataset, assemble_dataset, kochin_data_array

LOG = logging.getLogger(__name__)

class BEMSolver:

    def __init__(self,
                 green_function=Delhommeau(),
                 engine=BasicEngine(),
                 ):
        self.green_function = green_function
        self.engine = engine

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

        S, K = self.engine.build_matrices(problem, self.green_function)
        sources = self.engine.linear_solver(K, problem.boundary_condition)
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

    def exportable_settings(self):
        return {**self.green_function.exportable_settings,
                **self.engine.exportable_settings}

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
            They probably have not been stored by the solver because the option keep_details=True have not been set.
            Please re-run the resolution with this option.""")

        # Legacy
        self.engine.S_chunk_size = chunk_size

        S = self.engine.build_S_matrix_for_reconstruction(result.problem, mesh, self.green_function)
        phi = S @ result.sources

        LOG.debug(f"Done computing potential on {mesh.name} for {result}.")

        return phi

    def get_free_surface_elevation(self, result, free_surface, keep_details=False):
        """Compute the elevation of the free surface on a mesh for a previously solved problem.

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
        fs_elevation = 1j*result.omega/result.g * self.get_potential_on_mesh(result, free_surface.mesh)
        if keep_details:
            result.fs_elevation[free_surface] = fs_elevation
        return fs_elevation


# Legacy interface

def _arguments(f):
    return f.__code__.co_varnames[:f.__code__.co_argcount]

class Nemoh(BEMSolver):
    """Solver for the BEM problem based on Nemoh's Green function. Legacy API.
    Parameters are dispatched to the Delhommeau class and to the engine
    (BasicEngine or HierarchicalToeplitzMatrices).
    """

    def __init__(self, **params):
        self.green_function = Delhommeau(
           **{key: params[key] for key in params if key in _arguments(Delhommeau.__init__)}
	)
        if 'hierarchical_matrices' in params and params['hierarchical_matrices']:
            self.engine = HierarchicalToeplitzMatrices(
               **{key: params[key] for key in params if key in _arguments(HierarchicalToeplitzMatrices.__init__)}
            )
        else:
            self.engine = BasicEngine(
               **{key: params[key] for key in params if key in _arguments(BasicEngine.__init__)}
            )

    def build_matrices(self, *args, **kwargs):
        """Legacy API."""
        return self.engine.build_matrices(*args, **kwargs)

