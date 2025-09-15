"""Definition of the methods to build influence matrices, using possibly some sparse structures."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging
from abc import ABC, abstractmethod
from typing import Tuple, Union, Optional, Callable

import numpy as np
from scipy.linalg import lu_factor
from scipy.sparse import coo_matrix
from scipy.sparse import linalg as ssl

from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.meshes.symmetric import ReflectionSymmetricMesh, TranslationalSymmetricMesh, AxialSymmetricMesh

from capytaine.matrices import linear_solvers
from capytaine.matrices.block import BlockMatrix
from capytaine.matrices.low_rank import LowRankMatrix, NoConvergenceOfACA
from capytaine.matrices.block_toeplitz import BlockSymmetricToeplitzMatrix, BlockToeplitzMatrix, BlockCirculantMatrix

from capytaine.green_functions.abstract_green_function import AbstractGreenFunction
from capytaine.green_functions.delhommeau import Delhommeau

from capytaine.tools.block_circulant_matrices import (
        BlockCirculantMatrix as NewBlockCirculantMatrix,
        lu_decompose, has_been_lu_decomposed,
        MatrixLike, LUDecomposedMatrixLike
        )
from capytaine.tools.lru_cache import lru_cache_with_strict_maxsize

LOG = logging.getLogger(__name__)


####################
#  ABSTRACT CLASS  #
####################

class MatrixEngine(ABC):
    """Abstract method to build a matrix."""

    @abstractmethod
    def build_matrices(self, mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer):
        pass

    @abstractmethod
    def build_S_matrix(self, mesh1, mesh2, free_surface, water_depth, wavenumber):
        pass

    @abstractmethod
    def build_fullK_matrix(self, mesh1, mesh2, free_surface, water_depth, wavenumber):
        pass


##################
#  BASIC ENGINE  #
##################

class Counter:
    def __init__(self):
        self.nb_iter = 0

    def __call__(self, *args, **kwargs):
        self.nb_iter += 1


def solve_gmres(A, b):
    LOG.debug(f"Solve with GMRES for {A}.")

    if LOG.isEnabledFor(logging.INFO):
        counter = Counter()
        x, info = ssl.gmres(A, b, atol=1e-6, callback=counter)
        LOG.info(f"End of GMRES after {counter.nb_iter} iterations.")

    else:
        x, info = ssl.gmres(A, b, atol=1e-6)

    if info > 0:
        raise RuntimeError(f"No convergence of the GMRES after {info} iterations.\n"
                            "This can be due to overlapping panels or irregular frequencies.")

    return x


LUDecomposedMatrixOrNot = Union[MatrixLike, LUDecomposedMatrixLike]


class BasicMatrixEngine(MatrixEngine):
    """
    Default matrix engine.

    Features:
        - Caching of the last computed matrices.
        - Supports plane symmetries and nested plane symmetries.
        - Linear solver can be customized. Default is `lu_decomposition` with caching of the LU decomposition.

    Parameters
    ----------
    green_function: AbstractGreenFunction
        the low level implementation used to compute the coefficients of the matrices.
    linear_solver: str or function, optional
        Setting of the numerical solver for linear problems Ax = b.
        It can be set with the name of a preexisting solver
        (available: "lu_decomposition" and "gmres", the former is the default choice)
        or by passing directly a solver function.
    """

    green_function: AbstractGreenFunction
    _linear_solver: Union[str, Callable]
    last_computed_matrices: Optional[Tuple[MatrixLike, LUDecomposedMatrixOrNot]]

    def __init__(self, *, green_function=None, linear_solver='lu_decomposition'):

        self.green_function = Delhommeau() if green_function is None else green_function

        self._linear_solver = linear_solver

        self.last_computed_inputs = None
        self.last_computed_matrices = None

        self.exportable_settings = {
            'engine': 'BasicMatrixEngine',
            'linear_solver': str(linear_solver),
            **self.green_function.exportable_settings,
        }

    def __str__(self):
        params= [f"green_function={self.green_function}", f"linear_solver={repr(self._linear_solver)}"]
        return f"BasicMatrixEngine({', '.join(params)})"

    def __repr__(self):
        return self.__str__()

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    def build_S_matrix(
            self, mesh1, mesh2, free_surface, water_depth, wavenumber
            ) -> np.ndarray:
        """Similar to :code:`build_matrices`, but returning only :math:`S`"""
        # Calls directly evaluate instead of build_matrices because the caching
        # mechanism of build_matrices is not compatible with giving mesh1 as a
        # list of points, but we need that for post-processing
        S, _ = self.green_function.evaluate(
                mesh1, mesh2, free_surface, water_depth, wavenumber
            )
        return S

    def build_fullK_matrix(
            self, mesh1, mesh2, free_surface, water_depth, wavenumber
            ) -> np.ndarray:
        """Similar to :code:`build_matrices`, but returning only full :math:`K`
        (that is the three components of the gradient, not just the normal one)"""
        # TODO: could use symmetries. In particular for forward, we compute the
        # full velocity on the same mesh so symmetries could be used.
        _, fullK = self.green_function.evaluate(
                mesh1, mesh2, free_surface, water_depth, wavenumber,
                adjoint_double_layer=True, early_dot_product=False,
            )
        return fullK

    def _build_matrices_with_symmetries(
            self, mesh1, mesh2, *, free_surface, water_depth, wavenumber, adjoint_double_layer
            ) -> Tuple[MatrixLike, MatrixLike]:
        if (isinstance(mesh1, ReflectionSymmetricMesh)
                and isinstance(mesh2, ReflectionSymmetricMesh)
                and mesh1.plane == mesh2.plane):

            S_a, K_a = self._build_matrices_with_symmetries(
                mesh1[0], mesh2[0],
                free_surface=free_surface, water_depth=water_depth,
                wavenumber=wavenumber, adjoint_double_layer=adjoint_double_layer
                )
            S_b, K_b = self._build_matrices_with_symmetries(
                mesh1[0], mesh2[1],
                free_surface=free_surface, water_depth=water_depth,
                wavenumber=wavenumber, adjoint_double_layer=adjoint_double_layer
                )

            return NewBlockCirculantMatrix([S_a, S_b]), NewBlockCirculantMatrix([K_a, K_b])

        else:
            return self.green_function.evaluate(
                mesh1, mesh2, free_surface, water_depth, wavenumber,
                adjoint_double_layer=adjoint_double_layer, early_dot_product=True,
            )

    def _build_and_cache_matrices_with_symmetries(
            self, mesh1, mesh2, *, free_surface, water_depth, wavenumber, adjoint_double_layer
            ) -> Tuple[MatrixLike, LUDecomposedMatrixOrNot]:
        if (mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer) == self.last_computed_inputs:
            LOG.debug("%s: reading cache.", self.__class__.__name__)
            return self.last_computed_matrices
        else:
            LOG.debug("%s: computing new matrices.", self.__class__.__name__)
            self.last_computed_matrices = None  # Unlink former cached values, so the memory can be freed to compute new matrices.
            S, K = self._build_matrices_with_symmetries(
                    mesh1, mesh2,
                    free_surface=free_surface, water_depth=water_depth,
                    wavenumber=wavenumber, adjoint_double_layer=adjoint_double_layer
                    )
            self.last_computed_inputs = (mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer)
            self.last_computed_matrices = (S, K)
            return self.last_computed_matrices

    # Main interface for compliance with AbstractGreenFunction interface
    def build_matrices(self, mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer=True):
        r"""Build the influence matrices between mesh1 and mesh2.

        Parameters
        ----------
        mesh1: MeshLike or list of points
            mesh of the receiving body (where the potential is measured)
        mesh2: MeshLike
            mesh of the source body (over which the source distribution is integrated)
        free_surface: float
            position of the free surface (default: :math:`z = 0`)
        water_depth: float
            position of the sea bottom (default: :math:`z = -\infty`)
        wavenumber: float
            wavenumber (default: 1.0)
        adjoint_double_layer: bool, optional
            compute double layer for direct method (F) or adjoint double layer for indirect method (T) matrices (default: True)

        Returns
        -------
        tuple of matrix-like (Numpy arrays or BlockCirculantMatrix)
            the matrices :math:`S` and :math:`K`
        """
        return self._build_and_cache_matrices_with_symmetries(
            mesh1, mesh2,
            free_surface=free_surface, water_depth=water_depth,
            wavenumber=wavenumber, adjoint_double_layer=adjoint_double_layer
        )

    def linear_solver(self, A: LUDecomposedMatrixOrNot, b: np.ndarray) -> np.ndarray:
        """Solve a linear system with left-hand side A and right-hand-side b

        Parameters
        ----------
        A: matrix-like
            Expected to be the second output of `build_matrices`
        b: np.ndarray
            Vector of the correct length

        Returns
        -------
        x: np.ndarray
            Vector such that A@x = b
        """
        if not isinstance(self._linear_solver, str):
            # If not a string, it is expected to be a custom function that can
            # be called to solve the system
            x = self._linear_solver(A, b)

            if not x.shape == b.shape:
                raise ValueError(f"Error in linear solver of {self}: the shape of the output ({x.shape}) "
                                 f"does not match the expected shape ({b.shape})")

            return x

        elif self._linear_solver == "lu_decomposition":
            if not has_been_lu_decomposed(A):
                luA = lu_decompose(A)
                if A is self.last_computed_matrices[1]:
                    # In normal operation of Capytaine, `A` is always the $D$
                    # or $K$ matrix stored in the cache of the solver.
                    # Here we replace the matrix by its LU decomposition in the
                    # cache to avoid doing the decomposition again.
                    self.last_computed_matrices = (self.last_computed_matrices[0], luA)
            else:
                luA: LUDecomposedMatrixLike = A
            return luA.solve(b)

        elif self._linear_solver == "gmres":
            return solve_gmres(A, b)

        else:
            raise NotImplementedError(
                f"Unknown `linear_solver` in BasicMatrixEngine: {self._linear_solver}"
            )


###################################
#  HIERARCHIAL TOEPLITZ MATRICES  #
###################################

class HierarchicalToeplitzMatrixEngine(MatrixEngine):
    """An experimental matrix engine that build a hierarchical matrix with
     some block-Toeplitz structure.

    Parameters
    ----------
    ACA_distance: float, optional
        Above this distance, the ACA is used to approximate the matrix with a low-rank block.
    ACA_tol: float, optional
        The tolerance of the ACA when building a low-rank matrix.
    matrix_cache_size: int, optional
        number of matrices to keep in cache
    """

    def __init__(self, *, green_function=None, ACA_distance=8.0, ACA_tol=1e-2, matrix_cache_size=1):

        self.green_function = Delhommeau() if green_function is None else green_function

        if matrix_cache_size > 0:
            self.build_matrices = lru_cache_with_strict_maxsize(maxsize=matrix_cache_size)(self.build_matrices)

        self.ACA_distance = ACA_distance
        self.ACA_tol = ACA_tol

        self.linear_solver = linear_solvers.solve_gmres

        self.exportable_settings = {
            'engine': 'HierarchicalToeplitzMatrixEngine',
            'ACA_distance': ACA_distance,
            'ACA_tol': ACA_tol,
            'matrix_cache_size': matrix_cache_size,
            **self.green_function.exportable_settings
        }

    def __str__(self):
        params = [f"green_function={self.green_function}", f"ACA_distance={self.ACA_distance}"]
        params.append(f"ACA_tol={self.ACA_tol}")
        if self.exportable_settings['matrix_cache_size'] != 1:
            params.append(f"matrix_cache_size={self.exportable_settings['matrix_cache_size']}")
        return f"HierarchicalToeplitzMatrixEngine({','.join(params)})"

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    def build_S_matrix(self, mesh1, mesh2, free_surface, water_depth, wavenumber):
        """Similar to :code:`build_matrices`, but returning only :math:`S`"""
        S, _ = self.green_function.evaluate(
                mesh1, mesh2, free_surface, water_depth, wavenumber
            )
        return S

    def build_fullK_matrix(self, mesh1, mesh2, free_surface, water_depth, wavenumber):
        """Similar to :code:`build_matrices`, but returning only full :math:`K`
        (that is the three components of the gradient, not just the normal one)"""
        _, fullK = self.green_function.evaluate(
                mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer=True, early_dot_product=False,
            )
        return fullK

    def build_matrices(self,
                       mesh1, mesh2, free_surface, water_depth, wavenumber,
                       adjoint_double_layer=True):

        return self._build_matrices(
                           mesh1, mesh2, free_surface, water_depth, wavenumber,
                           adjoint_double_layer, _rec_depth=1)


    def _build_matrices(self,
                       mesh1, mesh2, free_surface, water_depth, wavenumber,
                       adjoint_double_layer, _rec_depth=1):
        """Recursively builds a hierarchical matrix between mesh1 and mesh2.

        Same arguments as :func:`BasicMatrixEngine.build_matrices`.

        :code:`_rec_depth` keeps track of the recursion depth only for pretty log printing.
        """

        if logging.getLogger().isEnabledFor(logging.DEBUG):
            log_entry = (
                "\t" * (_rec_depth+1) +
                "Build the S and K influence matrices between {mesh1} and {mesh2}"
                .format(mesh1=mesh1.name, mesh2=(mesh2.name if mesh2 is not mesh1 else 'itself'))
            )
        else:
            log_entry = ""  # will not be used

        # Distance between the meshes (for ACA).
        distance = np.linalg.norm(mesh1.center_of_mass_of_nodes - mesh2.center_of_mass_of_nodes)

        # I) SPARSE COMPUTATION
        # I-i) BLOCK TOEPLITZ MATRIX

        if (isinstance(mesh1, ReflectionSymmetricMesh)
                and isinstance(mesh2, ReflectionSymmetricMesh)
                and mesh1.plane == mesh2.plane):

            LOG.debug(log_entry + " using mirror symmetry.")

            S_a, V_a = self._build_matrices(
                mesh1[0], mesh2[0], free_surface, water_depth, wavenumber,
                adjoint_double_layer=adjoint_double_layer, _rec_depth=_rec_depth+1)
            S_b, V_b = self._build_matrices(
                mesh1[0], mesh2[1], free_surface, water_depth, wavenumber,
                adjoint_double_layer=adjoint_double_layer, _rec_depth=_rec_depth+1)

            return BlockSymmetricToeplitzMatrix([[S_a, S_b]]), BlockSymmetricToeplitzMatrix([[V_a, V_b]])

        elif (isinstance(mesh1, TranslationalSymmetricMesh)
              and isinstance(mesh2, TranslationalSymmetricMesh)
              and np.allclose(mesh1.translation, mesh2.translation)
              and mesh1.nb_submeshes == mesh2.nb_submeshes):

            LOG.debug(log_entry + " using translational symmetry.")

            S_list, V_list = [], []
            for submesh in mesh2:
                S, V = self._build_matrices(
                    mesh1[0], submesh, free_surface, water_depth, wavenumber,
                    adjoint_double_layer=adjoint_double_layer, _rec_depth=_rec_depth+1)
                S_list.append(S)
                V_list.append(V)
            for submesh in mesh1[1:][::-1]:
                S, V = self._build_matrices(
                    submesh, mesh2[0], free_surface, water_depth, wavenumber,
                    adjoint_double_layer=adjoint_double_layer, _rec_depth=_rec_depth+1)
                S_list.append(S)
                V_list.append(V)

            return BlockToeplitzMatrix([S_list]), BlockToeplitzMatrix([V_list])

        elif (isinstance(mesh1, AxialSymmetricMesh)
              and isinstance(mesh2, AxialSymmetricMesh)
              and mesh1.axis == mesh2.axis
              and mesh1.nb_submeshes == mesh2.nb_submeshes):

            LOG.debug(log_entry + " using rotation symmetry.")

            S_line, V_line = [], []
            for submesh in mesh2[:mesh2.nb_submeshes]:
                S, V = self._build_matrices(
                    mesh1[0], submesh, free_surface, water_depth, wavenumber,
                    adjoint_double_layer=adjoint_double_layer, _rec_depth=_rec_depth+1)
                S_line.append(S)
                V_line.append(V)

            return BlockCirculantMatrix([S_line]), BlockCirculantMatrix([V_line])

        # I-ii) LOW-RANK MATRIX WITH ACA

        elif distance > self.ACA_distance*mesh1.diameter_of_nodes or distance > self.ACA_distance*mesh2.diameter_of_nodes:

            LOG.debug(log_entry + " using ACA.")

            def get_row_func(i):
                s, v = self.green_function.evaluate(
                    mesh1.extract_one_face(i), mesh2,
                    free_surface, water_depth, wavenumber,
                    adjoint_double_layer=adjoint_double_layer
                )
                return s.flatten(), v.flatten()

            def get_col_func(j):
                s, v = self.green_function.evaluate(
                    mesh1, mesh2.extract_one_face(j),
                    free_surface, water_depth, wavenumber,
                    adjoint_double_layer=adjoint_double_layer
                )
                return s.flatten(), v.flatten()

            try:
                return LowRankMatrix.from_rows_and_cols_functions_with_multi_ACA(
                    get_row_func, get_col_func, mesh1.nb_faces, mesh2.nb_faces,
                    nb_matrices=2, id_main=1,  # Approximate V and get an approximation of S at the same time
                    tol=self.ACA_tol, dtype=np.complex128)
            except NoConvergenceOfACA:
                pass  # Continue with non sparse computation

        # II) NON-SPARSE COMPUTATIONS
        # II-i) BLOCK MATRIX

        if (isinstance(mesh1, CollectionOfMeshes)
              and isinstance(mesh2, CollectionOfMeshes)):

            LOG.debug(log_entry + " using block matrix structure.")

            S_matrix, V_matrix = [], []
            for submesh1 in mesh1:
                S_line, V_line = [], []
                for submesh2 in mesh2:
                    S, V = self._build_matrices(
                        submesh1, submesh2, free_surface, water_depth, wavenumber,
                        adjoint_double_layer=adjoint_double_layer, _rec_depth=_rec_depth+1)

                    S_line.append(S)
                    V_line.append(V)
                S_matrix.append(S_line)
                V_matrix.append(V_line)

            return BlockMatrix(S_matrix), BlockMatrix(V_matrix)

        # II-ii) PLAIN NUMPY ARRAY

        else:
            LOG.debug(log_entry)

            S, V = self.green_function.evaluate(
                mesh1, mesh2, free_surface, water_depth, wavenumber, adjoint_double_layer=adjoint_double_layer
            )
            return S, V

class HierarchicalPrecondMatrixEngine(HierarchicalToeplitzMatrixEngine):
    """An experimental matrix engine that build a hierarchical matrix with
     some block-Toeplitz structure.

    Parameters
    ----------
    ACA_distance: float, optional
        Above this distance, the ACA is used to approximate the matrix with a low-rank block.
    ACA_tol: float, optional
        The tolerance of the ACA when building a low-rank matrix.
    matrix_cache_size: int, optional
        number of matrices to keep in cache
    """

    def __init__(self, *, ACA_distance=8.0, ACA_tol=1e-2, matrix_cache_size=1):
        super().__init__(ACA_distance=ACA_distance, ACA_tol=ACA_tol, matrix_cache_size=matrix_cache_size)
        self.linear_solver = linear_solvers.solve_precond_gmres

    def build_matrices(self,
                       mesh1, mesh2, free_surface, water_depth, wavenumber,
                       green_function, adjoint_double_layer=True):
        """Recursively builds a hierarchical matrix between mesh1 and mesh2,
        and precomputes some of the quantities needed for the preconditioner.

        Same arguments as :func:`BasicMatrixEngine.build_matrices`, except for rec_depth
        """
        # Build the matrices using the method of the parent class
        S, K = super().build_matrices(mesh1, mesh2, free_surface, water_depth,
                                      wavenumber, green_function,
                                      adjoint_double_layer=adjoint_double_layer)

        path_to_leaf = mesh1.path_to_leaf()

        n = len(path_to_leaf)
        N = K.shape[0]

        # Navigate to the diagonal blocks and compute their LU decompositions
        DLU = []
        diag_shapes = []
        for leaf in range(n):
            # Navigate to the block containing the one we need
            # (one layer above in the dendrogram)
            #upper_block = self.access_block_by_path(K, path_to_leaf[leaf][:-1])
            upper_block = K.access_block_by_path(path_to_leaf[leaf][:-1])
            # find the local index in the full path
            ind = path_to_leaf[leaf][-1]
            # compute the LU decomposition and add to the list
            DLU.append(lu_factor(upper_block.all_blocks[ind, ind]))
            diag_shapes.append(upper_block.all_blocks[ind, ind].shape[0])

        # Build the restriction and precompute its multiplication by K
        R = np.zeros((n, N), dtype=complex)
        RA = np.zeros((n, N), dtype=complex)
        for ii in range(n):
            row_slice = slice(sum(diag_shapes[:ii]), sum(diag_shapes[:ii+1]))
            R[ii, row_slice] = 1
            # Compute the multiplication using only the relevant slices of K
            # The slices are found by navigating the tree
            #RA[ii, :] = self.slice_rmatvec(R[ii, :], ii)
            Aloc = K
            v = R[ii, :]
            va = np.zeros(N, dtype=complex)
            free = [0, N]

            for lvl, jj in enumerate(path_to_leaf[ii]):

                Nrows = Aloc.all_blocks[jj, jj].shape[0]

                if jj==0:
                    v = v[:Nrows]
                    w = v @ Aloc.all_blocks[0,1]
                    va[free[1]-len(w) : free[1]] = w
                    free[1] = free[1] - len(w)
                else:
                    v = v[-Nrows:]
                    w = v @ Aloc.all_blocks[1, 0]
                    va[free[0] : free[0]+len(w)] = w
                    free[0] = free[0] + len(w)

                Aloc = Aloc.all_blocks[jj, jj]

                if lvl == len(path_to_leaf[ii])-1:
                    w = v@Aloc
                    va[free[0] : free[1]] = w
                    free[0] = free[0] + len(w)

            RA[ii, :] = va

        Ac = RA @ R.T
        AcLU = lu_factor(Ac)

        # Now navigate again to the diagonal blocks and set them to zero
        for leaf in range(n):
            upper_block = K.access_block_by_path(path_to_leaf[leaf][:-1])
            ind = path_to_leaf[leaf][-1]
            # turn the diagonal block into a zero sparse matrix
            upper_block.all_blocks[ind, ind] = coo_matrix(upper_block.all_blocks[ind, ind].shape)

        def PinvA_mv(v):
            v = v + 1j*np.zeros(N)
            return v - linear_solvers._block_Jacobi_coarse_corr(
                             K, np.zeros(N, dtype=complex), v,
                             R, RA, AcLU, DLU, diag_shapes, n)

        PinvA = ssl.LinearOperator((N, N), matvec=PinvA_mv)

        return S, (K, R, RA, AcLU, DLU, diag_shapes, n, PinvA)
