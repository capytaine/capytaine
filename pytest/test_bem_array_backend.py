# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Check that the solver works with a Green function returning matrices from another array library.
`array_api_strict` is used as the second library: any left-over call to NumPy on the matrices fails loudly."""


import pytest

import numpy as np
import scipy.linalg as sl
import xarray as xr
import capytaine as cpt
import capytaine_test_helpers as helpers
from capytaine.green_functions.abstract_green_function import AbstractGreenFunction, GreenFunctionEvaluationError
from capytaine.tools.block_circulant_matrices import AbstractLUDecomposedMatrix, lu_decompose
from capytaine.tools.lazy_matrices import LazyMatrix

xps = pytest.importorskip("array_api_strict")  # Not available on Python 3.8


class StrictArrayGreenFunction(AbstractGreenFunction):
    """Delhommeau Green function returning `array_api_strict` arrays instead of NumPy arrays."""
    _default_parameters = {}
    matrices_array_backend = xps
    matrices_device = xps.asarray(0.0).device

    def __init__(self, *, nan=False, floating_point_precision="float64"):
        self.nan = nan
        self.floating_point_precision = floating_point_precision
        self.wrapped = helpers.green_function(floating_point_precision=floating_point_precision)
        self.exportable_settings = {
            "green_function": "StrictArrayGreenFunction", "nan": nan,
            "floating_point_precision": floating_point_precision,
        }

    def evaluate(self, *args, **kwargs):
        S, K = self.wrapped.evaluate(*args, **kwargs)
        if self.nan:
            S[0, 0] = np.nan
        return xps.asarray(S), xps.asarray(K)


class LUDecomposedStrictArray(AbstractLUDecomposedMatrix):
    """LU decomposition of an `array_api_strict` matrix, as a plugin would register for its own array library."""
    def __init__(self, A, *, overwrite_a=False):
        self._lu = sl.lu_factor(np.asarray(A))
        self.shape = A.shape
        self.dtype = A.dtype
        self._device = A.device

    def __array_namespace__(self, *, api_version=None):
        return xps

    @property
    def device(self):
        return self._device

    def solve(self, b):
        return xps.asarray(sl.lu_solve(self._lu, np.asarray(b)), device=self._device)


lu_decompose.register(type(xps.asarray(0.0)), LUDecomposedStrictArray)


def solvers(method):
    return (
        helpers.solver(method=method),
        cpt.BEMSolver(green_function=StrictArrayGreenFunction(), method=method),
    )


def unwrap(x):
    while hasattr(x, "symbol"):
        x = x.value
    return x


def assert_same_results(res, res_strict):
    for dof in res.forces:
        force, force_strict = res.forces[dof], res_strict.forces[dof]
        assert getattr(force_strict, "symbol", None) == getattr(force, "symbol", None)  # Zero and infinite frequency
        force, force_strict = unwrap(force), unwrap(force_strict)
        assert isinstance(force_strict, (complex, np.number)), "The forces should not be arrays of the other library"
        assert force_strict == pytest.approx(force, rel=1e-8)


@pytest.mark.parametrize("method", ["direct", "indirect"])
@pytest.mark.parametrize("omega", [1.0, 0.0, np.inf])
def test_radiation_problem(method, omega):
    solver, solver_strict = solvers(method)
    body = helpers.small_sphere_body()
    pb = cpt.RadiationProblem(body=body, omega=omega, radiating_dof="Surge")
    res = solver.solve(pb, keep_details=True)
    res_strict = solver_strict.solve(pb, keep_details=True)
    assert_same_results(res, res_strict)
    # Everything stored in the results is NumPy only
    for name in ["potential", "pressure"]:
        actual, expected = unwrap(getattr(res_strict, name)), unwrap(getattr(res, name))
        assert isinstance(actual, np.ndarray)
        np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-10)


@pytest.mark.parametrize("method", ["direct", "indirect"])
def test_diffraction_problem(method):
    solver, solver_strict = solvers(method)
    pb = cpt.DiffractionProblem(body=helpers.small_sphere_body(), omega=1.0, wave_direction=0.3)
    assert_same_results(solver.solve(pb), solver_strict.solve(pb))


def test_forward_speed():
    solver, solver_strict = solvers("indirect")
    pb = cpt.RadiationProblem(body=helpers.small_sphere_body(), wavelength=3.0, radiating_dof="Surge",
                              forward_speed=1.0, wave_direction=np.pi)
    assert_same_results(solver.solve(pb), solver_strict.solve(pb))


@pytest.mark.parametrize("nb_points", [3, 600])  # Above 500 points, the matrix is a LazyMatrix
def test_potential_and_velocity_at_points(nb_points):
    solver, solver_strict = solvers("indirect")
    pb = cpt.RadiationProblem(body=helpers.small_sphere_body(), omega=1.0, radiating_dof="Surge")
    res = solver.solve(pb, keep_details=True)
    res_strict = solver_strict.solve(pb, keep_details=True)
    rng = np.random.default_rng(0)
    points = np.stack([
        rng.uniform(-4, 4, nb_points), rng.uniform(-4, 4, nb_points), rng.uniform(-5, -3.5, nb_points)
    ], axis=-1)
    for compute in ["compute_potential", "compute_velocity"]:
        expected = getattr(solver, compute)(points, res)
        actual = getattr(solver_strict, compute)(points, res_strict)
        assert isinstance(actual, np.ndarray)
        np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=1e-10)


def test_lazy_matrix_of_strict_arrays_is_used_for_many_points():
    engine = cpt.DefaultMatrixEngine(green_function=StrictArrayGreenFunction())
    points = np.zeros((600, 3)) + [3.0, 0.0, -1.0]
    points[:, 0] += np.linspace(0, 1, 600)
    mesh = helpers.small_sphere_mesh()
    S = engine.build_S_matrix(points, mesh, free_surface=0.0, water_depth=np.inf, wavenumber=1.0)
    assert isinstance(S, LazyMatrix)
    assert S.__array_namespace__() is xps
    assert S.dtype == xps.complex128
    assert S.shape == (600, mesh.nb_faces)
    # Conversion to a NumPy array does not use the dtype of the other library
    S_numpy = cpt.DefaultMatrixEngine(green_function=helpers.green_function()).build_S_matrix(points, helpers.small_sphere_mesh(), free_surface=0.0, water_depth=np.inf, wavenumber=1.0)
    np.testing.assert_allclose(np.asarray(S), np.asarray(S_numpy), rtol=1e-8)


def test_fill_dataset():
    solver, solver_strict = solvers("indirect")
    kwargs = dict(wave_direction=[0.0], radiating_dof=["Surge"], omega=[0.5, 1.0])
    dataset = xr.Dataset(coords={k: v for k, v in kwargs.items()})
    ds = solver.fill_dataset(dataset, helpers.small_sphere_body())
    ds_strict = solver_strict.fill_dataset(dataset, helpers.small_sphere_body())
    np.testing.assert_allclose(ds_strict["added_mass"].values, ds["added_mass"].values, rtol=1e-8)
    np.testing.assert_allclose(ds_strict["radiation_damping"].values, ds["radiation_damping"].values, rtol=1e-8, atol=1e-12)
    np.testing.assert_allclose(ds_strict["excitation_force"].values, ds["excitation_force"].values, rtol=1e-8)


def test_nan_in_matrix_of_another_array_library():
    solver = cpt.BEMSolver(green_function=StrictArrayGreenFunction(nan=True))
    pb = cpt.RadiationProblem(body=helpers.small_sphere_body(), omega=1.0, radiating_dof="Surge")
    with pytest.raises(GreenFunctionEvaluationError):
        solver.solve(pb)


def test_gmres_with_matrix_of_another_array_library():
    solver = cpt.BEMSolver(engine=cpt.DefaultMatrixEngine(green_function=StrictArrayGreenFunction(), linear_solver="gmres"))
    pb = cpt.RadiationProblem(body=helpers.small_sphere_body(), omega=1.0, radiating_dof="Surge")
    with pytest.raises(NotImplementedError):
        solver.solve(pb)


@pytest.mark.parametrize("method", ["direct", "indirect"])
def test_several_dofs_reuse_the_cached_lu_decomposition(method):
    solver, solver_strict = solvers(method)
    problems = [cpt.RadiationProblem(body=helpers.small_sphere_body("rigid"), omega=1.0, radiating_dof=dof)
                for dof in ["Surge", "Surge", "Pitch"]]
    for res, res_strict in zip(solver.solve_all(problems), solver_strict.solve_all(problems)):
        assert_same_results(res, res_strict)
    assert isinstance(solver_strict.engine.last_computed_matrices[1], LUDecomposedStrictArray)


def test_potential_without_free_surface_keeps_the_imaginary_part_of_the_sources():
    # Without free surface, the matrix S is real while the sources are complex.
    solver, solver_strict = (cpt.BEMSolver(green_function=gf) for gf in [helpers.green_function(), StrictArrayGreenFunction()])
    pb = cpt.RadiationProblem(body=helpers.small_sphere_body(), free_surface=np.inf, omega=1.5, radiating_dof="Surge")
    res = solver.solve(pb, keep_details=True)
    res_strict = solver_strict.solve(pb, keep_details=True)
    points = np.array([[2.0, 0.0, -1.0], [0.0, 3.0, -2.0]])
    assert np.abs(res.sources.imag).max() > 0
    for compute in ["compute_potential", "compute_velocity"]:
        actual = getattr(solver_strict, compute)(points, res_strict)
        assert np.abs(actual.imag).max() > 0  # Casting the sources to real used to discard the imaginary part
        np.testing.assert_allclose(actual, getattr(solver, compute)(points, res), rtol=1e-8)


#############################
#  Meshes with symmetries   #
#############################

@pytest.mark.parametrize("name", helpers.symmetric_meshes_of_single_panel().keys())
@pytest.mark.parametrize("linear_solver", ["lu_decomposition", "lu_decomposition_with_overwrite"])
def test_matrices_and_linear_solver_with_symmetries(name, linear_solver):
    sym_mesh = helpers.symmetric_meshes_of_single_panel()[name]
    params = dict(free_surface=0.0, water_depth=np.inf, wavenumber=1.0, diagonal_term_in_double_layer=True)
    S_ref, K_ref = cpt.DefaultMatrixEngine(green_function=helpers.green_function()).build_matrices(sym_mesh.merged(), sym_mesh.merged(), **params)

    engine = cpt.DefaultMatrixEngine(green_function=StrictArrayGreenFunction(), linear_solver=linear_solver)
    S, K = engine.build_matrices(sym_mesh, sym_mesh, **params)
    assert K.__array_namespace__() is xps
    np.testing.assert_allclose(np.array(S), S_ref, atol=1e-12)
    np.testing.assert_allclose(np.array(K), K_ref, atol=1e-12)

    rng = np.random.default_rng(0)
    b = rng.normal(size=K.shape[0]) + 1j*rng.normal(size=K.shape[0])
    x = engine.linear_solver(K, xps.asarray(b))
    assert isinstance(x, type(xps.asarray(b)))
    np.testing.assert_allclose(np.asarray(x), np.linalg.solve(K_ref, b), atol=1e-10)
    # The matrix-vector product is also available for the block matrices
    np.testing.assert_allclose(np.asarray(S @ xps.asarray(b)), S_ref @ b, atol=1e-10)


def symmetric_sphere_bodies():
    from capytaine.meshes import ReflectionSymmetricMesh, RotationSymmetricMesh
    parallelepiped = cpt.mesh_parallelepiped(resolution=(4, 4, 4), center=(0, 0, -2))
    half = parallelepiped.clipped(origin=(0, 0, 0), normal=(1, 0, 0))
    quarter = half.clipped(origin=(0, 0, 0), normal=(0, 1, 0))
    sphere_wedge = cpt.mesh_sphere(radius=1, resolution=(6, 6), center=(0, 0, -2)).extract_wedge(n=3)
    return {
        "reflection": ReflectionSymmetricMesh(half=half, plane="yOz"),
        "nested_reflections": ReflectionSymmetricMesh(ReflectionSymmetricMesh(half=quarter, plane="xOz"), plane="yOz"),
        "rotation": RotationSymmetricMesh(wedge=sphere_wedge, n=3),
    }


@pytest.mark.parametrize("name", ["reflection", "nested_reflections", "rotation"])
@pytest.mark.parametrize("linear_solver", ["lu_decomposition", "lu_decomposition_with_overwrite"])
@pytest.mark.parametrize("method", ["direct", "indirect"])
def test_solve_with_symmetries(name, linear_solver, method):
    sym_mesh = symmetric_sphere_bodies()[name]
    dofs = cpt.rigid_body_dofs(rotation_center=(0, 0, -2))
    body = cpt.FloatingBody(mesh=sym_mesh, dofs=dofs)
    ref_body = cpt.FloatingBody(mesh=sym_mesh.merged(), dofs=dofs)
    solver_strict = cpt.BEMSolver(method=method, engine=cpt.DefaultMatrixEngine(
        green_function=StrictArrayGreenFunction(), linear_solver=linear_solver))
    solver_ref = helpers.solver(method=method)
    for dof in ["Surge", "Surge"]:
        res_strict = solver_strict.solve(cpt.RadiationProblem(body=body, omega=1.0, radiating_dof=dof))
        res_ref = solver_ref.solve(cpt.RadiationProblem(body=ref_body, omega=1.0, radiating_dof=dof))
        assert res_strict.forces[dof] == pytest.approx(res_ref.forces[dof], rel=1e-6)


def test_float32_precision_with_symmetries():
    sym_mesh = helpers.symmetric_meshes_of_single_panel()["dihedral"]
    params = dict(free_surface=0.0, water_depth=np.inf, wavenumber=1.0, diagonal_term_in_double_layer=True)
    _, K_ref = cpt.DefaultMatrixEngine(green_function=helpers.green_function()).build_matrices(sym_mesh.merged(), sym_mesh.merged(), **params)
    engine = cpt.DefaultMatrixEngine(green_function=StrictArrayGreenFunction(floating_point_precision="float32"))
    _, K = engine.build_matrices(sym_mesh, sym_mesh, **params)
    assert K.dtype == xps.complex64
    b = np.ones(K.shape[0], dtype=np.complex64)
    x = engine.linear_solver(K, xps.asarray(b))
    assert x.dtype == xps.complex64
    np.testing.assert_allclose(np.asarray(x), np.linalg.solve(K_ref, b), rtol=1e-3, atol=1e-4)
