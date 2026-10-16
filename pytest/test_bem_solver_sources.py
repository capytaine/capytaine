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
"""Evaluation of fields from bare sources, batched over several source distributions."""
import numpy as np

import capytaine as cpt
import capytaine_test_helpers as helpers
from capytaine.bem.engines import DefaultMatrixEngine
from capytaine.meshes.symmetric_meshes import ReflectionSymmetricMesh
from capytaine.tools.lazy_matrices import LazyMatrix


class CountingEngine(DefaultMatrixEngine):
    """Engine counting the calls to the matrix builders used in post-processing."""
    def __init__(self):
        super().__init__(green_function=helpers.green_function())
        self.nb_S_builds = 0
        self.nb_fullK_builds = 0

    def build_S_matrix(self, *args, **kwargs):
        self.nb_S_builds += 1
        return super().build_S_matrix(*args, **kwargs)

    def build_fullK_matrix(self, *args, **kwargs):
        self.nb_fullK_builds += 1
        return super().build_fullK_matrix(*args, **kwargs)


def symmetric_body():
    half = cpt.mesh_sphere(radius=1.0, resolution=(4, 4), center=(0, 1.5, 0)).immersed_part()
    mesh = ReflectionSymmetricMesh(half, plane="xOz")
    return cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(only=("Heave",)))


def bodies():
    return {
        "no lid": helpers.small_sphere_body(("Surge", "Heave")),
        "lid": helpers.small_sphere_body(("Surge", "Heave"), lid=True),
        "symmetric": symmetric_body(),
    }


def solve_several(solver, body, water_depth):
    problems = [cpt.DiffractionProblem(body=body, wave_direction=beta, omega=1.3, water_depth=water_depth) for beta in (0.0, 0.7, 2.0)]
    problems += [cpt.RadiationProblem(body=body, radiating_dof=dof, omega=1.3, water_depth=water_depth) for dof in body.dofs]
    return solver.solve_all(problems, keep_details=True)


def gf_params_of(result):
    return dict(free_surface=result.free_surface, water_depth=result.water_depth, wavenumber=result.encounter_wavenumber)


def test_batched_evaluation_matches_per_result_evaluation():
    points = np.array([[0.3, 0.2, -1.5], [2.0, 1.0, -0.5], [-1.0, 0.5, -0.2], [0.0, 3.0, -2.0]])
    for name, body in bodies().items():
        for water_depth in [np.inf, 5.0]:
            solver = helpers.solver()
            results = solve_several(solver, body, water_depth)
            mesh = results[0].body.mesh_including_lid
            sources = np.stack([res.sources for res in results])  # (nb_problems, nb_faces)
            gf_params = gf_params_of(results[0])

            potential = solver._compute_potential_from_sources(points, mesh, sources, **gf_params)
            assert potential.shape == (len(results), 4)
            gradient = solver._compute_potential_gradient_from_sources(points, mesh, sources, **gf_params)
            assert gradient.shape == (len(results), 4, 3)

            for i, res in enumerate(results):
                np.testing.assert_allclose(potential[i], solver.compute_potential(points, res), rtol=1e-10, atol=1e-12, err_msg=f"{name}, depth {water_depth}")
                np.testing.assert_allclose(gradient[i], solver._compute_potential_gradient(points, res), rtol=1e-10, atol=1e-12, err_msg=f"{name}, depth {water_depth}")


def test_extra_batch_axes_and_points_shapes():
    body = helpers.small_sphere_body(("Surge",), lid=True)
    solver = helpers.solver()
    res = solver.solve(cpt.DiffractionProblem(body=body, wave_direction=0.3, omega=1.0), keep_details=True)
    mesh = res.body.mesh_including_lid
    gf_params = gf_params_of(res)
    factors = np.array([1.0, 2.0, -1j])
    sources = factors[:, None, None] * np.stack([res.sources, res.sources])[None, :, :]  # (3, 2, nb_faces)
    X, Y, Z = np.meshgrid(np.linspace(-2, 2, 3), np.linspace(-2, 2, 4), [-1.0, -0.5])
    points = np.array([X, Y, Z])
    potential = solver._compute_potential_from_sources(points, mesh, sources, **gf_params)
    assert potential.shape == (3, 2, 4, 3, 2)
    gradient = solver._compute_potential_gradient_from_sources(points, mesh, sources, **gf_params)
    assert gradient.shape == (3, 2, 4, 3, 2, 3)
    ref = solver.compute_potential(points, res)
    np.testing.assert_allclose(potential[1, 1], 2*ref, rtol=1e-10)
    np.testing.assert_allclose(potential[2, 0], -1j*ref, rtol=1e-10)
    # Single distribution: no batch axis
    assert solver._compute_potential_from_sources(points, mesh, res.sources, **gf_params).shape == ref.shape
    # Points given as a mesh
    assert solver._compute_potential_gradient_from_sources(body.mesh, mesh, sources, **gf_params).shape == (3, 2, body.mesh.nb_faces, 3)


def test_public_methods_unchanged():
    body = helpers.small_sphere_body(("Heave",), lid=True)
    solver = helpers.solver()
    res = solver.solve(cpt.RadiationProblem(body=body, radiating_dof="Heave", omega=1.0), keep_details=True)
    points = np.array([[0.5, 0.5, -2.0], [1.5, 0.0, -0.3]])
    mesh = res.body.mesh_including_lid
    S = np.asarray(solver.engine.build_S_matrix(points, mesh, **gf_params_of(res)))
    np.testing.assert_allclose(solver.compute_potential(points, res), S @ res.sources, rtol=1e-12)
    np.testing.assert_allclose(solver.compute_pressure(points, res), 1j*res.omega*res.rho*(S @ res.sources), rtol=1e-12)
    assert solver.compute_velocity(points, res).shape == (2, 3)
    assert solver.compute_free_surface_elevation(points[:, :2], res).shape == (2,)
    assert solver.compute_potential(points[0], res).shape == (1,)


def test_batched_evaluation_builds_matrix_once():
    body = helpers.small_sphere_body(("Surge", "Heave"), lid=True)
    engine = CountingEngine()
    solver = cpt.BEMSolver(engine=engine)
    results = solve_several(solver, body, np.inf)
    mesh = results[0].body.mesh_including_lid
    sources = np.stack([res.sources for res in results])
    points = body.mesh.faces_centers

    solver._compute_potential_from_sources(points, mesh, sources, **gf_params_of(results[0]))
    assert (engine.nb_S_builds, engine.nb_fullK_builds) == (1, 0)
    solver._compute_potential_gradient_from_sources(points, mesh, sources, **gf_params_of(results[0]))
    assert (engine.nb_S_builds, engine.nb_fullK_builds) == (1, 1)

    for res in results:  # One build per result with the public methods
        solver.compute_velocity(points, res)
    assert engine.nb_fullK_builds == 1 + len(results)


def test_lazy_matrix_with_2d_right_hand_side():
    full = np.random.default_rng(0).normal(size=(30, 7)) + 0j
    lazy = LazyMatrix(lambda sl: full[sl, :], (30, 7), chunk_size=4, dtype=complex)
    rhs = np.random.default_rng(1).normal(size=(7, 3)) + 1j
    np.testing.assert_allclose(lazy @ rhs, full @ rhs)
    np.testing.assert_allclose(lazy @ rhs[:, 0], full @ rhs[:, 0])


def test_many_points_use_lazy_matrix_and_batching():
    body = helpers.small_sphere_body(("Surge",))
    solver = helpers.solver()
    results = [solver.solve(cpt.DiffractionProblem(body=body, wave_direction=beta, omega=1.0), keep_details=True) for beta in (0.0, 1.0)]
    rng = np.random.default_rng(0)
    points = np.stack([rng.uniform(-3, 3, 600), rng.uniform(-3, 3, 600), rng.uniform(-3, -0.2, 600)], axis=1)
    mesh = results[0].body.mesh_including_lid
    sources = np.stack([res.sources for res in results])
    assert isinstance(solver.engine.build_S_matrix(points, mesh, **gf_params_of(results[0])), LazyMatrix)
    potential = solver._compute_potential_from_sources(points, mesh, sources, **gf_params_of(results[0]))
    for i, res in enumerate(results):
        np.testing.assert_allclose(potential[i], solver.compute_potential(points, res), rtol=1e-10)


def test_zero_and_infinite_frequency_batched():
    body = helpers.small_sphere_body(("Heave",))
    solver = helpers.solver()
    for omega in [0.0, np.inf]:
        res = solver.solve(cpt.RadiationProblem(body=body, radiating_dof="Heave", omega=omega), keep_details=True)
        points = np.array([[0.5, 0.5, -2.0], [1.5, 0.0, -0.3]])
        mesh = res.body.mesh_including_lid
        batched = solver._compute_potential_from_sources(points, mesh, res.sources.reshape(1, -1), **gf_params_of(res))
        single = solver.compute_potential(points, res)
        assert batched.symbol == single.symbol  # The result stays symbolic
        assert batched.value.shape == (1, 2)
        np.testing.assert_allclose(batched.value[0], single.value, rtol=1e-10)
        gradient = solver._compute_potential_gradient_from_sources(points, mesh, res.sources.reshape(1, -1), **gf_params_of(res))
        assert gradient.symbol == single.symbol
        np.testing.assert_allclose(gradient.value[0], solver._compute_potential_gradient(points, res).value, rtol=1e-10)
