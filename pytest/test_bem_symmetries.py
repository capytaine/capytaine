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
"""Resolution of BEM problems on meshes with symmetries."""

from functools import lru_cache

import pytest

import numpy as np
import capytaine as cpt
import capytaine_test_helpers as helpers
from capytaine.meshes import ReflectionSymmetricMesh
from capytaine.meshes.predefined import mesh_parallelepiped, mesh_rectangle, mesh_vertical_cylinder, mesh_horizontal_cylinder
from capytaine.tools.block_circulant_matrices import BlockCirculantMatrix, NestedBlockCirculantMatrix


@lru_cache
def symmetric_meshes():
    return {
        **helpers.symmetric_meshes(),
        # A sphere piercing the free surface, mirrored into two separate bodies
        "two_floating_spheres": ReflectionSymmetricMesh(cpt.mesh_sphere(center=(0, 2, 0)).immersed_part(), plane="xOz"),
    }


@pytest.mark.parametrize("name", symmetric_meshes().keys())
@pytest.mark.parametrize("linear_solver", ["lu_decomposition", "gmres"])
@pytest.mark.parametrize("method", ["direct", "indirect"])
def test_same_results_with_and_without_symmetries(name, linear_solver, method):
    sym_mesh = symmetric_meshes()[name]
    solver = cpt.BEMSolver(method=method, engine=cpt.DefaultMatrixEngine(
        green_function=helpers.green_function(), linear_solver=linear_solver))
    dofs = cpt.rigid_body_dofs(rotation_center=(0, 0, -1))
    results = {}
    for mesh in [sym_mesh, sym_mesh.merged()]:
        body = cpt.FloatingBody(mesh=mesh, dofs=dofs)
        problems = [cpt.RadiationProblem(body=body, omega=1.0, radiating_dof="Pitch"),
                    cpt.DiffractionProblem(body=body, omega=1.0, wave_direction=0.3)]
        results[mesh is sym_mesh] = solver.solve_all(problems)
    for res_sym, res_ref in zip(results[True], results[False]):
        for dof in ["Surge", "Pitch"]:
            assert res_sym.forces[dof] == pytest.approx(res_ref.forces[dof], rel=1e-4)


@pytest.mark.parametrize("name", ["nested_reflections", "rotation_4", "dihedral"])
def test_same_matrices_with_and_without_symmetries(name):
    sym_mesh = helpers.symmetric_meshes_of_single_panel()[name]
    ref_mesh = sym_mesh.merged()
    engine = cpt.DefaultMatrixEngine()
    params = dict(free_surface=0.0, water_depth=np.inf, wavenumber=1.0, diagonal_term_in_double_layer=True)
    S_ref, K_ref = engine.build_matrices(ref_mesh, ref_mesh, **params)
    S, K = engine.build_matrices(sym_mesh, sym_mesh, **params)
    np.testing.assert_allclose(np.array(S), S_ref, atol=1e-12)
    np.testing.assert_allclose(np.array(K), K_ref, atol=1e-12)


def test_same_post_processing_matrices_with_and_without_symmetries():
    engine = cpt.DefaultMatrixEngine()
    sym_mesh = symmetric_meshes()["two_floating_spheres"]
    mesh = sym_mesh.merged()
    assert np.allclose(engine.build_S_matrix(sym_mesh.faces_centers, sym_mesh), engine.build_S_matrix(mesh.faces_centers, mesh))
    assert np.allclose(engine.build_fullK_matrix(sym_mesh.faces_centers, sym_mesh), engine.build_fullK_matrix(mesh.faces_centers, mesh))


@pytest.mark.parametrize("name", helpers.symmetric_meshes_of_single_panel().keys())
def test_gmres_linear_solver_with_symmetries(name):
    sym_mesh = helpers.symmetric_meshes_of_single_panel()[name]
    params = dict(free_surface=0.0, water_depth=np.inf, wavenumber=1.0, diagonal_term_in_double_layer=True)
    engine = cpt.DefaultMatrixEngine(green_function=helpers.green_function(), linear_solver="gmres")
    _, K = engine.build_matrices(sym_mesh, sym_mesh, **params)
    K_ref = np.array(K)
    rng = np.random.default_rng(0)
    b = rng.normal(size=K.shape[0]) + 1j*rng.normal(size=K.shape[0])
    np.testing.assert_allclose(engine.linear_solver(K, b), np.linalg.solve(K_ref, b), atol=1e-5)


def test_lid_with_plane_symmetry():
    mesh = mesh_horizontal_cylinder(reflection_symmetry=True).immersed_part()
    lid_mesh = ReflectionSymmetricMesh(
            mesh_rectangle(size=(1.0, 10,), faces_max_radius=0.5, center=(0, 0.5, -0.05,)),
            plane=mesh.plane
            )
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, radiating_dof="Heave")
    solver = helpers.solver()
    S, K = solver.engine.build_matrices(pb.body.mesh_including_lid, pb.body.mesh_including_lid,
                                        free_surface=pb.free_surface, water_depth=pb.water_depth,
                                        wavenumber=pb.wavenumber)
    assert isinstance(S, BlockCirculantMatrix)
    assert isinstance(K, BlockCirculantMatrix)

    res = solver.solve(pb)

    ref_body = cpt.FloatingBody(mesh=mesh.merged(), lid_mesh=lid_mesh.merged(), dofs=cpt.rigid_body_dofs())
    ref_pb = cpt.RadiationProblem(body=ref_body, wavelength=1.0, radiating_dof="Heave")
    ref_res = solver.solve(ref_pb)
    assert res.force["Heave"] == pytest.approx(ref_res.force["Heave"])


def test_lid_with_nested_plane_symmetry():
    mesh = mesh_parallelepiped(center=(0, 0, -1.0), reflection_symmetry=True)
    hull_mesh, lid_mesh = mesh.extract_lid()
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, radiating_dof="Heave")
    solver = helpers.solver()
    S, K = solver.engine.build_matrices(pb.body.mesh_including_lid, pb.body.mesh_including_lid,
                                        free_surface=pb.free_surface, water_depth=pb.water_depth,
                                        wavenumber=pb.wavenumber)
    assert isinstance(S, NestedBlockCirculantMatrix)
    assert isinstance(K, NestedBlockCirculantMatrix)

    res = solver.solve(pb)

    ref_body = cpt.FloatingBody(mesh=mesh.merged(), lid_mesh=lid_mesh.merged(), dofs=cpt.rigid_body_dofs())
    ref_pb = cpt.RadiationProblem(body=ref_body, wavelength=1.0, radiating_dof="Heave")
    ref_res = solver.solve(ref_pb)
    assert res.force["Heave"] == pytest.approx(ref_res.force["Heave"])


def test_lid_with_rotation_symmetry():
    mesh = mesh_vertical_cylinder(length=2.0, center=(0, 0, -1.0), axial_symmetry=True)
    hull_mesh, lid_mesh = mesh.extract_lid()
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    pb = cpt.RadiationProblem(body=body, wavelength=1.0, radiating_dof="Heave")
    solver = helpers.solver()
    S, K = solver.engine.build_matrices(pb.body.mesh_including_lid, pb.body.mesh_including_lid,
                                        free_surface=pb.free_surface, water_depth=pb.water_depth,
                                        wavenumber=pb.wavenumber)
    assert isinstance(S, BlockCirculantMatrix)
    assert isinstance(K, BlockCirculantMatrix)
    res = solver.solve(pb)

    ref_body = cpt.FloatingBody(mesh=mesh.merged(), lid_mesh=lid_mesh.merged(), dofs=cpt.rigid_body_dofs())
    ref_pb = cpt.RadiationProblem(body=ref_body, wavelength=1.0, radiating_dof="Heave")
    ref_res = solver.solve(ref_pb)
    assert res.force["Heave"] == pytest.approx(ref_res.force["Heave"])
