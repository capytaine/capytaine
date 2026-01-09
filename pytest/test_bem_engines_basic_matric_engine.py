from functools import lru_cache

import pytest

import numpy as np
import capytaine as cpt
from capytaine.new_meshes import Mesh, RotationSymmetricMesh, ReflectionSymmetricMesh
from capytaine.new_meshes.meshes import to_new_mesh


def test_engine_repr():
    r = repr(cpt.BasicMatrixEngine())
    assert "Delhommeau" in r


def test_cache_matrices():
    """Test how the BasicMatrixEngine caches the interaction matrices."""
    mesh = cpt.mesh_sphere(radius=1.0, resolution=(4, 3)).immersed_part()
    params_1 = dict(free_surface=0.0, water_depth=np.inf, wavenumber=1.0, adjoint_double_layer=True)
    params_2 = dict(free_surface=0.0, water_depth=np.inf, wavenumber=2.0, adjoint_double_layer=True)

    # No cache
    engine = cpt.BasicMatrixEngine()
    S, K             = engine._build_matrices_with_symmetries(mesh, mesh, **params_1)
    S_again, K_again = engine._build_matrices_with_symmetries(mesh, mesh, **params_1)
    assert S is not S_again
    assert K is not K_again

    # Cache
    engine = cpt.BasicMatrixEngine()
    S, K                     = engine._build_and_cache_matrices_with_symmetries(mesh, mesh, **params_1)
    S_again, K_again         = engine._build_and_cache_matrices_with_symmetries(mesh, mesh, **params_1)
    _, _                     = engine._build_and_cache_matrices_with_symmetries(mesh, mesh, **params_2)
    S_once_more, K_once_more = engine._build_and_cache_matrices_with_symmetries(mesh, mesh, **params_2)
    assert S is S_again
    assert S is not S_once_more
    assert K is K_again
    assert K is not K_once_more


@pytest.mark.parametrize("method", ['indirect', 'direct'])
def test_custom_linear_solver(method):
    """Solve a simple problem with a custom linear solver."""
    sphere = cpt.FloatingBody(mesh=cpt.mesh_sphere(radius=1.0, resolution=(4, 3)).immersed_part())
    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")
    problem = cpt.RadiationProblem(body=sphere, omega=1.0, water_depth=np.inf)

    reference_solver = cpt.BEMSolver(engine=cpt.BasicMatrixEngine(), method=method)
    reference_result = reference_solver.solve(problem)

    def my_linear_solver(A, b):
        """A dumb solver for testing."""
        return np.linalg.inv(A) @ b

    my_bem_solver = cpt.BEMSolver(
        engine=cpt.BasicMatrixEngine(linear_solver=my_linear_solver),
        method=method
    )
    assert 'my_linear_solver' in my_bem_solver.exportable_settings['linear_solver']

    result = my_bem_solver.solve(problem)
    assert np.isclose(reference_result.added_masses['Surge'], result.added_masses['Surge'])


@pytest.mark.parametrize("method", ['indirect', 'direct'])
def test_custom_linear_solver_returning_wrong_shape(method):
    sphere = cpt.FloatingBody(mesh=cpt.mesh_sphere(radius=1.0, resolution=(4, 3)).immersed_part())
    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")
    problem = cpt.RadiationProblem(body=sphere, omega=1.0, water_depth=np.inf)

    def my_linear_solver(A, b):
        """A dumb solver for testing."""
        return np.linalg.inv(A) @ b.reshape(-1, 1)

    my_bem_solver = cpt.BEMSolver(
        method=method,
        engine=cpt.BasicMatrixEngine(linear_solver=my_linear_solver)
    )
    with pytest.raises(ValueError):
        my_bem_solver.solve(problem)


def test_all_linear_solvers_work():
    """Test that all built-in linear solvers work i.e. "lu_decomposition", "lu_decomposition_with_overwrite" and "gmres"."""
    sphere = cpt.FloatingBody(mesh=cpt.mesh_sphere(radius=1.0, resolution=(4, 3)).immersed_part())
    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")
    problem = cpt.RadiationProblem(body=sphere, omega=1.0, water_depth=np.inf)

    bem_solver_lu = cpt.BEMSolver(
        engine=cpt.BasicMatrixEngine(linear_solver="lu_decomposition")
    )

    result_lu = bem_solver_lu.solve(problem)

    bem_solver_lu_with_overwrite = cpt.BEMSolver(
        engine=cpt.BasicMatrixEngine(linear_solver="lu_decomposition_with_overwrite")
    )

    result_lu_with_overwrite = bem_solver_lu_with_overwrite.solve(problem)

    bem_solver_gmres = cpt.BEMSolver(
        engine=cpt.BasicMatrixEngine(linear_solver="lu_decomposition")
    )

    result_gmres = bem_solver_gmres.solve(problem)

    assert np.isclose(
        result_lu.added_masses["Surge"], result_lu_with_overwrite.added_masses["Surge"]
    ) and np.isclose(
        result_lu.added_masses["Surge"], result_gmres.added_masses["Surge"]
    )


def test_lu_overwrite():
    K = np.array(np.random.rand(100,100), order='F')
    K_copy = np.copy(K)
    engine = cpt.BasicMatrixEngine(linear_solver="lu_decomposition_with_overwrite")
    engine.last_computed_matrices = (None, K)
    engine.linear_solver(K, np.random.rand(100))
    assert np.any(K_copy != K)


def test_lu_no_overwrite():
    K = np.array(np.random.rand(100,100), order='F')
    K_copy = np.copy(K)
    engine = cpt.BasicMatrixEngine(linear_solver="lu_decomposition")
    engine.last_computed_matrices = (None, K)
    engine.linear_solver(K, np.random.rand(100))
    assert np.all(K_copy == K)


def test_ram_estimation():
    sphere = cpt.FloatingBody(
        mesh=cpt.mesh_sphere(radius=1.0, resolution=(4, 3)).immersed_part()
    )
    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")
    problem = cpt.RadiationProblem(body=sphere, omega=1.0, water_depth=np.inf)
    nb_faces = problem.body.mesh.nb_faces
    green_function = cpt.green_functions.delhommeau.Delhommeau(
        floating_point_precision="float32"
    )

    reference_estimation = cpt.BasicMatrixEngine(
        linear_solver="gmres"
    ).compute_ram_estimation(problem)
    float_estimation = cpt.BasicMatrixEngine(
        linear_solver="gmres", green_function=green_function
    ).compute_ram_estimation(problem)
    lu_estimation = cpt.BasicMatrixEngine(
        linear_solver="lu_decomposition"
    ).compute_ram_estimation(problem)
    lu_and_float_estimation = cpt.BasicMatrixEngine(
        linear_solver="lu_decomposition", green_function=green_function
    ).compute_ram_estimation(problem)

    assert reference_estimation == nb_faces**2 * 16 * 2 / 1e9
    assert float_estimation == nb_faces**2 * 8 * 2 / 1e9
    assert lu_estimation == nb_faces**2 * 16 * 3 / 1e9
    assert lu_and_float_estimation == nb_faces**2 * 8 * 3 / 1e9


def test_ram_1_reflection_symmetry_estimation():
    reference_mesh = cpt.mesh_parallelepiped(resolution=(5, 5, 5), center =(0, 0, -2,))
    half_mesh = reference_mesh.clipped(cpt.Plane(point=(0, 0, 0), normal=(1, 0, 0)))
    mesh = ReflectionSymmetricMesh(half=to_new_mesh(half_mesh), plane="yOz")
    body = cpt.FloatingBody(
            mesh=mesh,
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.25,)),
            name=mesh.name,
            )
    problem = cpt.RadiationProblem(body=body, omega=1.0, water_depth=np.inf)
    nb_faces = problem.body.mesh.nb_faces

    green_function = cpt.green_functions.delhommeau.Delhommeau(
        floating_point_precision="float32"
    )

    gmres_estimation = cpt.BasicMatrixEngine(
        linear_solver="gmres"
    ).compute_ram_estimation(problem)
    gmres_estimation_float = cpt.BasicMatrixEngine(
        linear_solver="gmres",  green_function=green_function
    ).compute_ram_estimation(problem)
    lu_estimation = cpt.BasicMatrixEngine(
        linear_solver="lu_decomposition"
    ).compute_ram_estimation(problem)
    lu_overwrite_estimation = cpt.BasicMatrixEngine(
        linear_solver="lu_decomposition_with_overwrite"
    ).compute_ram_estimation(problem)

    assert gmres_estimation == nb_faces**2 * 2 * 16 / 1e9 / 2
    assert gmres_estimation_float == nb_faces**2 * 2 * 8 / 1e9 / 2
    assert round(lu_estimation,7) == nb_faces**2 * 16 * (2 * 1/2 + 1/2 + 1/2) / 1e9
    assert lu_overwrite_estimation == nb_faces**2 * 16 * (2 * 1/2 + 1/2) / 1e9


def test_ram_2_reflection_symmetries_estimation():
    reference_mesh = cpt.mesh_parallelepiped(resolution=(5, 5, 5), center =(0, 0, -2,))
    half_mesh = reference_mesh.clipped(cpt.Plane(point=(0, 0, 0), normal=(1, 0, 0)))
    quarter_mesh = half_mesh.clipped(cpt.Plane(point=(0, 0, 0), normal=(0, 1, 0)))
    mesh = ReflectionSymmetricMesh(half=ReflectionSymmetricMesh(half=to_new_mesh(quarter_mesh), plane="xOz"), plane="yOz")
    body = cpt.FloatingBody(
            mesh=mesh,
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.25,)),
            )
    problem = cpt.RadiationProblem(body=body, omega=1.0, water_depth=np.inf)
    nb_faces = problem.body.mesh.nb_faces

    green_function = cpt.green_functions.delhommeau.Delhommeau(
        floating_point_precision="float32"
    )

    gmres_estimation = cpt.BasicMatrixEngine(
        linear_solver="gmres"
    ).compute_ram_estimation(problem)
    gmres_estimation_float = cpt.BasicMatrixEngine(
        linear_solver="gmres",  green_function=green_function
    ).compute_ram_estimation(problem)
    lu_estimation = cpt.BasicMatrixEngine(
        linear_solver="lu_decomposition"
    ).compute_ram_estimation(problem)
    lu_overwrite_estimation = cpt.BasicMatrixEngine(
        linear_solver="lu_decomposition_with_overwrite"
    ).compute_ram_estimation(problem)

    assert gmres_estimation == nb_faces**2 * 2 * 16 / 1e9 / 4
    assert gmres_estimation_float == nb_faces**2 * 2 * 8 / 1e9 / 4
    assert lu_estimation == nb_faces**2 * 16 * (2 * 1/4 + 1/4 + 1/2) / 1e9
    assert lu_overwrite_estimation == nb_faces**2 * 16 * (2 * 1/4 + 1/2) / 1e9


@pytest.mark.parametrize("n", [2,3,4])
def test_ram_rotation_symmetries_estimation(n):
    reference_mesh = cpt.mesh_sphere(radius = 1, resolution=(5, 5), center=(0,0,-2))
    new_reference_mesh = to_new_mesh(reference_mesh)
    wedge = new_reference_mesh.extract_wedge(n=n)
    rotation_mesh = RotationSymmetricMesh(wedge=wedge, n=n)

    body = cpt.FloatingBody(
            mesh=rotation_mesh,
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.25,)),
            )
    problem = cpt.RadiationProblem(body=body, omega=1.0, water_depth=np.inf)
    nb_faces = problem.body.mesh.nb_faces

    green_function = cpt.green_functions.delhommeau.Delhommeau(
        floating_point_precision="float32"
    )

    gmres_estimation = cpt.BasicMatrixEngine(
        linear_solver="gmres"
    ).compute_ram_estimation(problem)
    gmres_estimation_float = cpt.BasicMatrixEngine(
        linear_solver="gmres",  green_function=green_function
    ).compute_ram_estimation(problem)
    lu_estimation = cpt.BasicMatrixEngine(
        linear_solver="lu_decomposition"
    ).compute_ram_estimation(problem)
    lu_overwrite_estimation = cpt.BasicMatrixEngine(
        linear_solver="lu_decomposition_with_overwrite"
    ).compute_ram_estimation(problem)

    assert gmres_estimation == nb_faces**2 * 2 * 16 / 1e9 / n
    assert gmres_estimation_float == nb_faces**2 * 2 * 8 / 1e9 / n
    assert lu_estimation == nb_faces**2 * 16 * (2 * 1/n + 1/n + 1/n) / 1e9
    assert lu_overwrite_estimation == nb_faces**2 * 16 * (2 * 1/n + 1/n) / 1e9


@lru_cache
def single_panel():
    vertices = np.array(
        [[0.5, 0.0, 0.0], [0.5, 0.0, -0.5], [0.5, 0.5, -0.3], [0.5, 0.5, -0.2]]
    )
    faces = np.array([[0, 1, 2, 3]])
    single_panel = Mesh(vertices=vertices, faces=faces)
    return single_panel


@pytest.mark.parametrize("sym_mesh", [
    ReflectionSymmetricMesh(ReflectionSymmetricMesh(single_panel(), plane="xOz"), plane="yOz"),
    RotationSymmetricMesh(single_panel(), n=4, axis='z+'),
    ], ids=["nested_reflections", "rotation+"])
def test_symmetry(sym_mesh):
    ref_mesh = sym_mesh.merged()
    engine = cpt.BasicMatrixEngine()
    params = dict(free_surface=0.0, water_depth=np.inf, wavenumber=1.0)
    S_ref, K_ref = engine.build_matrices(ref_mesh, ref_mesh, **params)
    S, K = engine.build_matrices(sym_mesh, sym_mesh, **params)
    assert np.allclose(np.array(S), S_ref)
    assert np.allclose(np.array(K), K_ref)