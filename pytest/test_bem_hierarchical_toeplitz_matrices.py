"""Tests for the resolution of the BEM problems using advanced techniques."""

import pytest

import numpy as np

import capytaine as cpt
from capytaine.meshes.meshes import Mesh
from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.meshes.symmetric import AxialSymmetricMesh, ReflectionSymmetricMesh, TranslationalSymmetricMesh

from capytaine.bodies import FloatingBody
from capytaine.bodies.predefined.spheres import Sphere
from capytaine.bodies.predefined.cylinders import HorizontalCylinder, VerticalCylinder

from capytaine.bem.problems_and_results import RadiationProblem
from capytaine.io.xarray import assemble_dataset

from capytaine.meshes.geometry import xOz_Plane, yOz_Plane

from capytaine.matrices.low_rank import LowRankMatrix

solver_with_sym = cpt.BEMSolver(engine=cpt.HierarchicalToeplitzMatrixEngine(ACA_distance=8, matrix_cache_size=0))
solver_without_sym = cpt.BEMSolver(engine=cpt.BasicMatrixEngine(matrix_cache_size=0))
method = ['indirect','direct']
adjoint_double_layer = [True,False]

# Use a single solver in the whole module to avoid reinitialisation of the solver (0.5 second).
# Do not use a matrix cache in order not to risk influencing a test with another.


@pytest.mark.parametrize("method", method)
@pytest.mark.parametrize("depth", [10.0, np.inf])
@pytest.mark.parametrize("omega", [0.1, 10.0])
def test_floating_sphere(depth, omega, method):
    """Comparison of the added mass and radiation damping
    for a heaving sphere described using several symmetries
    in finite and infinite depth.
    """
    reso = 2

    full_sphere = Sphere(radius=1.0, ntheta=reso, nphi=4*reso, axial_symmetry=False, clip_free_surface=True)
    full_sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=full_sphere, omega=omega, water_depth=depth)
    result1 = solver_with_sym.solve(problem,method=method)

    half_sphere_mesh = full_sphere.mesh.extract_faces(
        np.where(full_sphere.mesh.faces_centers[:, 1] > 0)[0],
        name="half_sphere_mesh")
    two_halves_sphere = FloatingBody(ReflectionSymmetricMesh(half_sphere_mesh, xOz_Plane))
    two_halves_sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=two_halves_sphere, omega=omega, water_depth=depth)
    result2 = solver_with_sym.solve(problem,method=method)

    quarter_sphere_mesh = half_sphere_mesh.extract_faces(
        np.where(half_sphere_mesh.faces_centers[:, 0] > 0)[0],
        name="quarter_sphere_mesh")
    four_quarter_sphere = FloatingBody(ReflectionSymmetricMesh(ReflectionSymmetricMesh(quarter_sphere_mesh, yOz_Plane), xOz_Plane))
    assert 'None' not in four_quarter_sphere.mesh.tree_view()
    four_quarter_sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=four_quarter_sphere, omega=omega, water_depth=depth)
    result3 = solver_with_sym.solve(problem,method=method)

    clever_sphere = Sphere(radius=1.0, ntheta=reso, nphi=4*reso, axial_symmetry=True, clip_free_surface=True)
    clever_sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=clever_sphere, omega=omega, water_depth=depth)
    result4 = solver_with_sym.solve(problem,method=method)

    # (quarter_sphere + half_sphere + full_sphere + clever_sphere).show()

    volume = 4/3*np.pi
    assert np.isclose(result1.added_masses["Heave"], result2.added_masses["Heave"], atol=1e-4*volume*problem.rho)
    assert np.isclose(result1.added_masses["Heave"], result3.added_masses["Heave"], atol=1e-4*volume*problem.rho)
    assert np.isclose(result1.added_masses["Heave"], result4.added_masses["Heave"], atol=1e-4*volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result2.radiation_dampings["Heave"], atol=1e-4*volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result3.radiation_dampings["Heave"], atol=1e-4*volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result4.radiation_dampings["Heave"], atol=1e-4*volume*problem.rho)


def test_two_vertical_cylinders():
    distance = 5

    buoy = VerticalCylinder(length=3, radius=0.5, center=(-distance/2, -1, 0), nx=8, nr=3, ntheta=8)
    buoy.mesh = buoy.mesh.merged()
    buoy.mesh = buoy.mesh.keep_immersed_part()
    buoy.add_translation_dof(name="Sway")

    two_buoys = FloatingBody.join_bodies(buoy, buoy.translated_x(distance))
    two_buoys.mesh = buoy.mesh.symmetrized(yOz_Plane)  # Use a ReflectionSymmetry as mesh...

    problems = [RadiationProblem(body=two_buoys, omega=1.0, radiating_dof=dof) for dof in two_buoys.dofs]

    results = assemble_dataset(solver_without_sym.solve_all(problems))
    # Check that the resulting matrix is symmetric
    assert np.isclose(results['added_mass'].data[0, 0, 0], results['added_mass'].data[0, 1, 1])
    assert np.isclose(results['added_mass'].data[0, 1, 0], results['added_mass'].data[0, 0, 1])
    assert np.isclose(results['radiation_damping'].data[0, 0, 0], results['radiation_damping'].data[0, 1, 1])
    assert np.isclose(results['radiation_damping'].data[0, 1, 0], results['radiation_damping'].data[0, 0, 1])

    results_with_sym = assemble_dataset(solver_with_sym.solve_all(problems))
    assert np.allclose(results['added_mass'].data, results_with_sym['added_mass'].data, rtol=1e-2)
    assert np.allclose(results['radiation_damping'].data, results_with_sym['radiation_damping'].data, rtol=1e-3)


@pytest.mark.parametrize("method", method)
def test_odd_axial_symmetry(method):
    """Buoy with odd number of slices."""
    def shape(z):
            return 0.1*(-(z+1)**2 + 16)
    buoy = FloatingBody(AxialSymmetricMesh.from_profile(shape, z_range=np.linspace(-5.0, 0.0, 9), nphi=5))
    buoy.add_translation_dof(direction=(0, 0, 1), name="Heave")

    problem = RadiationProblem(body=buoy, omega=2.0)
    result1 = solver_with_sym.solve(problem,method=method)

    full_buoy = FloatingBody(buoy.mesh.merged())
    full_buoy.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=full_buoy, omega=2.0)
    result2 = solver_with_sym.solve(problem,method=method)

    volume = buoy.mesh.volume
    assert np.isclose(result1.added_masses["Heave"], result2.added_masses["Heave"], atol=1e-4*volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result2.radiation_dampings["Heave"], atol=1e-4*volume*problem.rho)


@pytest.mark.parametrize("method", method)
@pytest.mark.parametrize("depth", [10.0, np.inf])
def test_horizontal_cylinder(depth,method):
    cylinder = HorizontalCylinder(length=10.0, radius=1.0, reflection_symmetry=False, translation_symmetry=False, nr=2, ntheta=10, nx=10)
    assert isinstance(cylinder.mesh, Mesh)
    cylinder.translate_z(-3.0)
    cylinder.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=cylinder, omega=1.0, water_depth=depth)
    result1 = solver_with_sym.solve(problem,method=method)

    trans_cylinder = HorizontalCylinder(length=10.0, radius=1.0, reflection_symmetry=False, translation_symmetry=True, nr=2, ntheta=10, nx=10)
    assert isinstance(trans_cylinder.mesh, CollectionOfMeshes)
    assert isinstance(trans_cylinder.mesh[0], TranslationalSymmetricMesh)
    trans_cylinder.translate_z(-3.0)
    trans_cylinder.add_translation_dof(direction=(0, 0, 1), name="Heave")
    problem = RadiationProblem(body=trans_cylinder, omega=1.0, water_depth=depth)
    result2 = solver_with_sym.solve(problem,method=method)

    # S, V = solver_with_sym.build_matrices(trans_cylinder.mesh, trans_cylinder.mesh)
    # S.plot_shape()

    assert np.isclose(result1.added_masses["Heave"], result2.added_masses["Heave"], atol=1e-4*cylinder.volume*problem.rho)
    assert np.isclose(result1.radiation_dampings["Heave"], result2.radiation_dampings["Heave"], atol=1e-4*cylinder.volume*problem.rho)


# HIERARCHICAL MATRICES
@pytest.mark.parametrize("method", method)
def test_low_rank_matrices(method):
    radius = 1.0
    resolution = 2
    perimeter = 2*np.pi*radius
    buoy = Sphere(radius=radius, center=(0.0, 0.0, 0.0),
                  ntheta=int(perimeter*resolution/2), nphi=int(perimeter*resolution),
                  clip_free_surface=True, axial_symmetry=False, name=f"buoy")
    buoy.add_translation_dof(name="Heave")
    two_distant_buoys = FloatingBody.join_bodies(buoy, buoy.translated_x(20))
    two_distant_buoys.mesh._meshes[1].name = "other_buoy_mesh"

    S, V = solver_with_sym.engine.build_matrices(two_distant_buoys.mesh, two_distant_buoys.mesh, 0.0, np.inf, 1.0, solver_with_sym.green_function)
    assert isinstance(S.all_blocks[0, 1], LowRankMatrix)
    assert isinstance(S.all_blocks[1, 0], LowRankMatrix)
    # S.plot_shape()

    problem = RadiationProblem(body=two_distant_buoys, omega=1.0, radiating_dof="buoy__Heave")
    result = solver_with_sym.solve(problem,method=method)
    result2 = solver_without_sym.solve(problem,method=method)

    assert np.isclose(result.added_masses['buoy__Heave'], result2.added_masses['buoy__Heave'], atol=10.0)
    assert np.isclose(result.radiation_dampings['buoy__Heave'], result2.radiation_dampings['buoy__Heave'], atol=10.0)


@pytest.mark.parametrize("method,adjoint_double_layer",zip(method,adjoint_double_layer))
def test_array_of_spheres(method,adjoint_double_layer):
    radius = 1.0
    resolution = 2
    perimeter = 2*np.pi*radius
    buoy = Sphere(radius=radius, center=(0.0, 0.0, 0.0),
                  ntheta=int(perimeter*resolution/2), nphi=int(perimeter*resolution),
                  clip_free_surface=True, axial_symmetry=False, name=f"buoy")
    buoy.add_translation_dof(name="Surge")
    buoy.add_translation_dof(name="Sway")
    buoy.add_translation_dof(name="Heave")

    # Corner case
    dumb_array = buoy.assemble_regular_array(distance=5.0, nb_bodies=(1, 1))
    assert dumb_array.mesh == buoy.mesh

    # Main case
    array = buoy.assemble_regular_array(distance=4.0, nb_bodies=(3, 3))

    assert isinstance(array.mesh, TranslationalSymmetricMesh)
    assert isinstance(array.mesh[0], TranslationalSymmetricMesh)
    assert array.mesh[0][0] == buoy.mesh

    assert len(array.dofs) == 3*3*3
    assert "2_0__Heave" in array.dofs

    #
    array = buoy.assemble_regular_array(distance=4.0, nb_bodies=(3, 1))

    fullS, fullV = solver_without_sym.engine.build_matrices(array.mesh, array.mesh, 0.0, np.inf,
			1.0, solver_without_sym.green_function, adjoint_double_layer=adjoint_double_layer)
    S, V = solver_with_sym.engine.build_matrices(array.mesh, array.mesh, 0.0, np.inf,
			1.0, solver_with_sym.green_function, adjoint_double_layer=adjoint_double_layer)

    assert isinstance(S, cpt.matrices.block.BlockMatrix)
    assert np.allclose(S.full_matrix(), fullS)
    assert np.allclose(V.full_matrix(), fullV)

    problem = RadiationProblem(body=array, omega=1.0, radiating_dof="2_0__Heave", water_depth=np.inf)

    result = solver_with_sym.solve(problem, method=method)
    result2 = solver_without_sym.solve(problem, method=method)

    assert np.isclose(result.added_masses['2_0__Heave'], result2.added_masses['2_0__Heave'], atol=15.0)
    assert np.isclose(result.radiation_dampings['2_0__Heave'], result2.radiation_dampings['2_0__Heave'], atol=15.0)
