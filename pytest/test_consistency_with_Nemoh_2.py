"""Quantitatively compare the results of Capytaine with the results from Nemoh 2."""

import pytest
import numpy as np
import capytaine as cpt

from capytaine.bodies.predefined.spheres import Sphere
from capytaine.bodies.predefined.cylinders import HorizontalCylinder
from capytaine.post_pro.free_surfaces import FreeSurface
from capytaine.post_pro.kochin import compute_kochin


@pytest.fixture
def nemoh2_solver():
    gf = cpt.Delhommeau(
            tabulation_nr=328, tabulation_nz=46,
            tabulation_grid_shape='legacy', tabulation_nb_integration_points=251,
            gf_singularities="high_freq",
            )
    solver = cpt.BEMSolver(
            engine=cpt.BasicMatrixEngine(matrix_cache_size=0),
            green_function=gf
            )
    return solver


def test_immersed_sphere(nemoh2_solver):
    """Compare with Nemoh 2.0 for a sphere in infinite fluid.

    The test is ran for two degrees of freedom; due to the symmetries of the problem, the results should be the same.
    They are actually slightly different due to the meshing of the sphere.
    """
    sphere = Sphere(radius=1.0, ntheta=10, nphi=40, clip_free_surface=False)
    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")
    sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")

    problem = cpt.RadiationProblem(body=sphere, radiating_dof="Heave", free_surface=np.inf, water_depth=np.inf)
    result = nemoh2_solver.solve(problem)
    assert np.isclose(result.added_masses["Heave"],       2187, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(result.added_masses["Surge"],        0.0, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(result.radiation_dampings["Heave"],  0.0, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(result.radiation_dampings["Surge"],  0.0, atol=1e-3*sphere.volume*problem.rho)

    problem = cpt.RadiationProblem(body=sphere, radiating_dof="Surge", free_surface=np.inf, water_depth=np.inf)
    result = nemoh2_solver.solve(problem)
    assert np.isclose(result.added_masses["Surge"],       2194, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(result.added_masses["Heave"],        0.0, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(result.radiation_dampings["Surge"],  0.0, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(result.radiation_dampings["Heave"],  0.0, atol=1e-3*sphere.volume*problem.rho)


def test_build_matrix_of_rankine_and_reflected_rankine(nemoh2_solver):
    gf = nemoh2_solver.green_function
    sphere = Sphere(radius=1.0, ntheta=2, nphi=3, clip_free_surface=True)

    S, V = gf.evaluate(sphere.mesh, sphere.mesh, 0.0, np.inf, 0.0)
    S_ref = np.array([[-0.15413386, -0.21852682, -0.06509213, -0.16718431, -0.06509213, -0.16718431],
                      [-0.05898834, -0.39245688, -0.04606661, -0.18264734, -0.04606661, -0.18264734],
                      [-0.06509213, -0.16718431, -0.15413386, -0.21852682, -0.06509213, -0.16718431],
                      [-0.04606661, -0.18264734, -0.05898834, -0.39245688, -0.04606661, -0.18264734],
                      [-0.06509213, -0.16718431, -0.06509213, -0.16718431, -0.15413386, -0.21852682],
                      [-0.04606661, -0.18264734, -0.04606661, -0.18264734, -0.05898834, -0.39245688]])
    assert np.allclose(S, S_ref)

    S, V = gf.evaluate(sphere.mesh, sphere.mesh, 0.0, np.inf, np.inf)
    S_ref = np.array([[-0.12666269, -0.07804937, -0.03845837, -0.03993999, -0.03845837, -0.03993999],
                      [-0.02106031, -0.16464793, -0.01169102, -0.02315146, -0.01169102, -0.02315146],
                      [-0.03845837, -0.03993999, -0.12666269, -0.07804937, -0.03845837, -0.03993999],
                      [-0.01169102, -0.02315146, -0.02106031, -0.16464793, -0.01169102, -0.02315146],
                      [-0.03845837, -0.03993999, -0.03845837, -0.03993999, -0.12666269, -0.07804937],
                      [-0.01169102, -0.02315146, -0.01169102, -0.02315146, -0.02106031, -0.16464793]])
    assert np.allclose(S, S_ref)


def test_floating_sphere_finite_freq(nemoh2_solver):
    """Compare with Nemoh 2.0 for some cases of a heaving sphere at the free surface in infinite depth."""
    sphere = Sphere(radius=1.0, ntheta=3, nphi=12, clip_free_surface=True)
    sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")

    # omega = 1, radiation
    problem = cpt.RadiationProblem(body=sphere, omega=1.0, water_depth=np.inf)
    result = nemoh2_solver.solve(problem, keep_details=True)
    assert np.isclose(result.added_masses["Heave"],       1819.6, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(result.radiation_dampings["Heave"], 379.39, atol=1e-3*sphere.volume*problem.rho)

    # omega = 1, free surface
    grid = np.meshgrid(np.linspace(-50.0, 50.0, 5), np.linspace(-50.0, 50.0, 5))
    eta = nemoh2_solver.compute_free_surface_elevation(grid, result)
    ref = np.array(
            [[-0.4340802E-02-0.4742809E-03j, -0.7986111E-03+0.4840984E-02j, 0.2214827E-02+0.4700642E-02j, -0.7986111E-03+0.4840984E-02j, -0.4340803E-02-0.4742807E-03j],
             [-0.7986111E-03+0.4840984E-02j, 0.5733187E-02-0.2179381E-02j, 0.9460892E-03-0.7079404E-02j, 0.5733186E-02-0.2179381E-02j, -0.7986110E-03+0.4840984E-02j],
             [0.2214827E-02+0.4700643E-02j, 0.9460892E-03-0.7079403E-02j, -0.1381670E-01+0.6039315E-01j, 0.9460892E-03-0.7079405E-02j, 0.2214827E-02+0.4700643E-02j],
             [-0.7986111E-03+0.4840984E-02j, 0.5733186E-02-0.2179381E-02j, 0.9460891E-03-0.7079404E-02j, 0.5733187E-02-0.2179380E-02j, -0.7986113E-03+0.4840984E-02j],
             [-0.4340803E-02-0.4742807E-03j, -0.7986111E-03+0.4840984E-02j, 0.2214827E-02+0.4700643E-02j, -0.7986113E-03+0.4840983E-02j, -0.4340803E-02-0.4742809E-03j]]
        )
    assert np.allclose(eta/(-1j*problem.omega), ref, rtol=1e-2)

    # omega = 1, diffraction
    problem = cpt.DiffractionProblem(body=sphere, omega=1.0, water_depth=np.inf)
    result = nemoh2_solver.solve(problem, keep_details=True)
    assert np.isclose(result.forces["Heave"], 1834.9 * np.exp(-2.933j), rtol=1e-3)

    # omega = 1, Kochin function of diffraction problem

    kochin = compute_kochin(result, np.linspace(0, np.pi, 10))

    ref_kochin = np.array([
        0.20229*np.exp(-1.5872j), 0.20369*np.exp(-1.5871j),
        0.20767*np.exp(-1.5868j), 0.21382*np.exp(-1.5863j),
        0.22132*np.exp(-1.5857j), 0.22931*np.exp(-1.5852j),
        0.23680*np.exp(-1.5847j), 0.24291*np.exp(-1.5843j),
        0.24688*np.exp(-1.5841j), 0.24825*np.exp(-1.5840j),
        ])

    assert np.allclose(kochin, ref_kochin, rtol=1e-3)

    # omega = 2, radiation
    problem = cpt.RadiationProblem(body=sphere, omega=2.0, water_depth=np.inf)
    result = nemoh2_solver.solve(problem)
    assert np.isclose(result.added_masses["Heave"],       1369.3, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(result.radiation_dampings["Heave"], 1425.6, atol=1e-3*sphere.volume*problem.rho)

    # omega = 2, diffraction
    problem = cpt.DiffractionProblem(body=sphere, omega=2.0, water_depth=np.inf)
    result = nemoh2_solver.solve(problem)
    assert np.isclose(result.forces["Heave"], 5846.6 * np.exp(-2.623j), rtol=1e-3)


def test_alien_sphere(nemoh2_solver):
    """Compare with Nemoh 2.0 for some cases of a heaving sphere at the free surface in infinite depth
    for a non-usual gravity and density."""
    sphere = Sphere(radius=1.0, ntheta=3, nphi=12, clip_free_surface=True)
    sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")

    # radiation
    problem = cpt.RadiationProblem(body=sphere, rho=450.0, g=1.625, omega=1.0, radiating_dof="Heave")
    result = nemoh2_solver.solve(problem)
    assert np.isclose(result.added_masses["Heave"],       515, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(result.radiation_dampings["Heave"], 309, atol=1e-3*sphere.volume*problem.rho)

    # diffraction
    problem = cpt.DiffractionProblem(body=sphere, rho=450.0, g=1.625, omega=1.0, water_depth=np.inf)
    result = nemoh2_solver.solve(problem)
    assert np.isclose(result.forces["Heave"], 548.5 * np.exp(-2.521j), rtol=1e-2)


def test_floating_sphere_finite_depth(nemoh2_solver):
    """Compare with Nemoh 2.0 for some cases of a heaving sphere at the free surface in finite depth."""
    sphere = Sphere(radius=1.0, ntheta=3, nphi=12, clip_free_surface=True)
    sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")

    # omega = 1, radiation
    problem = cpt.RadiationProblem(body=sphere, omega=1.0, radiating_dof="Heave", water_depth=10.0)
    result = nemoh2_solver.solve(problem, keep_details=True)
    assert np.isclose(result.added_masses["Heave"],       1740.6, atol=1e-3*sphere.volume*problem.rho)
    assert np.isclose(result.radiation_dampings["Heave"], 380.46, atol=1e-3*sphere.volume*problem.rho)

    kochin = compute_kochin(result, np.linspace(0, np.pi, 3))
    assert np.allclose(kochin, np.roll(kochin, 1))  # The far field is the same in all directions.
    assert np.isclose(kochin[0]/(-1j*problem.omega), -0.2267+3.49e-3j, rtol=1e-3)

    # omega = 1, diffraction
    problem = cpt.DiffractionProblem(body=sphere, omega=1.0, wave_direction=0.0, water_depth=10.0)
    result = nemoh2_solver.solve(problem)
    assert np.isclose(result.forces["Heave"], 1749.4 * np.exp(-2.922j), rtol=1e-3)

    # omega = 2, radiation
    problem = cpt.RadiationProblem(body=sphere, omega=2.0, radiating_dof="Heave", water_depth=10.0)
    result = nemoh2_solver.solve(problem)
    assert np.isclose(result.added_masses["Heave"],       1375.0, atol=1e-2*sphere.volume*problem.rho)
    assert np.isclose(result.radiation_dampings["Heave"], 1418.0, atol=1e-2*sphere.volume*problem.rho)

    # omega = 2, diffraction
    problem = cpt.DiffractionProblem(body=sphere, omega=2.0, wave_direction=0.0, water_depth=10.0)
    result = nemoh2_solver.solve(problem)
    assert np.isclose(result.forces["Heave"], 5872.8 * np.exp(-2.627j), rtol=1e-2)


def test_two_distant_spheres_in_finite_depth(nemoh2_solver):
    radius = 0.5
    resolution = 4
    perimeter = 2*np.pi*radius
    buoy = Sphere(radius=radius, center=(0.0, 0.0, 0.0),
                  ntheta=int(perimeter*resolution/2), nphi=int(perimeter*resolution),
                  clip_free_surface=True, axial_symmetry=False, name="buoy")
    other_buoy = buoy.translated_x(20, name="other_buoy")
    both_buoys = buoy.join_bodies(other_buoy)
    both_buoys.add_translation_dof(name="Surge")
    problem = cpt.RadiationProblem(body=both_buoys, radiating_dof="Surge", water_depth=10, omega=7.0)
    result = nemoh2_solver.solve(problem)

    total_volume = 2*4/3*np.pi*radius**3
    assert np.isclose(result.added_masses['Surge'], 124.0, atol=1e-3*total_volume*problem.rho)
    assert np.isclose(result.radiation_dampings['Surge'], 913.3, atol=1e-2*total_volume*problem.rho)


def test_multibody(nemoh2_solver):
    """Compare with Nemoh 2.0 for two bodies."""
    sphere = Sphere(radius=1.0, ntheta=5, nphi=20)
    sphere.translate_z(-2.0)
    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")
    sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")

    cylinder = HorizontalCylinder(length=5.0, radius=1.0,
                                  nx=10, nr=1, ntheta=10)
    cylinder.translate([+1.5, 3.0, -3.0])
    cylinder.add_translation_dof(direction=(1, 0, 0), name="Surge")
    cylinder.add_translation_dof(direction=(0, 0, 1), name="Heave")

    both = cylinder + sphere
    total_volume = cylinder.volume + sphere.volume
    # both.show()

    problems = [cpt.RadiationProblem(body=both, radiating_dof=dof, omega=1.0) for dof in both.dofs]
    problems += [cpt.DiffractionProblem(body=both, wave_direction=0.0, omega=1.0)]
    results = [nemoh2_solver.solve(problem) for problem in problems]
    data = cpt.assemble_dataset(results)

    data_from_nemoh_2 = np.array([
        [3961.86548, 50.0367661, -3.32347107, 6.36901855E-02, 172.704819, 19.2018471, -5.67303181, -2.98873377],
        [-3.08301544, 5.72392941E-02, 14522.1689, 271.796814, 128.413834, 6.03351116, 427.167358, 64.1587067],
        [161.125534, 17.8332844, 126.392113, 5.88006783, 2242.47412, 7.17850924, 1.29002571, 0.393169671],
        [-5.02560759, -2.75930357, 419.927460, 63.3179016, 1.23501396, 0.416424811, 2341.57593, 15.8266096],
    ])

    dofs_names = list(both.dofs.keys())
    assert np.allclose(
        data['added_mass'].sel(omega=1.0, radiating_dof=dofs_names, influenced_dof=dofs_names).values,
        data_from_nemoh_2[:, ::2],
        atol=1e-3*total_volume*problems[0].rho
    )
    assert np.allclose(
        data['radiation_damping'].sel(omega=1.0, radiating_dof=dofs_names, influenced_dof=dofs_names).values,
        data_from_nemoh_2[:, 1::2],
        atol=1e-3*total_volume*problems[0].rho
    )
