import os
import logging

import numpy as np
import pytest

import capytaine as cpt
from capytaine.bem.problems_and_results import LinearPotentialFlowProblem, \
        LinearPotentialFlowResult, DiffractionResult, RadiationResult
from capytaine.io.legacy import import_cal_file


@pytest.fixture
def sphere():
    sphere = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(0, 0, -2), radius=1.0, resolution=(10, 20)),
            name="sphere",
            )
    sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    return sphere


@pytest.fixture
def solver():
    solver = cpt.BEMSolver()
    return solver


def test_LinearPotentialFlowProblem(sphere):
    # Without a body
    pb = LinearPotentialFlowProblem(omega=1.0)
    assert pb.omega == 1.0
    assert pb.period == 2*np.pi
    assert pb.wavenumber == 1.0/9.81
    assert pb.wavelength == 9.81*2*np.pi

    assert LinearPotentialFlowProblem(free_surface=np.inf, water_depth=np.inf).water_depth == np.inf
    assert LinearPotentialFlowProblem(free_surface=0.0, water_depth=np.inf).water_depth == np.inf

    pb = LinearPotentialFlowProblem(free_surface=0.0, water_depth=1.0, omega=1.0)
    assert pb.water_depth == 1.0
    assert np.isclose(pb.omega**2, pb.g*pb.wavenumber*np.tanh(pb.wavenumber*pb.water_depth))

    with pytest.raises(NotImplementedError):
        LinearPotentialFlowProblem(free_surface=2.0)

    with pytest.raises(NotImplementedError):
        LinearPotentialFlowProblem(free_surface=np.inf, water_depth=2.0)

    with pytest.raises(ValueError):
        LinearPotentialFlowProblem(free_surface=0.0, water_depth=-1.0)

    with pytest.raises(TypeError):
        LinearPotentialFlowProblem(radiating_dof="Heave")

    with pytest.raises(ValueError):
        LinearPotentialFlowProblem(body=cpt.FloatingBody(mesh=cpt.Mesh([], [])))

    # With a body
    pb = LinearPotentialFlowProblem(body=sphere, omega=1.0)
    pb.boundary_condition = sphere.mesh.faces_normals @ (1, 1, 1)
    assert list(pb.influenced_dofs.keys()) == ['Heave']

    pb2 = LinearPotentialFlowProblem(body=sphere, omega=2.0)
    pb2.boundary_condition = sphere.mesh.faces_normals @ (1, 1, 1)
    assert pb < pb2

    # Test transformation to result class
    res = pb.make_results_container()
    assert isinstance(res, LinearPotentialFlowResult)
    assert res.problem is pb
    assert res.omega == pb.omega
    assert res.period == pb.period
    assert res.body is pb.body


def test_mesh_inconsistent_with_dofs(sphere):
    n = sphere.mesh.nb_faces
    sphere.dofs["Wrong_dof"] = np.ones((n//2, 3))
    with pytest.raises(ValueError):
        cpt.RadiationProblem(body=sphere, wavenumber=1.0)
    with pytest.raises(ValueError):
        cpt.DiffractionProblem(body=sphere, wavenumber=1.0)


def test_mesh_above_the_free_surface(caplog):
    with caplog.at_level(logging.WARNING):
        body = cpt.FloatingBody(mesh=cpt.mesh_sphere(), dofs=cpt.rigid_body_dofs())
        pb = cpt.RadiationProblem(body=body, omega=1.0)
    assert '50 panels above the free surface' in caplog.text
    assert np.all(pb.body.mesh.vertices[:, 2].max() <= 0.0)


def test_mesh_below_the_sea_bottom(caplog):
    with caplog.at_level(logging.WARNING):
        body = cpt.FloatingBody(
                mesh=cpt.mesh_sphere(radius=0.5, center=(0, 0, -1.0)),
                dofs=cpt.rigid_body_dofs()
                )
        pb = cpt.RadiationProblem(body=body, omega=1.0, water_depth=1.0)
    assert '50 panels below the sea bottom' in caplog.text
    assert np.all(pb.body.mesh.vertices[:, 2].min() >= -1.0)


def test_backward_compatibility_with_sea_bottom_argument(caplog):
    with caplog.at_level(logging.WARNING):
        pb = cpt.DiffractionProblem(sea_bottom=-10.0)
        assert pb.water_depth == 10.0
    assert 'water_depth' in caplog.text


@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_setting_wavelength(water_depth):
    λ = 10*np.random.rand()
    assert np.isclose(cpt.DiffractionProblem(wavelength=λ, water_depth=water_depth).wavelength, λ)


@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_setting_wavenumber(water_depth):
    k = 10*np.random.rand()
    assert np.isclose(cpt.DiffractionProblem(wavenumber=k, water_depth=water_depth).wavenumber, k)


@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_setting_period(water_depth):
    T = 10*np.random.rand()
    assert np.isclose(cpt.DiffractionProblem(period=T, water_depth=water_depth).period, T)


def test_setting_too_many_frequencies():
    with pytest.raises(ValueError, match="at most one"):
        cpt.DiffractionProblem(omega=1.0, wavelength=1.0)


def test_diffraction_problem(sphere):
    assert cpt.DiffractionProblem(omega=1.0).body is None

    pb = cpt.DiffractionProblem(body=sphere, wave_direction=1.0, wavenumber=1.0)
    assert len(pb.boundary_condition) == sphere.mesh.nb_faces

    with pytest.raises(TypeError):
        cpt.DiffractionProblem(boundary_conditions=[0, 0, 0])

    assert "DiffractionProblem" in str(cpt.DiffractionProblem(wavenumber=1.0, g=10, rho=1025, free_surface=np.inf))

    res = pb.make_results_container()
    assert isinstance(res, DiffractionResult)


def test_wave_direction_radians_warning(sphere, caplog):
    with caplog.at_level(logging.WARNING):
        cpt.DiffractionProblem(body=sphere, omega=1.0, wave_direction=180)
    assert 'in radians and not in degrees' in caplog.text


def test_radiation_problem(sphere, caplog):
    pb = cpt.RadiationProblem(body=sphere, omega=1.0)
    assert len(pb.boundary_condition) == sphere.mesh.nb_faces

    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")
    pb2 = cpt.RadiationProblem(body=sphere, radiating_dof="Heave")
    assert np.all(pb.boundary_condition == pb2.boundary_condition)

    assert "RadiationProblem" in str(cpt.RadiationProblem(g=10, rho=1025, free_surface=np.inf))

    res = pb.make_results_container()
    assert isinstance(res, RadiationResult)
    assert res.added_mass == res.added_masses == {}
    assert res.radiation_damping == res.radiation_dampings == {}

    res = pb.make_results_container(forces={"Heave": 1.0 + 2.0j})
    assert res.added_mass == res.added_masses == {"Heave": 1.0}
    assert res.radiation_damping == res.radiation_dampings == {"Heave": 2.0}


def test_radiation_problem_with_wrong_dof_name():
    mesh = cpt.mesh_sphere()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    with pytest.raises(ValueError, match="'Flurp'"):
        cpt.RadiationProblem(body=body, wavelength=1.0, radiating_dof="Flurp")


@pytest.mark.parametrize("cal_file", ["Nemoh.cal", "Nemoh_v3.cal"])
def test_import_cal_file(cal_file):
    """Test the importation of legacy Nemoh.cal files."""
    current_file_path = os.path.dirname(os.path.abspath(__file__))

    # Non symmetrical body
    cal_file_path = os.path.join(current_file_path, "Nemoh_verification_cases", "NonSymmetrical", cal_file)
    problems = import_cal_file(cal_file_path)

    assert len(problems) == 6*41+41
    for problem in problems:
        assert problem.rho == 1000.0
        assert problem.g == 9.81
        assert problem.water_depth == np.inf
        assert isinstance(problem.body, cpt.FloatingBody)
        assert problem.body.nb_dofs == 6
        assert problem.body.mesh.nb_vertices == 299  # Duplicate vertices are removed during import.
        assert problem.body.mesh.nb_faces == 280
        assert problem.omega in np.linspace(0.1, 2.0, 41)
        if isinstance(problem, cpt.DiffractionProblem):
            assert problem.wave_direction == 0.0

    # Symmetrical cylinder
    cal_file_path = os.path.join(current_file_path, "Nemoh_verification_cases", "Cylinder", cal_file)
    problems = import_cal_file(cal_file_path)

    assert len(problems) == 6*2+2
    for problem in problems:
        assert problem.rho == 1000.0
        assert problem.g == 9.81
        assert problem.water_depth == np.inf
        assert isinstance(problem.body.mesh, cpt.ReflectionSymmetricMesh)
        assert isinstance(problem.body.mesh[0], cpt.Mesh)
        assert problem.body.nb_dofs == 6
        # assert problem.body.mesh.nb_vertices == 2*540
        assert problem.body.mesh.nb_faces == 2*300
        assert problem.omega == 0.1 or problem.omega == 2.0
        if isinstance(problem, cpt.DiffractionProblem):
            assert problem.wave_direction == 0.0


def test_results():
    assert isinstance(LinearPotentialFlowProblem().make_results_container(), LinearPotentialFlowResult)

    pb = cpt.DiffractionProblem(g=10, rho=1023, free_surface=np.inf)
    res = DiffractionResult(pb)
    assert res.g == pb.g == 10
    assert "DiffractionResult" in str(res)

    pb = cpt.RadiationProblem(g=10, rho=1023, free_surface=np.inf)
    res = RadiationResult(pb)
    assert res.g == pb.g == 10
    assert "RadiationResult" in str(res)


@pytest.fixture
def broken_bem_solver():
    ref_gf = cpt.Delhommeau()
    class BrokenGreenFunction:
        def evaluate(self, m1, m2, fs, wd, wavenumber, *args, **kwargs):
            if wavenumber < 2.0:
                raise NotImplementedError("I'm potato")
            else:
                return ref_gf.evaluate(m1, m2, fs, wd, wavenumber, *args, **kwargs)
    broken_bem_solver = cpt.BEMSolver(green_function=BrokenGreenFunction())
    return broken_bem_solver


def test_failed_resolution_failing(broken_bem_solver, sphere):
    pb = cpt.DiffractionProblem(body=sphere, wavenumber=1.0, wave_direction=0.0)
    with pytest.raises(NotImplementedError):
        broken_bem_solver.solve(pb)


def test_failed_resolution_catched(broken_bem_solver, sphere):
    from capytaine.bem.problems_and_results import FailedDiffractionResult
    pb = cpt.DiffractionProblem(body=sphere, wavenumber=1.0, wave_direction=0.0)
    failed_res = broken_bem_solver.solve_all([pb])[0]
    assert isinstance(failed_res, FailedDiffractionResult)
