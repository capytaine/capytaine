import pytest
from pytest import approx
import numpy as np
import capytaine as cpt
import xarray as xr

# PROBLEM DEFINITION

@pytest.fixture
def body():
    mesh = cpt.mesh_vertical_cylinder(radius=1.0, length=1.01, center=(0.0, 0.0, -0.5), resolution=(2, 20, 10)).immersed_part(water_depth=1.0)
    body = cpt.FloatingBody(mesh=mesh)
    body.add_translation_dof(name="Surge")
    return body

@pytest.fixture
def solver():
    return cpt.BEMSolver()

def test_encouter_frequency_along_waves(body):
    pb = cpt.DiffractionProblem(body=body, omega=2.0, forward_speed=1.0, wave_direction=0.0)
    assert pb.forward_speed == 1.0
    assert pb.encounter_omega < pb.omega  # Object is moving in the same direction as the waves

def test_encouter_frequency_against_waves(body):
    pb = cpt.DiffractionProblem(body=body, omega=2.0, forward_speed=-1.0, wave_direction=0.0)
    assert pb.forward_speed == -1.0
    assert pb.encounter_omega > pb.omega  # Object is moving against the wave

def test_encouter_frequency_orthogonal_waves(body):
    pb = cpt.DiffractionProblem(body=body, omega=2.0, forward_speed=1.0, wave_direction=np.pi/2)
    assert pb.forward_speed == 1.0
    assert pb.encounter_omega == pytest.approx(pb.omega)

def test_encounter_frequency_radiation_problem(body):
    pb = cpt.RadiationProblem(body=body, omega=2.0, forward_speed=1.0, wave_direction=0.0, radiating_dof="Surge")
    assert pb.forward_speed == 1.0
    assert pb.encounter_omega < pb.omega

def test_solve_diffraction_problem(body, solver):
    pb = cpt.DiffractionProblem(body=body, omega=2.0, forward_speed=1.0, wave_direction=0.0)
    res = solver.solve(pb)
    print(res.forces)

# POST-PROCESSING

@pytest.fixture
def result(body, solver):
    pb = cpt.DiffractionProblem(body=body, omega=2.0, forward_speed=1.0, wave_direction=0.0)
    return solver.solve(pb)

def test_pressure_reconstruction(result, solver):
    pressure = solver.compute_pressure(result.body.mesh, result)
    assert result.body.integrate_pressure(pressure) == approx(result.forces)

def test_velocity_reconstruction(result, solver):
    points = np.meshgrid(np.linspace(-10.0, 10.0, 3), np.linspace(-10.0, 10.0, 3), np.linspace(-10.0, -1.0, 2))
    velocity = solver.compute_velocity(points, result)

def test_free_surface_elevation(result, solver):
    points = np.meshgrid(np.linspace(-10.0, 10.0, 3), np.linspace(-10.0, 10.0, 3))
    fse = solver.compute_free_surface_elevation(points, result)

# DATASETS

def test_problem_from_dataset(body, solver):
    from capytaine.io.xarray import problems_from_dataset
    test_matrix = xr.Dataset(coords={
        "omega": [1.0],
        "forward_speed": [0.0, 1.0],
        "wave_direction": [0.0, np.pi/2],
        "radiating_dof": ["Surge"],
        })
    pbs = problems_from_dataset(test_matrix, body)
    assert len(pbs) == 7
    # Four diffraction problems + Three radiation problems, since the wave_direction is only relevant when forward_speed is non 0.0

def test_fill_dataset(body, solver):
    test_matrix = xr.Dataset(coords={
        "omega": [1.0],
        "forward_speed": [0.0, 1.0],
        "wave_direction": [0.0, np.pi/2],
        "radiating_dof": ["Surge"],
        })
    ds = solver.fill_dataset(test_matrix, body)
    assert "wave_direction" in ds.added_mass.coords

# VALIDATION

def test_malenica_single_test_case(body, solver):
    rho = 1025
    froude = 0.05
    pb = cpt.DiffractionProblem(body=body, wavenumber=2.0, forward_speed=froude*np.sqrt(9.81), wave_direction=0.0, rho=rho)
    res = solver.solve(pb)
    # assert res.forces["Surge"] == approx(0.66*rho**2)

