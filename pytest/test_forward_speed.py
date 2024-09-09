import pytest
from pytest import approx
import numpy as np
import pandas as pd
import capytaine as cpt
import xarray as xr

# PROBLEM DEFINITION

@pytest.fixture
def body():
    mesh = cpt.mesh_vertical_cylinder(radius=1.0, length=1.01, center=(0.0, 0.0, -0.5), resolution=(2, 20, 10)).immersed_part(water_depth=1.0)
    body = cpt.FloatingBody(mesh=mesh, name="body")
    body.add_translation_dof(name="Surge")
    return body

@pytest.fixture
def solver():
    return cpt.BEMSolver()

def test_encouter_frequency_along_waves(body):
    pb = cpt.DiffractionProblem(body=body, omega=2.0, forward_speed=1.0, wave_direction=0.0)
    assert pb.forward_speed == 1.0
    assert pb.encounter_omega < pb.omega  # Object is moving in the same direction as the waves
    assert pb.encounter_wave_direction == pb.wave_direction

def test_encouter_frequency_against_waves(body):
    pb = cpt.DiffractionProblem(body=body, omega=2.0, forward_speed=-1.0, wave_direction=0.0)
    assert pb.forward_speed == -1.0
    assert pb.encounter_omega > pb.omega  # Object is moving against the wave
    assert pb.encounter_wave_direction == pb.wave_direction

def test_encouter_frequency_faster_than_waves(body):
    pb = cpt.DiffractionProblem(body=body, omega=2.0, forward_speed=10.0, wave_direction=0.0)
    assert pb.forward_speed == 10.0
    assert pb.encounter_omega >= 0.0
    assert pb.encounter_wave_direction == approx(pb.wave_direction + np.pi)

def test_encouter_frequency_orthogonal_waves(body):
    pb = cpt.DiffractionProblem(body=body, omega=2.0, forward_speed=1.0, wave_direction=np.pi/2)
    assert pb.forward_speed == 1.0
    assert pb.encounter_omega == approx(pb.omega)

def test_encounter_frequency_radiation_problem(body):
    pb = cpt.RadiationProblem(body=body, omega=2.0, forward_speed=1.0, wave_direction=0.0, radiating_dof="Surge")
    assert pb.forward_speed == 1.0
    assert pb.encounter_omega < pb.omega

def test_multibody(body, solver):
    two_bodies = body + body.translated_x(5.0, name="other_body")
    with pytest.raises(NotImplementedError):
        pb = cpt.RadiationProblem(body=two_bodies, omega=2.0, forward_speed=1.0, radiating_dof="body__Surge")

def test_non_rigid_body(body, solver):
    body = body.copy(name="body")
    body.dofs["Shear"] = np.array([[z, 0, 0] for (x, y, z) in body.mesh.faces_centers])
    with pytest.raises(NotImplementedError):
        pb = cpt.RadiationProblem(body=body, omega=2.0, forward_speed=1.0, radiating_dof="Shear")

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

def test_fill_dataset_with_forward_speed(body, solver):
    test_matrix = xr.Dataset(coords={
        "omega": [1.0],
        "forward_speed": [0.0, 1.0],
        "wave_direction": [0.0, np.pi/2],
        "radiating_dof": ["Surge"],
        })
    ds = solver.fill_dataset(test_matrix, body)
    assert "wave_direction" in ds.added_mass.coords
    assert "encounter_omega" in ds.coords
    assert "encounter_wave_direction" in ds.coords

def test_fill_dataset_without_forward_speed(body, solver):
    test_matrix = xr.Dataset(coords={
        "omega": [1.0],
        "wave_direction": [0.0, np.pi/2],
        "radiating_dof": ["Surge"],
        })
    ds = solver.fill_dataset(test_matrix, body)
    assert "wave_direction" not in ds.added_mass.coords
    assert "encounter_omega" not in ds.coords
    assert "encounter_wave_direction" not in ds.coords

def test_no_explicit_wave_direction(body):
    test_matrix = xr.Dataset(coords={
        "wavelength": np.linspace(0.1, 10.0, 2),
        "forward_speed": [0.0, 10.0],
        "radiating_dof": ["Surge"],
        })
    ds = cpt.BEMSolver().fill_dataset(test_matrix, body)
    assert ds.wave_direction.values == np.array([0.0])


# VALIDATION

rho = 1025
g = 9.81

# Based on Figure 1 from Malenica 1995
MALENICA_EXCITATION_FORCE = pd.DataFrame([
    (0.269, 1.662, 0.05), (0.615, 3.192, 0.05), (1.086, 3.168, 0.05), (1.568, 2.225, 0.05), (1.955, 1.608, 0.05),
    (0.201, 1.281, 0.0), (0.506, 2.928, 0.0), (1.308, 2.733, 0.0), (1.840, 1.890, 0.0),
    (0.353, 2.263, -0.05), (0.568, 3.224, -0.05), (1.171, 2.963, -0.05), (1.977, 1.893, -0.05),
    ], columns=["wavenumber", "normalized_force", "froude_number"])

@pytest.mark.parametrize("ref_data", MALENICA_EXCITATION_FORCE.iterrows())
def test_malenica_excitation_force(body, solver, ref_data):
    from capytaine.bem.airy_waves import froude_krylov_force
    pb = cpt.DiffractionProblem(
        body=body,
        wavenumber=ref_data[1].wavenumber,
        forward_speed=ref_data[1].froude_number*np.sqrt(g),
        wave_direction=0.0,
        water_depth=1.0,
        rho=rho
    )
    res = solver.solve(pb)
    excitation_force = res.forces["Surge"] + froude_krylov_force(pb)["Surge"]
    assert np.abs(excitation_force)/(g*rho) == approx(ref_data[1].normalized_force, rel=5e-2)


MALENICA_ADDED_MASS = pd.DataFrame([
    (0.620, 3.302, 0.05), (0.855, 2.634, 0.05), (1.303, 1.407, 0.05),
    (0.510, 3.433, 0.0), (0.916, 2.217, 0.0), (1.495, 0.916, 0.0),
    (0.602, 3.172, -0.05), (1.069, 1.518, -0.05),
    ], columns=["wavenumber", "normalized_added_mass", "froude_number"])

@pytest.mark.parametrize("ref_data", MALENICA_ADDED_MASS.iterrows())
def test_malenica_added_mass(body, solver, ref_data):
    from capytaine.bem.airy_waves import froude_krylov_force
    pb = cpt.RadiationProblem(
        body=body,
        wavenumber=ref_data[1].wavenumber,
        forward_speed=ref_data[1].froude_number*np.sqrt(g),
        wave_direction=0.0,
        water_depth=1.0,
        radiating_dof="Surge",
        rho=rho
    )
    res = solver.solve(pb)
    assert np.abs(res.added_masses["Surge"])/(rho) == approx(ref_data[1].normalized_added_mass, rel=1e-1)

MALENICA_RADIATION_DAMPING = pd.DataFrame([
    (0.456, 0.965, 0.05), (0.727, 1.917, 0.05), (1.129, 2.284, 0.05),
    (0.434, 0.978, 0.0), (0.674, 1.884, 0.0), (1.125, 2.218, 0.0),
    (0.408, 0.971, -0.05), (0.642, 1.906, -0.05), (1.042, 2.234, -0.05),
    ], columns=["wavenumber", "normalized_radiation_damping", "froude_number"])

@pytest.mark.parametrize("ref_data", MALENICA_RADIATION_DAMPING.iterrows())
def test_malenica_radiation_damping(body, solver, ref_data):
    from capytaine.bem.airy_waves import froude_krylov_force
    pb = cpt.RadiationProblem(
        body=body,
        wavenumber=ref_data[1].wavenumber,
        forward_speed=ref_data[1].froude_number*np.sqrt(g),
        wave_direction=0.0,
        water_depth=1.0,
        radiating_dof="Surge",
        rho=rho
    )
    res = solver.solve(pb)
    assert np.abs(res.radiation_dampings["Surge"])/(rho*pb.encounter_omega) == approx(ref_data[1].normalized_radiation_damping, rel=1e-1)


def test_near_zero_encounter_frequency_radiation(body, solver):
    pb = cpt.RadiationProblem(body=body, omega=1.0, forward_speed=9.805, radiating_dof="Surge")
    assert pb.encounter_omega > 0.0
    solver.solve(pb)


def test_near_zero_encounter_frequency_diffraction(body, solver):
    pb = cpt.DiffractionProblem(body=body, omega=1.0, forward_speed=9.805)
    assert pb.encounter_omega > 0.0
    solver.solve(pb)


def test_zero_encounter_frequency_radiation(body, solver):
    pb = cpt.RadiationProblem(body=body, omega=1.0, forward_speed=9.81, radiating_dof="Surge")
    assert float(pb.encounter_omega) == 0.0
    res = solver.solve(pb)
    assert not np.isinf(res.forces["Surge"])
    assert not np.isnan(res.forces["Surge"])
    assert np.isnan(res.added_mass["Surge"])
    assert np.isnan(res.radiation_dampings["Surge"])


def test_zero_encounter_frequency_diffraction(body, solver):
    pb = cpt.DiffractionProblem(body=body, omega=1.0, forward_speed=9.81)
    assert float(pb.encounter_omega) == 0.0
    res = solver.solve(pb)
    assert not np.isinf(res.forces["Surge"])
    assert not np.isnan(res.forces["Surge"])
