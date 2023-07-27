import pytest
import numpy as np
import capytaine as cpt

@pytest.fixture
def body():
    mesh = cpt.mesh_sphere().immersed_part()
    return cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())

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
    pb = cpt.RadiationProblem(body=body, omega=2.0, forward_speed=1.0, wave_direction=0.0, radiating_dof="Heave")
    assert pb.forward_speed == 1.0
    assert pb.encounter_omega < pb.omega

def test_solve(body, solver):
    pb = cpt.DiffractionProblem(body=body, omega=2.0, forward_speed=1.0, wave_direction=0.0)
    res = solver.solve(pb)
    print(res.forces)

