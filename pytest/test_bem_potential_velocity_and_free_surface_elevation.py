"""Tests related to post-processing computation of potential, velocity and free surface elevation."""

import pytest
from pytest import approx
import numpy as np
import capytaine as cpt

#######################################################################
#                  Test shapes of inputs and outputs                  #
#######################################################################

@pytest.fixture
def solver():
    return cpt.BEMSolver()

@pytest.fixture
def result(solver):
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh)
    body.add_translation_dof(name="Heave")
    result = solver.solve(cpt.DiffractionProblem(body=body, omega=1.0, wave_direction=np.pi/4), keep_details=True)
    return result


##################
#  Single point  #
##################

def test_airy_waves_potential_at_point(result):
    from capytaine.bem.airy_waves import airy_waves_potential
    point = (0.0, 0.0, -3.0)
    phi = airy_waves_potential(point, result)
    assert phi.shape == (1,)

def test_compute_potential_at_point(solver, result):
    point = (0.0, 0.0, -3.0)
    phi = solver.compute_potential(point, result)
    assert phi.shape == (1,)

def test_airy_waves_velocity_at_point(result):
    from capytaine.bem.airy_waves import airy_waves_velocity
    point = (0.0, 0.0, -3.0)
    u = airy_waves_velocity(point, result)
    assert u.shape == (1, 3)

def test_compute_velocity_at_point(solver, result):
    point = (0.0, 0.0, -3.0)
    u = solver.compute_velocity(point, result)
    assert u.shape == (1, 3)

def test_airy_waves_free_surface_elevation_at_point(result):
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    point = (0.0, 3.0)
    fse = airy_waves_free_surface_elevation(point, result)
    assert fse.shape == (1,)
    point_ = (0.0, 3.0, 0.0)
    fse_ = airy_waves_free_surface_elevation(point_, result)
    assert fse_.shape == (1,)
    assert np.allclose(fse, fse_)

def test_compute_free_surface_elevation_at_point(solver, result):
    point = (0.0, 3.0)
    fse = solver.compute_free_surface_elevation(point, result)
    assert fse.shape == (1,)
    point_ = (0.0, 3.0, 0.0)
    fse_ = solver.compute_free_surface_elevation(point_, result)
    assert fse_.shape == (1,)
    assert np.allclose(fse, fse_)


#################
#  Points list  #
#################

def test_airy_waves_potential_at_points_list(result):
    from capytaine.bem.airy_waves import airy_waves_potential
    points = [(0.0, 0.0, -3.0), (0.0, 1.0, -2.0), (1.0, 1.0, -1.0)]
    phi = airy_waves_potential(points, result)
    assert phi.shape == (3,)

def test_compute_potential_at_points_list(solver, result):
    points = [(0.0, 0.0, -3.0), (0.0, 1.0, -2.0), (1.0, 1.0, -1.0)]
    phi = solver.compute_potential(points, result)
    assert phi.shape == (3,)

def test_airy_waves_velocity_at_points_list(result):
    from capytaine.bem.airy_waves import airy_waves_velocity
    points = [(0.0, 0.0, -3.0), (0.0, 1.0, -2.0), (1.0, 1.0, -1.0)]
    u = airy_waves_velocity(points, result)
    assert u.shape == (3, 3)

def test_compute_velocity_at_points_list(solver, result):
    points = [(0.0, 0.0, -3.0), (0.0, 1.0, -2.0), (1.0, 1.0, -1.0)]
    u = solver.compute_velocity(points, result)
    assert u.shape == (3, 3)

def test_airy_waves_free_surface_elevation_at_points_list(result):
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    points = [(0.0, 3.0), (1.0, -2.0), (1.0, -1.0)]
    fse = airy_waves_free_surface_elevation(points, result)
    assert fse.shape == (3,)
    points_ = [(0.0, 3.0, 0.0), (1.0, -2.0, 0.0), (1.0, -1.0, 0.0)]
    fse_ = airy_waves_free_surface_elevation(points_, result)
    assert fse_.shape == (3,)
    assert np.allclose(fse, fse_)

def test_compute_free_surface_elevation_at_points_list(solver, result):
    points = [(0.0, 3.0), (1.0, -2.0), (1.0, -1.0)]
    fse = solver.compute_free_surface_elevation(points, result)
    assert fse.shape == (3,)
    points_ = [(0.0, 3.0, 0.0), (1.0, -2.0, 0.0), (1.0, -1.0, 0.0)]
    fse_ = solver.compute_free_surface_elevation(points_, result)
    assert fse_.shape == (3,)
    assert np.allclose(fse, fse_)


##############
#  Meshgrid  #
##############

def test_airy_waves_potential_at_meshgrid(result):
    from capytaine.bem.airy_waves import airy_waves_potential
    points = np.meshgrid(np.linspace(0.0, 1.0, 2), np.linspace(-1.0, 1.0, 3), np.linspace(-2.0, 0.0, 4))
    phi = airy_waves_potential(points, result)
    assert phi.shape == points[0].shape

def test_compute_potential_at_meshgrid(solver, result):
    points = np.meshgrid(np.linspace(2.0, 3.0, 2), np.linspace(-1.0, 1.0, 3), np.linspace(-2.0, -1.0, 4))
    phi = solver.compute_potential(points, result)
    assert phi.shape == points[0].shape

def test_airy_waves_velocity_at_meshgrid(result):
    from capytaine.bem.airy_waves import airy_waves_velocity
    points = np.meshgrid(np.linspace(0.0, 1.0, 2), np.linspace(-1.0, 1.0, 3), np.linspace(-2.0, 0.0, 4))
    u = airy_waves_velocity(points, result)
    assert u.shape == (*points[0].shape, 3)

def test_compute_velocity_at_meshgrid(solver, result):
    points = np.meshgrid(np.linspace(2.0, 3.0, 2), np.linspace(4.0, 5.0, 3), np.linspace(-2.0, -1.0, 4))
    u = solver.compute_velocity(points, result)
    assert u.shape == (*points[0].shape, 3)

def test_airy_waves_free_surface_elevation_at_meshgrid(result):
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    points = np.meshgrid(np.linspace(0.0, 1.0, 2), np.linspace(-1.0, 1.0, 3))
    fse = airy_waves_free_surface_elevation(points, result)
    assert fse.shape == points[0].shape
    points_ = np.meshgrid(np.linspace(0.0, 1.0, 2), np.linspace(-1.0, 1.0, 3), np.linspace(0.0, 0.0, 1))
    fse_ = airy_waves_free_surface_elevation(points_, result)
    assert fse_.shape == points_[0].shape
    assert np.allclose(fse, fse_.squeeze())

def test_compute_free_surface_elevation_at_meshgrid(solver, result):
    points = np.meshgrid(np.linspace(2.0, 3.0, 2), np.linspace(3.0, 4.0, 3))
    fse = solver.compute_free_surface_elevation(points, result)
    assert fse.shape == points[0].shape
    points_ = np.meshgrid(np.linspace(2.0, 3.0, 2), np.linspace(3.0, 4.0, 3), np.linspace(0.0, 0.0, 1))
    fse_ = solver.compute_free_surface_elevation(points_, result)
    assert fse_.shape == points_[0].shape
    assert np.allclose(fse, fse_.squeeze())


##########
#  Mesh  #
##########

def test_airy_waves_potential_on_mesh(result):
    from capytaine.bem.airy_waves import airy_waves_potential
    mesh = cpt.mesh_vertical_cylinder(radius=2.0, length=5.0).immersed_part()
    phi = airy_waves_potential(mesh, result)
    assert phi.shape == (mesh.nb_faces,)

def test_compute_potential_on_mesh(solver, result):
    mesh = cpt.mesh_vertical_cylinder(radius=2.0, length=5.0).immersed_part()
    phi = solver.compute_potential(mesh, result)
    assert phi.shape == (mesh.nb_faces,)

def test_airy_waves_velocity_on_mesh(result):
    from capytaine.bem.airy_waves import airy_waves_velocity
    mesh = cpt.mesh_vertical_cylinder(radius=2.0, length=5.0).immersed_part()
    u = airy_waves_velocity(mesh, result)
    assert u.shape == (mesh.nb_faces, 3)

def test_compute_velocity_on_mesh(solver, result):
    mesh = cpt.mesh_vertical_cylinder(radius=2.0, length=5.0).immersed_part()
    u = solver.compute_velocity(mesh, result)
    assert u.shape == (mesh.nb_faces, 3)

def test_airy_waves_free_surface_elevation_on_mesh(result):
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    mesh = cpt.mesh_rectangle(center=(0, 0, 0), normal=(0, 0, -1), resolution=(3, 3), size=(2.0, 2.0))
    fse = airy_waves_free_surface_elevation(mesh, result)
    assert fse.shape == (mesh.nb_faces,)

def test_compute_free_surface_elevation_on_mesh(solver, result):
    mesh = cpt.mesh_rectangle(center=(0, 0, 0), normal=(0, 0, -1), resolution=(3, 3), size=(2.0, 2.0))
    fse = solver.compute_free_surface_elevation(mesh, result)
    assert fse.shape == (mesh.nb_faces,)


#################
#  FreeSurface  #
#################

def test_airy_waves_free_surface_elevation_on_free_surface(result):
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    fs = cpt.FreeSurface(nx=3, ny=3)
    fse = airy_waves_free_surface_elevation(fs, result)
    assert fse.shape == (fs.mesh.nb_faces,)

def test_compute_free_surface_elevation_on_free_surface(solver, result):
    fs = cpt.FreeSurface(nx=3, ny=3)
    fse = solver.compute_free_surface_elevation(fs, result)
    assert fse.shape == (fs.mesh.nb_faces,)


#######################################################################
#                            Check values                             #
#######################################################################

def test_pressure_integration(solver, result):
    f = result.body.integrate_pressure(solver.compute_pressure(result.body.mesh, result))
    assert f == pytest.approx(result.forces)

def test_reconstruction_of_given_boundary_condition(solver, result):
    velocities = solver.compute_velocity(result.body.mesh, result)
    normal_velocities = np.einsum('...k,...k->...', velocities, result.body.mesh.faces_normals)
    np.testing.assert_allclose(normal_velocities, result.problem.boundary_condition)

def test_airy_wave_free_surface_elevation_values():
    from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
    pb = cpt.DiffractionProblem(wave_direction=0.0, wavelength=1.0)

    assert np.isclose(np.real(airy_waves_free_surface_elevation([0, 0], pb)), 1.0)
    assert np.isclose(np.real(airy_waves_free_surface_elevation([0.25, 0], pb)), 0.0, atol=1e-5)
    assert np.isclose(np.real(airy_waves_free_surface_elevation([0.5, 0], pb)), -1.0)

def test_integrated_pressure(solver, result):
    pressure = solver.compute_pressure(result.body.mesh, result)
    forces = result.body.integrate_pressure(pressure)
    assert result.forces == approx(forces)

def test_fse_zero_frequency(solver):
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh)
    body.add_translation_dof(name="Heave")
    pb = cpt.RadiationProblem(body=body, omega=0.0, radiating_dof="Heave")
    res = solver.solve(pb, keep_details=True)
    points = -np.random.rand(10, 3)
    pot = solver.compute_potential(points, res)
    fse = solver.compute_free_surface_elevation(points, res)
    with pytest.raises(TypeError):
        vel = solver.compute_velocity(points, res)

def test_direct_solver(solver):
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh)
    body.add_translation_dof(name="Heave")
    res = solver.solve(cpt.DiffractionProblem(body=body, omega=1.0, wave_direction=np.pi/4), keep_details=True, method='direct')
    points = [(0.0, 0.0, -3.0), (0.0, 1.0, -2.0), (1.0, 1.0, -1.0)]
    with pytest.raises(Exception, match="direct method"):
        pot = solver.compute_potential(points, res)
