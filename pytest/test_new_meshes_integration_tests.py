import numpy as np

import capytaine as cpt
import xarray as xr

import capytaine.new_meshes as meshes
from capytaine.meshes.mesh_like_protocol import MeshLike

def test_meshes_are_meshlike():
    mesh = meshes.Mesh(
        vertices=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, -1.0]]),
        faces=np.array([[0, 1, 2, 2]])
    )
    assert isinstance(mesh, MeshLike)

def test_meshes_can_be_used_to_solve_pb():
    mesh = meshes.Mesh(
        vertices=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, -1.0]]),
        faces=np.array([[0, 1, 2, 2]])
    )
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb = cpt.RadiationProblem(body=body, omega=1.0, radiating_dof="Heave")
    solver.solve(pb, method="indirect")

def test_minimal_implementation_compute_velocity():
    mesh = meshes.Mesh(
        vertices=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, -1.0]]),
        faces=np.array([[0, 1, 2, 2]])
    )
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb = cpt.RadiationProblem(body=body, omega=1.0, radiating_dof="Heave")
    res = solver.solve(pb, method="indirect")
    solver.compute_velocity(mesh, res)

def test_minimal_implementation_solve_pb_with_lid():
    mesh = meshes.Mesh(
        vertices=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, -1.0]]),
        faces=np.array([[0, 1, 2, 2]])
    )
    lid_mesh = meshes.Mesh(
        vertices=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]),
        faces=np.array([[0, 1, 2, 2]])
    )
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb = cpt.RadiationProblem(body=body, omega=1.0, radiating_dof="Heave")
    solver.solve(pb, method="indirect")

def test_minimal_implementation_fill_dataset():
    mesh = meshes.Mesh(
        vertices=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, -1.0]]),
        faces=np.array([[0, 1, 2, 2]])
    )
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    test_matrix = xr.Dataset(
        coords={"wavenumber": 1.0, "radiating_dof": ["Heave"], "wave_direction": 1.0}
    )
    solver.fill_dataset(test_matrix, body, method="direct")
