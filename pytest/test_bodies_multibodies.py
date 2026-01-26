import pytest

import numpy as np
import capytaine as cpt
import xarray as xr

from capytaine.bodies.multibodies import Multibody


def test_multibody_without_names():
    body_1 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(0, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            name="body"
            )
    body_2 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(2, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(2, 0, 0)),
            name="body"
            )
    with pytest.raises(ValueError):
        Multibody([body_1, body_2],)


def test_multibody_resolution():
    body_1 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(0, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            name="body_1",
            )
    body_2 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(2, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(2, 0, 0)),
            name="body_2",
            )
    multi = Multibody(
        [body_1, body_2],
        # own_dofs={
        #     "Heave": np.array([[0.0, 0.0, 1.0] for _ in range(body_1.mesh.nb_faces + body_2.mesh.nb_faces)])
        # }
    )
    solver = cpt.BEMSolver()
    res = solver.solve(cpt.DiffractionProblem(body=multi, omega=1.0, wave_direction=0.0))
    ref_res = solver.solve(cpt.DiffractionProblem(body=multi.as_FloatingBody(), omega=1.0, wave_direction=0.0))
    assert all(np.isclose(ref_res.forces[k], res.forces[k]) for k in multi.dofs)

    res = solver.solve(cpt.RadiationProblem(body=multi, omega=1.0, radiating_dof="body_1__Heave"))
    ref_res = solver.solve(cpt.RadiationProblem(body=multi.as_FloatingBody(), omega=1.0, radiating_dof="body_1__Heave"))
    assert all(np.isclose(ref_res.forces[k], res.forces[k]) for k in multi.dofs)

def test_multibody_hydrostatics():
    body_1 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(0, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            center_of_mass=(0, 0, 0),
            mass=1000.0,
            name="body_1",
            )
    body_2 = cpt.FloatingBody(
            mesh=cpt.mesh_sphere(center=(2, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(2, 0, 0)),
            center_of_mass=(2, 0, 0),
            mass=1000.0,
            name="body_2",
            )
    multi = Multibody(
            [body_1, body_2],
            # own_dofs={
            #     "Heave": np.array([[0.0, 0.0, 1.0] for _ in range(body_1.mesh.nb_faces + body_2.mesh.nb_faces)])
            #     }
            )
    assert np.allclose(multi.center_of_mass['body_1'], body_1.center_of_mass)
    assert np.allclose(multi.center_of_mass['body_2'], body_2.center_of_mass)
    assert np.allclose(multi.center_of_buoyancy['body_1'], body_1.center_of_buoyancy)
    assert np.allclose(multi.center_of_buoyancy['body_2'], body_2.center_of_buoyancy)

def test_hydrostatic_stiffness():
    body_1 = cpt.FloatingBody(
            mesh=cpt.mesh_parallelepiped(center=(0, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            center_of_mass=(0, 0, 0),
            mass=1000.0,
            name="body_1",
            )
    body_2 = cpt.FloatingBody(
            mesh=cpt.mesh_parallelepiped(center=(2, 0, 0)).immersed_part(),
            dofs=cpt.rigid_body_dofs(rotation_center=(2, 0, 0)),
            center_of_mass=(2, 0, 0),
            mass=1000.0,
            name="body_2",
            )
    multi = Multibody(
            [body_1, body_2],
            # own_dofs={
            #     "Heave": np.array([[0.0, 0.0, 1.0] for _ in range(body_1.mesh.nb_faces + body_2.mesh.nb_faces)])
            #     }
            )
    h = multi.compute_hydrostatic_stiffness()
    assert np.allclose(h.values[:6, :6], body_1.compute_hydrostatic_stiffness().values)
    assert np.allclose(h.values[6:, 6:], body_2.compute_hydrostatic_stiffness().values)
    assert np.allclose(h.values[6:, :6], 0.0)
    assert np.allclose(h.values[:6, 6:], 0.0)

    m = multi.compute_rigid_body_inertia()
    assert np.allclose(m.values[:6, :6], body_1.compute_rigid_body_inertia().values)
    assert np.allclose(m.values[6:, 6:], body_2.compute_rigid_body_inertia().values)
    assert np.allclose(m.values[6:, :6], 0.0)
    assert np.allclose(m.values[:6, 6:], 0.0)
