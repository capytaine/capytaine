import pytest

import numpy as np
import capytaine as cpt


def test_froude_krylov_force():
    # Froude-Krylov force should be unchanged by the lid
    from capytaine.bem.airy_waves import froude_krylov_force
    mesh = cpt.mesh_sphere().immersed_part()
    body_without_lid = cpt.FloatingBody(
            mesh=mesh,
            dofs=cpt.rigid_body_dofs()
            )
    pb_without_lid = cpt.DiffractionProblem(body=body_without_lid)

    lid_mesh = mesh.generate_lid(z=-0.1, info=False)
    body_with_lid = cpt.FloatingBody(
            mesh=mesh,
            lid_mesh=lid_mesh,
            dofs=cpt.rigid_body_dofs()
            )
    pb_with_lid = cpt.DiffractionProblem(body=body_with_lid)

    fk_with_lid = froude_krylov_force(pb_with_lid)["Heave"]
    fk_without_lid = froude_krylov_force(pb_without_lid)["Heave"]
    assert fk_with_lid == pytest.approx(fk_without_lid)


def test_lid_below_free_surface():
    mesh = cpt.AxialSymmetricMesh.from_profile(lambda z: (z + 1.0)**2, np.linspace(-1.0, 0.0, 10)).merged()
    lid_mesh = mesh.generate_lid(z=-0.5)
    x, y, z = lid_mesh.faces_centers.T
    assert np.all(np.hypot(x, y) <= (z + 1.0)**2)
