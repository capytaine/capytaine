import numpy as np
import capytaine as cpt
import pytest

method = ['indirect', 'direct']

@pytest.mark.parametrize("method", method)
def test_sum_of_dofs(method):
    body1 = cpt.FloatingBody(mesh=cpt.mesh_sphere(radius=1.0, resolution=(3, 12), center=(0, 0, -3)), name="body1")
    body1.add_translation_dof(name="Heave")

    body2 = cpt.FloatingBody(mesh=cpt.mesh_sphere(radius=1.0, resolution=(3, 8), center=(5, 3, -1.5)), name="body2")
    body2.add_translation_dof(name="Heave")

    both = body1 + body2
    both.add_translation_dof(name="Heave")

    problems = [cpt.RadiationProblem(body=both, radiating_dof=dof, omega=1.0) for dof in both.dofs]
    solver = cpt.BEMSolver()
    results = solver.solve_all(problems, method=method)
    dataset = cpt.assemble_dataset(results)

    both_added_mass = dataset['added_mass'].sel(radiating_dof="Heave", influenced_dof="Heave").data
    body1_added_mass = dataset['added_mass'].sel(radiating_dof="body1__Heave", influenced_dof="body1__Heave").data
    body2_added_mass = dataset['added_mass'].sel(radiating_dof="body2__Heave", influenced_dof="body2__Heave").data

    assert np.allclose(both_added_mass, body1_added_mass + body2_added_mass, rtol=1e-2)


@pytest.mark.parametrize("method", method)
def test_rotation_axis(method):
    body = cpt.FloatingBody(mesh=cpt.mesh_parallelepiped(resolution=(4, 4, 4), center=(0, 0, -1), name="body"))
    body.add_translation_dof(name="Sway")
    body.add_rotation_dof(axis=cpt.Axis(point=(0, 0, 0), vector=(0, 0, 1)), name="Yaw")

    l = 2.0
    body.add_rotation_dof(axis=cpt.Axis(point=(l, 0, 0), vector=(0, 0, 1)), name="other_rotation")

    assert np.allclose(body.dofs['other_rotation'], (body.dofs['Yaw'] - l*body.dofs['Sway']))

    problems = [cpt.RadiationProblem(body=body, radiating_dof=dof, omega=1.0) for dof in body.dofs]
    solver = cpt.BEMSolver()
    results = solver.solve_all(problems, keep_details=True, method=method)
    dataset = cpt.assemble_dataset(results)

    if ( method == 'indirect' ):
      sources = {result.radiating_dof: result.sources for result in results}
      assert np.allclose(sources['other_rotation'],
                       sources['Yaw'] - l*sources['Sway'], atol=1e-4)

    potential = {result.radiating_dof: result.potential for result in results}
    assert np.allclose(potential['other_rotation'],
                       potential['Yaw'] - l*potential['Sway'], atol=1e-4)

    A_m = dataset['added_mass'].sel(radiating_dof="other_rotation", influenced_dof="other_rotation").data
    A = dataset['added_mass'].sel(radiating_dof=["Yaw", "Sway"], influenced_dof=["Yaw", "Sway"]).data
    P = np.array([1, -l])
    assert np.isclose(A_m, P.T @ A @ P)
