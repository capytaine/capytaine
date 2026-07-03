# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
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

    both = (body1 + body2).as_FloatingBody()
    both.add_translation_dof(name="Heave")

    problems = [cpt.RadiationProblem(body=both, radiating_dof=dof, omega=1.0) for dof in both.dofs]
    solver = cpt.BEMSolver(method=method)
    results = solver.solve_all(problems)
    dataset = cpt.assemble_dataset(results)

    both_added_mass = dataset['added_mass'].sel(radiating_dof="Heave", influenced_dof="Heave").data
    body1_added_mass = dataset['added_mass'].sel(radiating_dof="body1__Heave", influenced_dof="body1__Heave").data
    body2_added_mass = dataset['added_mass'].sel(radiating_dof="body2__Heave", influenced_dof="body2__Heave").data

    assert np.allclose(both_added_mass, body1_added_mass + body2_added_mass, rtol=1e-2)


@pytest.mark.parametrize("method", method)
def test_rotation_axis(method):
    body = cpt.FloatingBody(mesh=cpt.mesh_parallelepiped(resolution=(4, 4, 4), center=(0, 0, -1), name="body"))
    body.add_translation_dof(name="Sway")
    body.add_rotation_dof(rotation_center=(0, 0, 0), direction=(0, 0, 1), name="Yaw")

    l = 2.0
    body.add_rotation_dof(rotation_center=(l, 0, 0), direction=(0, 0, 1), name="other_rotation")

    assert np.allclose(
        body.dofs['other_rotation'].evaluate_motion(body.mesh),
        (body.dofs['Yaw'].evaluate_motion(body.mesh) - l*body.dofs['Sway'].evaluate_motion(body.mesh))
    )

    problems = [cpt.RadiationProblem(body=body, radiating_dof=dof, omega=1.0) for dof in body.dofs]
    solver = cpt.BEMSolver(method=method)
    results = solver.solve_all(problems, keep_details=True)
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
