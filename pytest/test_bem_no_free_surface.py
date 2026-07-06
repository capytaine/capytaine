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
import pytest

import capytaine as cpt


@pytest.mark.parametrize("z_center", [-10.0, 0.0, 10.0])
def test_analytical_solution(z_center):
    radius = 1.0
    mesh = cpt.mesh_sphere(center=(0, 0, z_center), radius=radius, resolution=(10, 10))
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    pb = cpt.RadiationProblem(body=body, free_surface=np.inf, radiating_dof="Surge")
    solver = cpt.BEMSolver(method="direct")
    res = solver.solve(pb)
    assert res.forces["Surge"] == pytest.approx(2/3*np.pi*pb.rho*radius**3, rel=1e-2)


def test_translation_invariance_of_no_free_surface_case():
    def force_on_body(z):
        mesh = cpt.mesh_parallelepiped(center=(0, 0, z))
        body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)))
        pb = cpt.RadiationProblem(body=body, free_surface=np.inf, water_depth=np.inf, radiating_dof="Surge")
        solver = cpt.BEMSolver(method="direct")
        res = solver.solve(pb)
        return res.force["Surge"]
    assert np.isclose(force_on_body(0.0), force_on_body(-1.0))
    assert np.isclose(force_on_body(0.0), force_on_body(1.0))
