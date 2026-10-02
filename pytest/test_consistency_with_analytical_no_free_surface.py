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


def test_potential_and_velocity_at_field_points_of_translating_sphere():
    # A sphere of radius a translating along x in an unbounded fluid at the
    # velocity U generates the dipole potential
    #     phi = -U a^3 x / (2 r^3).
    # (Capytaine's radiation potentials include an additional factor -i omega.)
    a = 1.0
    mesh = cpt.mesh_sphere(center=(0, 0, 0), radius=a, resolution=(40, 40))
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)))
    pb = cpt.RadiationProblem(body=body, free_surface=np.inf, radiating_dof="Surge", omega=1.5)
    solver = cpt.BEMSolver()
    res = solver.solve(pb, keep_details=True)

    points = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 1.0], [1.0, 1.0, 3.0], [-1.5, 0.5, 0.5]])
    x, y, z = points.T
    r = np.linalg.norm(points, axis=1)
    phi = -a**3 * x / (2*r**3)
    grad_phi = -a**3/2 * np.stack([1/r**3 - 3*x**2/r**5, -3*x*y/r**5, -3*x*z/r**5], axis=1)

    potential = solver.compute_potential(points, res)
    velocity = solver.compute_velocity(points, res)
    np.testing.assert_allclose(potential/(-1j*pb.omega), phi, rtol=5e-2, atol=1e-3)
    np.testing.assert_allclose(velocity/(-1j*pb.omega), grad_phi, rtol=5e-2, atol=1e-3)


@pytest.mark.parametrize("diagonal_term_in_double_layer, expected", [(False, 0.5), (True, 1.0)])
def test_known_values_of_double_layer_matrix_on_closed_sphere(diagonal_term_in_double_layer, expected):
    # Without free surface, the Green function is the Rankine source only.
    # By Gauss's law, the double layer matrix D applied to a constant density
    # on a closed surface is
    # - 1/2 for its principal value part (half of the solid angle),
    # - 1 when the diagonal term I/2 is included.
    mesh = cpt.mesh_sphere(radius=1.0, center=(0, 0, 0), resolution=(20, 20))
    _, D = cpt.Delhommeau().evaluate(
        mesh, mesh, free_surface=np.inf, adjoint_double_layer=False,
        diagonal_term_in_double_layer=diagonal_term_in_double_layer,
    )
    np.testing.assert_allclose(D @ np.ones(mesh.nb_faces), expected, atol=1e-2)
