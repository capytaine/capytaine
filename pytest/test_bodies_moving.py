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
"""Tests of the expansion of the motion of a body and of the integration of a pressure on the moving body."""

import numpy as np

import capytaine as cpt
import capytaine_test_helpers as helpers
from capytaine.tools.perturbation_expansion import DifferenceFrequencyExpansion as Expansion

RNG = np.random.default_rng(seed=0)


def skew(a):
    """Matrices S(a) such that S(a) @ v = a × v, for a of shape (..., 3)."""
    s = np.zeros(a.shape[:-1] + (3, 3), dtype=a.dtype)
    s[..., 0, 1] = -a[..., 2]
    s[..., 0, 2] = a[..., 1]
    s[..., 1, 0] = a[..., 2]
    s[..., 1, 2] = -a[..., 0]
    s[..., 2, 0] = -a[..., 1]
    s[..., 2, 1] = a[..., 0]
    return s


def difference_frequency_product(a, b):
    """Q_kl = 1/4 (a_k conj(b_l) + conj(a_l) b_k) for scalar amplitudes a and b of shape (nb_dir,)."""
    return (a[:, None] * np.conj(b[None, :]) + np.conj(a[None, :]) * b[:, None]) / 4


def second_order_rotation_matrix(rotation):
    """Difference-frequency QTF of the second-order term of R = R_roll R_pitch R_yaw (convention of WAMIT).

    For small angles θ, R = I + S(θ) + 1/2 Σ_j θ_j² S_j² + Σ_{j<k} θ_j θ_k S_j S_k + ..., with S_j = S(e_j).
    Each product θ_j θ_k is replaced by its difference-frequency QTF.

    rotation: array of shape (nb_dir, 3), complex amplitudes of roll, pitch and yaw.
    Returns an array of shape (nb_dir, nb_dir, 3, 3).
    """
    S = skew(np.eye(3))  # S[j] = S(e_j)
    nb_dir = rotation.shape[0]
    R2 = np.zeros((nb_dir, nb_dir, 3, 3), dtype=complex)
    for j in range(3):
        for k in range(j, 3):
            weight = 0.5 if j == k else 1.0
            R2 += weight * difference_frequency_product(rotation[:, j], rotation[:, k])[:, :, None, None] * (S[j] @ S[k])
    return R2


def rigid_body_and_motion(nb_dir=3):
    rotation_center = np.array([0.1, 0.0, -0.2])
    body = cpt.FloatingBody(mesh=helpers.small_sphere_mesh(), dofs=cpt.rigid_body_dofs(rotation_center=rotation_center))
    q = RNG.normal(size=(nb_dir, 6)) + 1j * RNG.normal(size=(nb_dir, 6))  # (Surge, Sway, Heave, Roll, Pitch, Yaw)
    return body, q, rotation_center


def test_second_order_rigid_body_motion_and_nanson_cofactor():
    # For rigid body motions, the order 2 of the displacement and the Nanson cofactor det(F) F^-T are identical.
    # This test check the values computed by the code for both against the reference `second_order_rotation_matrix`.
    body, q, rotation_center = rigid_body_and_motion()
    u, grad_u = body._displacement_expansions(q)
    deformation = np.eye(3) + grad_u
    cofactor = deformation.det().expand_dims(-1).expand_dims(-1) * deformation.inv().T
    R2 = second_order_rotation_matrix(q[:, 3:])
    np.testing.assert_allclose(u.coefficients[2], (body.mesh.faces_centers - rotation_center) @ np.swapaxes(R2, -1, -2), atol=1e-12)
    np.testing.assert_allclose(cofactor.coefficients[2], np.broadcast_to(R2[:, :, None], R2.shape[:2] + (body.mesh.nb_faces, 3, 3)), atol=1e-12)


def test_integrate_second_order_pressure_on_moving_body():
    # A pressure with only an order 2 is integrated on the mean body.
    body, q, _ = rigid_body_and_motion()
    nb_dir, nb_faces = q.shape[0], body.mesh.nb_faces
    p2 = RNG.normal(size=(nb_dir, nb_dir, nb_faces)) + 1j * RNG.normal(size=(nb_dir, nb_dir, nb_faces))
    p2 = (p2 + np.conj(np.swapaxes(p2, 0, 1))) / 2
    p = Expansion([0.0, np.zeros((nb_dir, nb_faces)), p2])
    second_order_force = body.integrate_pressure_on_moving_body(p, q).coefficients[2]
    for k in range(nb_dir):
        for l in range(nb_dir):
            np.testing.assert_allclose(second_order_force[k, l], list(body.integrate_pressure(p2[k, l]).values()), rtol=1e-12)
