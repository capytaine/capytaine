import pytest

import numpy as np
import capytaine as cpt

from capytaine.bem.solver import BEMSolver
from capytaine.green_functions.delhommeau import Delhommeau
from capytaine.bem.engines import BasicMatrixEngine
from capytaine.bem.problems_and_results import RadiationProblem
from capytaine.bodies.predefined.spheres import Sphere

sphere = Sphere(radius=1.0, ntheta=2, nphi=3, clip_free_surface=True)
sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")

def test_gradiant_G_shape():
    mesh = cpt.mesh_sphere(radius=1, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
    S, gradG_1 = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.infty, 1.0, early_dot_product=False, method='direct')
    assert gradG_1.shape == (mesh.nb_faces, mesh.nb_faces, 3)

def test_gradient_G_a_posteriori_scalar_product():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.infty, 1.0, early_dot_product=False, method='direct')
    S, K = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.infty, 1.0, early_dot_product=True, method='direct')

    # Scalar product with the faces normal vectors
    K_ = np.zeros_like(S)
    for i in range(mesh.nb_faces):
        for j in range(mesh.nb_faces):
            K_[i, j] = gradG[i, j, :] @ mesh.faces_normals[j, :]

    # Shorter implementiation of the exact same scalar products
    K__ = np.einsum('...k,...k->...', gradG, mesh.faces_normals)

    np.testing.assert_allclose(K, K_)
    np.testing.assert_allclose(K, K__)
    np.testing.assert_allclose(K_, K__)

def test_gradient_G_diagonal_term():
    mesh = cpt.mesh_sphere(radius=1, center=(0, 0, 0), resolution=(10, 10)).immersed_part()
    # Passing two different Python objects
    _, gradG_1 = cpt.Delhommeau().evaluate(mesh.copy(), mesh, 0.0, np.infty, 1.0, early_dot_product=False, method='direct')
    # Passing the same Python object
    _, gradG_2 = cpt.Delhommeau().evaluate(mesh, mesh, 0.0, np.infty, 1.0, early_dot_product=False, method='direct')

    diag_normal = np.zeros_like(gradG_2)
    for i in range(mesh.nb_faces):
        diag_normal[i, i, :] = mesh.faces_normals[i, :]

    np.testing.assert_allclose(gradG_1, gradG_2 - 0.5*diag_normal, atol=1e-10, rtol=1e-10)

