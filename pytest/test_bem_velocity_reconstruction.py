
import pytest
import numpy as np
import capytaine as cpt

def test_full_K():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    mesh2 = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    S, K = cpt.Delhommeau().evaluate(mesh, mesh2, 0.0, -np.infty, 1.0, K_dim=1)
    S, gradG = cpt.Delhommeau().evaluate(mesh, mesh2, 0.0, -np.infty, 1.0, K_dim=3)
    K_ = np.zeros_like(K)
    D_ = np.zeros_like(K)
    for i in range(K.shape[0]):
        for j in range(K.shape[1]):
            K_[i, j] = gradG[i, j, :] @ mesh.faces_normals[i, :]
            D_[i, j] = gradG[i, j, :] @ mesh.faces_normals[j, :]
    K - K_
