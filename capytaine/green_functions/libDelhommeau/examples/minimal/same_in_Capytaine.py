import numpy as np
import capytaine as cpt

vert = np.array([0.0, 0.0, -1.0,
                 1.0, 0.0, -1.0,
                 1.0, 1.0, -1.0,
                 0.0, 1.0, -1.0,
                 1.0, 0.0, -1.0,
                 2.0, 0.0, -1.0,
                 2.0, 1.0, -1.0,
                 1.0, 1.0, -1.0]).reshape((8, 3))
faces = np.arange(0, 8).reshape(2, 4)
mesh = cpt.Mesh(vert, faces)
print("Rankine")
S, K = cpt.Delhommeau().evaluate(mesh, mesh, sea_bottom=-np.infty, free_surface=np.infty)
print(S)
print(K)

print("k=1.0")
S, K = cpt.Delhommeau().evaluate(mesh, mesh, sea_bottom=-np.infty, wavenumber=1.0)
print(S)
print(K)

print("k=2.0")
S, K = cpt.Delhommeau().evaluate(mesh, mesh, sea_bottom=-np.infty, wavenumber=2.0)
print(S)
print(K)
