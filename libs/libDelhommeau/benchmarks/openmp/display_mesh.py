import numpy as np
import capytaine as cpt

vertices = np.loadtxt("vertices.dat")
normals = np.loadtxt("normal.dat")
mesh = cpt.Mesh(vertices, np.arange(0, vertices.shape[0]).reshape(-1, 4))
mesh.show()
