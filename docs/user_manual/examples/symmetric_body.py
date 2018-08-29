#!/usr/bin/env python

import numpy as np

from meshmagick.geometry import xOz_Plane

from capytaine import FloatingBody
from capytaine.mesh.symmetries import ReflectionSymmetry

# Load a mesh from a file.
mesh_from_file = FloatingBody.from_file("example.dat").mesh

# Get the indices of the cells such that the center is in the half-space y > 0.
indices_of_some_faces = np.where(mesh_from_file.faces_centers[:, 1] > 0)[0]

# Create a new mesh with those faces.
half_mesh = mesh_from_file.extract_faces(indices_of_some_faces)

# Define a new body by reflecting the half mesh accros the xOz plane.
body = FloatingBody(ReflectionSymmetry(half_mesh, xOz_Plane), name="symmetric_body")

# Display the body with VTK.
body.show()
