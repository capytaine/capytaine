#!/usr/bin/env python

import numpy as np

from capytaine import FloatingBody
from capytaine.tools.geometry import xOz_Plane
from capytaine.mesh.mesh_clipper import MeshClipper
from capytaine.mesh.symmetries import ReflectionSymmetry

# Load a mesh from a file.
original_body = FloatingBody.from_file("example.dat").mesh

# METHOD 1
# Get the indices of the panels such that the center is in the half-space y > 0.
indices_of_some_faces = np.where(original_body.mesh.faces_centers[:, 1] > 0)[0]
half_mesh = original_body.extract_faces(indices_of_some_faces)

# METHOD 2
# Clip the mesh
half_mesh = MeshClipper(original_body.mesh, plane=xOz_Plane).clipped_mesh

# Define a new body by reflecting the half mesh across the xOz plane.
body = FloatingBody(ReflectionSymmetry(half_mesh, xOz_Plane), name="symmetric_body")

# Display the body with VTK.
body.show()
