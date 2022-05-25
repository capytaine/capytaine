#!/usr/bin/env python

import numpy as np
import capytaine as cpt

# Load a mesh
original_mesh = cpt.Sphere(axial_symmetry=False).mesh

# METHOD 1
# Get the indices of the panels such that the center is in the half-space y > 0.
indices_of_some_faces = np.where(original_mesh.faces_centers[:, 1] > 0)[0]
half_mesh = original_mesh.extract_faces(indices_of_some_faces)

# METHOD 2
# Clip the mesh
half_mesh = original_mesh.clip(plane=cpt.xOz_Plane)

# Define a new body by reflecting the half mesh across the xOz plane.
body = cpt.FloatingBody(
    cpt.ReflectionSymmetricMesh(half_mesh, cpt.xOz_Plane),
    name="symmetric_body")

# Display the body with VTK.
body.show()
