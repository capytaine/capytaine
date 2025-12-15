# Copyright 2025 Mews Labs
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

from capytaine.new_meshes.meshes import Mesh


def test_clip_above_geometry():
    """Test clipping above the geometry (no faces should be removed)"""

    # Define a simple mesh (a square)
    vertices = [[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]]
    faces = [[0, 1, 2, 3]]
    mesh = Mesh(vertices=vertices, faces=faces)

    origin = np.array([1.5, 0, 1.5])
    normal = np.array([0, 0, 1])

    clipped_mesh = mesh.clipped(origin=origin, normal=normal)

    # Check that the faces have not been removed
    assert (
        len(clipped_mesh.faces) == 1
    ), "Expected 1 face to remain after clipping, got {}".format(
        len(clipped_mesh.faces)
    )  # The mesh should remain unchanged
    assert np.array_equal(
        clipped_mesh.faces, faces
    ), "Expected faces to remain unchanged, got {}".format(clipped_mesh.faces)
    print("Test 1 (Clipping above the geometry) passed.")


def test_clip_below_geometry():
    """Test clipping below the geometry (all faces should be removed)"""

    # Define a simple mesh (a square)
    vertices = [[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]]
    faces = [[0, 1, 2, 3]]
    mesh = Mesh(vertices=vertices, faces=faces)

    # Define the clipping plane below the mesh (normal = +z)
    origin = np.array([-0.5, 0, -0.5])  # A point above the mesh
    normal = np.array([0, 0, 1])

    clipped_mesh = mesh.clipped(origin=origin, normal=normal)

    # Check that all faces have been removed (the mesh should be empty)
    assert (
        clipped_mesh.nb_faces == 0
    ), "Expected no faces to remain after clipping, got {}".format(
        clipped_mesh.nb_faces
    )
    assert (
        clipped_mesh.nb_vertices == 0
    ), "Expected no vertices to remain after clipping, got {}".format(
        clipped_mesh.nb_vertices
    )
    print("Test 2 (Clipping below the geometry) passed.")


def test_clip_partial():
    """Test partial clipping (some faces should be cut by the plane)"""

    # Define a simple mesh (a square face)
    vertices = [[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]]
    faces = [[0, 1, 2, 3]]
    mesh = Mesh(vertices=vertices, faces=faces)
    print(mesh.faces)

    # Define the clipping plane that partially cuts the face (normal = +z)
    origin = np.array([0.5, 0.0, 0.5])
    normal = np.array([0, 0, 1])

    # Check that the face was split into two triangles
    clipped_mesh = mesh.clipped(origin=origin, normal=normal)

    # Check that the face was split into two triangles
    assert (
        clipped_mesh.nb_triangles == 2
    ), f"Expected 2 triangles, got {clipped_mesh.nb_triangles}"
    assert clipped_mesh.nb_quads == 0, f"Expected 0 quads, got {clipped_mesh.nb_quads}"
    print("Test 3 (Partial clipping) passed.")


def test_clip_partial_2():
    """Test partial clipping (some faces should be cut by the plane)"""

    # Define a simple mesh (a square face)
    vertices = [[0, 0, 0], [1, 0, 0], [1, 0, 2], [0, 0, 1]]
    faces = [[0, 1, 2, 3]]
    mesh = Mesh(vertices=vertices, faces=faces)

    # Define the clipping plane that partially cuts the face (normal = +z)
    origin = np.array([0.5, 0, 1.2])
    normal = np.array([0, 0, 1])

    clipped_mesh = mesh.clipped(origin=origin, normal=normal)

    # Check that the face was split into one triangle and one quad
    assert (
        clipped_mesh.nb_triangles == 1
    ), f"Expected 1 triangle, got {clipped_mesh.nb_triangles}"
    assert clipped_mesh.nb_quads == 1, f"Expected 1 quad, got {clipped_mesh.nb_quads}"
    print("Test 4 (Partial clipping 2) passed.")


def test_faces_on_the_free_surface():
    vertices = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]
    faces = [[0, 1, 2, 3]]
    mesh = Mesh(vertices=vertices, faces=faces)
    assert mesh.immersed_part().nb_faces == 1
    assert mesh.translated_z(-1e-10).immersed_part().nb_faces == 1
    assert mesh.translated_z(1e-9).immersed_part().nb_faces == 1  # Within tol of 1e-8
    assert mesh.translated_z(1e-7).immersed_part().nb_faces == 0


def test_clip_metadata():
    vertices = [[0, 0, -1], [0, 0, 1], [1, 0, 1], [1, 0, -1]]
    faces = [[0, 1, 2, 3]]
    mesh = Mesh(vertices=vertices, faces=faces, faces_metadata={'foo': ['a']})
    imm_mesh = mesh.immersed_part()
    assert np.all(imm_mesh.faces_metadata['foo'] == np.array(['a', 'a']))
