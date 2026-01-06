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
"""Tests for RotationSymmetricMesh class."""

from functools import lru_cache

import numpy as np
import pytest
import capytaine as cpt

from capytaine.new_meshes import (
    Mesh,
    RotationSymmetricMesh,
    ReflectionSymmetricMesh
)
from capytaine.new_meshes.meshes import to_new_mesh


@lru_cache
def single_panel():
    vertices = np.array(
        [[0.5, 0.0, 0.0], [0.5, 0.0, -0.5], [0.5, 0.5, -0.5], [0.5, 0.5, 0.0]]
    )
    faces = np.array([[0, 1, 2, 3]])
    single_panel = Mesh(vertices=vertices, faces=faces)
    return single_panel


@pytest.mark.parametrize("n", [3, 4, 6])
def test_create_rotation_symmetric_mesh_sphere(n):
    """Create a sphere with 3-fold rotation symmetry and reconstruct it."""
    sphere = to_new_mesh(cpt.mesh_sphere(radius=1.0, center=(0, 0, 0), resolution=(12, 12)))
    wedge = sphere.extract_wedge(n=n)
    rot_sym = RotationSymmetricMesh(wedge=wedge, n=n)

    assert rot_sym.n == n
    assert rot_sym.axis == "z+"
    assert rot_sym.wedge.nb_faces == wedge.nb_faces

    merged = rot_sym.merged()

    # Check that merged mesh has at least as many faces as original
    # (clipping creates boundary faces, so reconstructed may have more)
    assert (
        merged.nb_faces >= sphere.nb_faces * 0.9
    )  # Allow 10% tolerance for degenerate faces
    assert merged.nb_vertices > 0

    # Check that reconstruction creates a reasonable mesh
    assert merged.nb_faces > wedge.nb_faces  # Should have more than just the wedge

def test_invalid_rotation_order():
    """Test that n < 2 raises ValueError."""
    with pytest.raises(ValueError, match="Rotation order must be >= 2"):
        RotationSymmetricMesh(wedge=single_panel(), n=1)

def test_repr():
    rot_sym = RotationSymmetricMesh(wedge=single_panel(), n=3)

    repr_str = repr(rot_sym)
    assert "RotationSymmetricMesh" in repr_str
    assert "n=3" in repr_str

def test_create_from_profile_bottom_to_top():
    meridian_points = np.array([(np.sqrt(1-z**2), 0.0, z) for z in np.linspace(-1.0, 1.0, 5)])
    sphere = RotationSymmetricMesh.from_profile_points(meridian_points, n=5)
    assert sphere.wedge.nb_faces == 4
    assert sphere.nb_faces == 4*5
    assert np.all(np.array([sphere.faces_normals[i, :] @ sphere.faces_centers[i, :] for i in range(sphere.nb_faces)]) >= 0.0)  # Outwards normals

def test_create_from_profile_top_to_bottom():
    meridian_points = np.array([(np.sqrt(1-z**2), 0.0, z) for z in np.linspace(1.0, -1.0, 5)])
    sphere = RotationSymmetricMesh.from_profile_points(meridian_points, n=5)
    assert sphere.wedge.nb_faces == 4
    assert sphere.nb_faces == 4*5
    assert np.all(np.array([sphere.faces_normals[i, :] @ sphere.faces_centers[i, :] for i in range(sphere.nb_faces)]) >= 0.0)  # Outwards normals

def test_create_symmetric_mesh_with_inherited_metadata():
    mesh = single_panel().with_metadata(foo=[1.0])
    sym_mesh = RotationSymmetricMesh(wedge=mesh, n=4)
    assert np.allclose(sym_mesh.faces_metadata["foo"], [1.0, 1.0, 1.0, 1.0])
    merged = sym_mesh.merged()
    assert np.allclose(merged.faces_metadata["foo"], [1.0, 1.0, 1.0, 1.0])

def test_create_symmetric_mesh_with_own_metadata():
    mesh = single_panel().with_metadata(foo=[1.0])
    metadata = {'bar': [1.0, 2.0, 3.0, 4.0]}
    sym_mesh = RotationSymmetricMesh(wedge=mesh, n=4, faces_metadata=metadata)
    assert np.allclose(sym_mesh.faces_metadata["foo"], [1.0, 1.0, 1.0, 1.0])
    assert np.allclose(sym_mesh.faces_metadata["bar"], [1.0, 2.0, 3.0, 4.0])
    merged = sym_mesh.merged()
    assert np.allclose(merged.faces_metadata["foo"], [1.0, 1.0, 1.0, 1.0])
    assert np.allclose(merged.faces_metadata["bar"], [1.0, 2.0, 3.0, 4.0])

def test_transform_keeping_symmetry():
    sym_mesh = RotationSymmetricMesh(wedge=single_panel(), n=4)
    rotated = sym_mesh.rotated_z(np.pi/4)
    assert isinstance(rotated, RotationSymmetricMesh)

def test_merged_face_ordering():
    sym_mesh = RotationSymmetricMesh(wedge=single_panel(), n=4).with_metadata(foo=range(4))
    assert np.allclose(
            sym_mesh.merged().faces_centers,
            sym_mesh.faces_centers,
            )
    assert np.allclose(
            sym_mesh.merged().faces_metadata['foo'],
            sym_mesh.faces_metadata['foo'],
            )

@pytest.mark.xfail
def test_mirrored():
    sym_mesh = RotationSymmetricMesh(wedge=single_panel(), n=4).with_metadata(foo=range(4))
    assert np.allclose(
            sym_mesh.merged().mirrored("xOz").faces_centers,
            sym_mesh.mirrored("xOz").merged().faces_centers,
            )
    assert np.allclose(
            sym_mesh.merged().mirrored("xOz").faces_metadata['foo'],
            sym_mesh.mirrored("xOz").merged().faces_metadata['foo'],
            )

def test_join_same_symmetries():
    sym_mesh = RotationSymmetricMesh(wedge=single_panel(), n=4)
    joined = sym_mesh + sym_mesh.rotated_z(np.pi/4)
    assert isinstance(joined, RotationSymmetricMesh)

def test_transform_losing_symmetry():
    sym_mesh = RotationSymmetricMesh(wedge=single_panel(), n=4)
    rotated = sym_mesh.rotated_x(np.pi/4)
    assert isinstance(rotated, Mesh)

def test_join_with_regular_mesh():
    """Join ReflectionSymmetricMesh with regular Mesh."""
    sym1 = RotationSymmetricMesh(wedge=single_panel(), n=4)

    vertices2 = np.array([[2, 0, 0], [3, 0, 0], [3, 1, 0], [2, 1, 0]])
    faces2 = [[0, 1, 2, 3]]
    mesh2 = Mesh(vertices=vertices2, faces=faces2)

    result = sym1 + mesh2
    assert isinstance(result, Mesh)

def test_join_with_metadata():
    sym1 = RotationSymmetricMesh(
        wedge=single_panel().with_metadata(center=[-0.25]),
        n=4,
        faces_metadata={'is_ghost': [False, True, True, True]},
    )
    sym2 = RotationSymmetricMesh(
        wedge=single_panel().translated_z(3.0).with_metadata(center=[2.75]),
        n=4,
        faces_metadata={'is_ghost': [False, True, True, True]},
    )
    result = sym1 + sym2
    # First the metadata on result.half, then on result.other_half
    assert np.allclose(result.faces_metadata['center'], result.faces_centers[:, 2])
    assert np.allclose(result.faces_metadata['is_ghost'], [False, False, True, True, True, True, True, True])

def test_clipped_symmetric_mesh():
    mesh = single_panel().translated_z(0.25).with_metadata(foo=[1.0])
    sym = RotationSymmetricMesh(wedge=mesh, n=4, faces_metadata={'bar': [2.0, 2.0, 2.0, 2.0]})
    imm_sym = sym.immersed_part()
    assert np.isclose(imm_sym.faces_areas.sum(), 0.5*sym.faces_areas.sum())
    assert imm_sym.faces_metadata['foo'].shape[0] == imm_sym.nb_faces
    assert np.allclose(imm_sym.faces_metadata['bar'][imm_sym.faces_centers[:, 1] > 0.0], 2.0)

# NESTED SYMMETRIES

@lru_cache
def nested_symmetry_with_inner_reflection():
    inner_sym = ReflectionSymmetricMesh(half=single_panel(), plane="xOz")
    return RotationSymmetricMesh(wedge=inner_sym, n=4)

def test_nested_symmetry_with_inner_reflection():
    sym = nested_symmetry_with_inner_reflection()

    assert sym.nb_faces == 8  # 1 * 2 * 4
    assert sym.nb_vertices == 32  # 4 * 2 * 4

    # Merge to get full mesh
    merged = sym.merged()
    assert merged.nb_faces == 8

def test_transforming_symmetries_with_inner_reflection():
    sym = nested_symmetry_with_inner_reflection()
    t_sym = sym.translated_z(-1.0)
    assert isinstance(t_sym, RotationSymmetricMesh)
    assert isinstance(t_sym.wedge, ReflectionSymmetricMesh)
    assert np.all(t_sym.vertices[:, 2] <= 0.5)

    t_sym = sym.translated_x(1.0)
    assert isinstance(t_sym, Mesh)

def test_joining_nested_symmetries_with_inner_reflection():
    sym = nested_symmetry_with_inner_reflection()
    b_sym = sym + sym.translated_z(-1.0)
    assert isinstance(b_sym, RotationSymmetricMesh)
    assert b_sym.n == sym.n
    assert isinstance(b_sym.wedge, ReflectionSymmetricMesh)
    assert b_sym.wedge.plane == sym.wedge.plane

    b_sym = sym + sym.translated_x(1.0)
    assert isinstance(b_sym, Mesh)

def test_clipping_nested_symmetries_with_inner_reflection():
    sym = nested_symmetry_with_inner_reflection().translated_z(0.25)
    c_sym = sym.immersed_part()
    assert isinstance(c_sym, RotationSymmetricMesh)
    assert c_sym.n == sym.n
    assert isinstance(c_sym.wedge, ReflectionSymmetricMesh)
    assert c_sym.wedge.plane == sym.wedge.plane

def test_nested_symmetry_with_outer_reflection():
    inner_sym = RotationSymmetricMesh(wedge=single_panel(), n=4)
    outer_sym = ReflectionSymmetricMesh(half=inner_sym, plane="xOz")

    assert outer_sym.nb_faces == 8  # 1 * 2 * 4
    assert outer_sym.nb_vertices == 32  # 4 * 2 * 4

    # Merge to get full mesh
    merged = outer_sym.merged()
    assert merged.nb_faces == 8
