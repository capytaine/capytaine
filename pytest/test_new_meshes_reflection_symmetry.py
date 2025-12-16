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
"""Tests for ReflectionSymmetricMesh class."""

from functools import lru_cache

import numpy as np
import pytest

from capytaine.new_meshes import (
    Mesh,
    ReflectionSymmetricMesh,
)


@lru_cache
def single_panel():
    vertices = np.array(
        [[0.5, 0.5, 0.0], [1.0, 0.5, 0.0], [1.0, 1.0, 0.0], [0.5, 1.0, 0.0]]
    )
    faces = np.array([[0, 1, 2, 3]])
    single_panel = Mesh(vertices=vertices, faces=faces)
    return single_panel


def test_create_symmetric_mesh_xoz():
    sym_mesh = ReflectionSymmetricMesh(half=single_panel(), plane="xOz")
    assert sym_mesh.plane == "xOz"
    assert sym_mesh.nb_faces == 2
    assert sym_mesh.nb_vertices == 8

def test_create_symmetric_mesh_yoz():
    sym_mesh = ReflectionSymmetricMesh(half=single_panel(), plane="yOz")
    assert sym_mesh.plane == "yOz"
    assert sym_mesh.nb_faces == 2
    assert sym_mesh.nb_vertices == 8

def test_invalid_plane_raises_error():
    with pytest.raises(ValueError, match="Plane must be 'xOz' or 'yOz'"):
        ReflectionSymmetricMesh(half=single_panel(), plane="invalid")

def test_merged_mesh_xoz():
    sym_mesh = ReflectionSymmetricMesh(half=single_panel(), plane="xOz")
    merged = sym_mesh.merged()

    # Should have 8 vertices (4 original + 4 reflected)
    assert merged.nb_vertices == 8
    assert merged.nb_faces == 2

    # Check that y-coordinates are reflected
    assert np.allclose(merged.faces_centers[0, 1], 0.75)  # Original
    assert np.allclose(merged.faces_centers[1, 1], -0.75)  # Reflected

def test_merged_mesh_yoz():
    sym_mesh = ReflectionSymmetricMesh(half=single_panel(), plane="yOz")
    merged = sym_mesh.merged()

    # Check that x-coordinates are reflected
    assert np.allclose(merged.faces_centers[0, 0], 0.75)  # Original
    assert np.allclose(merged.faces_centers[1, 0], -0.75)  # Reflected

def test_nested_symmetry_both_planes():
    inner_sym = ReflectionSymmetricMesh(half=single_panel(), plane="xOz")
    outer_sym = ReflectionSymmetricMesh(half=inner_sym, plane="yOz")

    assert outer_sym.nb_faces == 4  # 1 * 2 * 2
    assert outer_sym.nb_vertices == 16  # 4 * 2 * 2

    # Merge to get full mesh
    merged = outer_sym.merged()
    assert merged.nb_faces == 4
    assert merged.nb_vertices == 16

def test_merged_faces_have_correct_winding():
    """Test that reflected faces have reversed winding order."""
    sym_mesh = ReflectionSymmetricMesh(half=single_panel(), plane="xOz")
    merged = sym_mesh.merged()

    assert np.isclose(merged.faces_normals[0, 2], 1.0)
    assert np.isclose(merged.faces_normals[1, 2], 1.0)

def test_create_symmetric_mesh_with_inherited_metadata():
    mesh = single_panel().with_metadata(foo=[1.0])
    sym_mesh = ReflectionSymmetricMesh(half=mesh, plane="xOz")
    assert np.allclose(sym_mesh.faces_metadata["foo"], [1.0, 1.0])
    merged = sym_mesh.merged()
    assert np.allclose(merged.faces_metadata["foo"], [1.0, 1.0])

def test_create_symmetric_mesh_with_own_metadata():
    mesh = single_panel().with_metadata(foo=[1.0])
    metadata = {'bar': [1.0, 2.0]}
    sym_mesh = ReflectionSymmetricMesh(half=mesh, plane="xOz", faces_metadata=metadata)
    assert np.allclose(sym_mesh.faces_metadata["foo"], [1.0, 1.0])
    assert np.allclose(sym_mesh.faces_metadata["bar"], [1.0, 2.0])
    merged = sym_mesh.merged()
    assert np.allclose(merged.faces_metadata["foo"], [1.0, 1.0])
    assert np.allclose(merged.faces_metadata["bar"], [1.0, 2.0])

def test_join_same_symmetry_xoz():
    """Join two ReflectionSymmetricMesh with same xOz symmetry."""
    # Create two half meshes
    vertices1 = np.array([[0, 0.5, 0], [1, 0.5, 0], [1, 1, 0], [0, 1, 0]])
    faces1 = [[0, 1, 2, 3]]
    half1 = Mesh(vertices=vertices1, faces=faces1)
    sym1 = ReflectionSymmetricMesh(half=half1, plane="xOz")

    vertices2 = np.array([[2, 0.5, 0], [3, 0.5, 0], [3, 1, 0], [2, 1, 0]])
    faces2 = [[0, 1, 2, 3]]
    half2 = Mesh(vertices=vertices2, faces=faces2)
    sym2 = ReflectionSymmetricMesh(half=half2, plane="xOz")

    # Join them
    result = sym1 + sym2

    # Should preserve symmetry
    assert isinstance(result, ReflectionSymmetricMesh)
    assert result.plane == "xOz"
    assert result.nb_faces == 4  # 2 half faces * 2 with symmetry

def test_join_different_symmetries():
    """Join two ReflectionSymmetricMesh with different symmetries."""
    vertices1 = np.array([[0, 0.5, 0], [1, 0.5, 0], [1, 1, 0], [0, 1, 0]])
    faces1 = [[0, 1, 2, 3]]
    half1 = Mesh(vertices=vertices1, faces=faces1)
    sym1 = ReflectionSymmetricMesh(half=half1, plane="xOz")

    vertices2 = np.array([[0.5, 0, 0], [1, 0, 0], [1, 1, 0], [0.5, 1, 0]])
    faces2 = [[0, 1, 2, 3]]
    half2 = Mesh(vertices=vertices2, faces=faces2)
    sym2 = ReflectionSymmetricMesh(half=half2, plane="yOz")

    result = sym1 + sym2

    # Should return regular Mesh (no symmetry)
    assert isinstance(result, Mesh)
    assert not isinstance(result, ReflectionSymmetricMesh)
    assert result.nb_faces == 4  # 2 half faces * 2 with symmetry

def test_join_nested_symmetries():
    """Join two nested ReflectionSymmetricMesh with same symmetries."""
    # Create quarter meshes with both xOz and yOz symmetries
    vertices1 = np.array([[0.5, 0.5, 0], [1, 0.5, 0], [1, 1, 0], [0.5, 1, 0]])
    faces1 = [[0, 1, 2, 3]]
    quarter1 = Mesh(vertices=vertices1, faces=faces1)
    inner1 = ReflectionSymmetricMesh(half=quarter1, plane="xOz")
    sym1 = ReflectionSymmetricMesh(half=inner1, plane="yOz")

    vertices2 = np.array([[2.5, 0.5, 0], [3, 0.5, 0], [3, 1, 0], [2.5, 1, 0]])
    faces2 = [[0, 1, 2, 3]]
    quarter2 = Mesh(vertices=vertices2, faces=faces2)
    inner2 = ReflectionSymmetricMesh(half=quarter2, plane="xOz")
    sym2 = ReflectionSymmetricMesh(half=inner2, plane="yOz")

    result = sym1 + sym2

    # Should preserve nested symmetry
    assert isinstance(result, ReflectionSymmetricMesh)
    assert result.nb_faces == 8  # 2 quarter faces * 4 with both symmetries

def test_join_with_regular_mesh_with_warning():
    """Join ReflectionSymmetricMesh with regular Mesh."""
    vertices1 = np.array([[0, 0.5, 0], [1, 0.5, 0], [1, 1, 0], [0, 1, 0]])
    faces1 = [[0, 1, 2, 3]]
    half1 = Mesh(vertices=vertices1, faces=faces1)
    sym1 = ReflectionSymmetricMesh(half=half1, plane="xOz")

    vertices2 = np.array([[2, 0, 0], [3, 0, 0], [3, 1, 0], [2, 1, 0]])
    faces2 = [[0, 1, 2, 3]]
    mesh2 = Mesh(vertices=vertices2, faces=faces2)

    result = sym1 + mesh2

    # Should return regular Mesh
    assert isinstance(result, Mesh)
    assert not isinstance(result, ReflectionSymmetricMesh)

def test_join_with_metadata():
    sym1 = ReflectionSymmetricMesh(
        half=single_panel().with_metadata(center=[0.75]),
        faces_metadata={'is_ghost': [False, True]},
        plane="xOz"
    )
    sym2 = ReflectionSymmetricMesh(
        half=single_panel().translated_x(3.0).with_metadata(center=[3.75]),
        faces_metadata={'is_ghost': [False, True]},
        plane="xOz"
    )
    result = sym1 + sym2
    # First the metadata on result.half, then on result.other_half
    assert np.allclose(result.faces_metadata['center'], result.faces_centers[:, 0])
    assert np.allclose(result.faces_metadata['is_ghost'], [False, False, True, True])

def test_clipped_symmetric_mesh():
    vertices = np.array([[0, 0.5, -0.5], [0, 0.5, 0.5], [0, 1.5, 0.5], [0, 1.5, -0.5]])
    faces = [[0, 1, 2, 3]]
    half = Mesh(vertices=vertices, faces=faces, faces_metadata={"foo": [1.0]})
    sym = ReflectionSymmetricMesh(half=half, plane="xOz", faces_metadata={'bar': [2.0, 3.0]})
    imm_sym = sym.immersed_part()
    assert np.isclose(imm_sym.faces_areas.sum(), 0.5*sym.faces_areas.sum())
    assert imm_sym.faces_metadata['foo'].shape[0] == imm_sym.nb_faces
    assert np.allclose(imm_sym.faces_metadata['bar'][imm_sym.faces_centers[:, 1] > 0.0], 2.0)
