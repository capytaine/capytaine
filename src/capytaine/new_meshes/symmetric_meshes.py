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

from __future__ import annotations

from typing import Optional
from functools import cached_property
from itertools import chain

import numpy as np

from .abstract_meshes import AbstractMesh
from .meshes import Mesh

class ReflectionSymmetricMesh(AbstractMesh):
    """A mesh with reflection symmetry across a plane.

    This class represents a mesh that has reflection symmetry across either
    the xOz plane (y=0) or yOz plane (x=0). Only half of the mesh is stored,
    and the full mesh can be reconstructed by reflecting across the symmetry plane.

    Supports nested symmetries: if the half mesh is itself a ReflectionSymmetricMesh,
    this represents a quarter mesh with symmetries across both planes.

    Attributes
    ----------
    half: AbstractMesh
        The half mesh
    plane: str
        The symmetry plane, either "xOz" or "yOz"
    faces_metadata: Dict[str, np.ndarray]
        Some arrays with the same first dimension (should be the number
        of faces of the whole mesh) storing some fields defined on all the
        faces of the mesh.
    name: str, optional
        Name for the mesh

    Examples
    --------
    >>> # Create a mesh with xOz symmetry (y=0 plane)
    >>> half_mesh = Mesh(vertices=..., faces=...)
    >>> symmetric_mesh = ReflectionSymmetricMesh(half=half_mesh, plane="xOz")
    >>>
    >>> # Create a mesh with both xOz and yOz symmetries (quarter mesh)
    >>> quarter_mesh = Mesh(vertices=..., faces=...)
    >>> sym_xOz = ReflectionSymmetricMesh(half=quarter_mesh, plane="xOz")
    >>> sym_both = ReflectionSymmetricMesh(half=sym_xOz, plane="yOz")
    >>>
    >>> # Get the full merged mesh
    >>> full_mesh = symmetric_mesh.merged()
    """

    def __init__(
            self,
            half: AbstractMesh, *,
            plane: str,
            faces_metadata: Optional[Dict[str, np.ndarray]] = None,
            name: Optional[str] = None
            ):

        if plane not in ["xOz", "yOz"]:
            raise ValueError(f"Plane must be 'xOz' or 'yOz', got '{plane}'")

        self.half = half
        self.plane = plane
        self.other_half = self.half.mirrored(plane)

        self.faces_metadata = {k: np.concatenate([v, v]) for k, v in half.faces_metadata.items()}
        if faces_metadata is not None:
            self.faces_metadata.update(**{k: np.asarray(faces_metadata[k]) for k in faces_metadata})

        for m in self.faces_metadata:
            assert self.faces_metadata[m].shape[0] == self.nb_faces

        self.name = str(name) if name is not None else None

    @cached_property
    def nb_vertices(self) -> int:
        return 2*self.half.nb_vertices

    @cached_property
    def nb_faces(self) -> int:
        return 2*self.half.nb_faces

    @cached_property
    def vertices(self) -> np.ndarray:
        return np.concatenate([self.half.vertices, self.other_half.vertices])

    @cached_property
    def faces(self) -> np.ndarray:
        return np.concatenate([self.half.faces, self.other_half.faces])

    @cached_property
    def faces_normals(self) -> np.ndarray:
        return np.concatenate([self.half.faces_normals, self.other_half.faces_normals])

    @cached_property
    def faces_areas(self) -> np.ndarray:
        return np.concatenate([self.half.faces_areas, self.other_half.faces_areas])

    @cached_property
    def faces_centers(self) -> np.ndarray:
        return np.concatenate([self.half.faces_centers, self.other_half.faces_centers])

    @cached_property
    def faces_radiuses(self) -> np.ndarray:
        return np.concatenate([self.half.faces_radiuses, self.other_half.faces_radiuses])

    @cached_property
    def quadrature_points(self) -> np.ndarray:
        return (
                np.concatenate([self.half.quadrature_points[0], self.other_half.quadrature_points[0]]),
                np.concatenate([self.half.quadrature_points[1], self.other_half.quadrature_points[1]]),
                )

    def __str__(self) -> str:
        return (f"ReflectionSymmetricMesh(half={str(self.half)}"
                + f", plane='{self.plane}'"
                + (f", name=\"{self.name}\")" if self.name is not None else ")"))

    def __short_str__(self) -> str:
        return self.__str__()

    def __repr__(self) -> str:
        return self.__str__()

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    def extract_faces(self, faces_id, *, name=None) -> Mesh:
        return self.merged().extract_faces(faces_id, name=name)

    def translated(self, shift, *, name=None) -> Union[ReflectionSymmetricMesh, Mesh]:
        if ((self.plane == 'xOz' and abs(shift[1]) < 1e-12)
            or(self.plane == 'yOz' and abs(shift[0]) < 1e-12)):
            return ReflectionSymmetricMesh(
                half=self.half.translated(shift),
                plane=self.plane,
                faces_metadata=self.faces_metadata,
                name=name
            )
        else:
            return self.merged().translated(shift, name=name)

    def rotated_with_matrix(self, R, *, name=None) -> Mesh:
        return self.merged().rotated_with_matrix(R, name=name)

    def mirrored(self, plane: Literal['xOz', 'yOz'], *, name=None) -> ReflectionSymmetricMesh:
        return ReflectionSymmetricMesh(
            half=self.half.mirrored(plane),
            plane=self.plane,
            faces_metadata=self.faces_metadata,
            name=name
        )

    def join_meshes(self, *meshes, return_masks=False, name=None) -> Union[ReflectionSymmetricMesh, Mesh]:
        if (all(isinstance(m, ReflectionSymmetricMesh) for m in meshes) and
            all(m.plane == self.plane for m in meshes)):
                if return_masks:
                    joined_halves, half_masks = self.half.join_meshes(
                        *[m.half for m in meshes],
                        return_masks=True
                    )
                    masks = [np.concatenate([half_mask, half_mask]) for half_mask in half_masks]
                else:
                    joined_halves = self.half.join_meshes(
                        *[m.half for m in meshes],
                        return_masks=False
                    )

                faces_metadata = {k: np.concatenate(
                    [m.faces_metadata[k][:m.nb_faces//2] for m in chain([self], meshes)]
                    + [m.faces_metadata[k][m.nb_faces//2:] for m in chain([self], meshes)],
                    axis=0
                    )
                                  for k in AbstractMesh._common_metadata_keys(*meshes)}

                joined_mesh = ReflectionSymmetricMesh(
                        half=joined_halves,
                        plane=self.plane,
                        faces_metadata=faces_metadata,
                        name=name,
                        )

                if return_masks:
                    return joined_mesh, masks
                else:
                    return joined_mesh

        else:
            return Mesh.join_meshes(
                self.merged(),
                *[m.merged() for m in meshes],
                return_masks=return_masks,
                name=name
            )

    def with_normal_vector_going_down(self, **kwargs) -> ReflectionSymmetricMesh:
        return ReflectionSymmetricMesh(
                half=self.half.with_normal_vector_going_down(),
                plane=self.plane,
                faces_metadata=self.faces_metadata,
                name=self.name)

    def copy(self, *, faces_metadata=None, name=None) -> ReflectionSymmetricMesh:
        if faces_metadata is None:
            faces_metadata = self.faces_metadata.copy()
        if name is None:
            name = self.name
        return ReflectionSymmetricMesh(
                half=self.half.copy(),
                plane=self.plane,
                faces_metadata=faces_metadata,
                name=self.name)

    def merged(self, name=None) -> Mesh:
        return Mesh.join_meshes(
            self.half.merged(),
            self.other_half.merged()
        ).with_metadata(
            **self.faces_metadata
        )

    def clipped(self, *, origin, normal, name=None) -> Union[ReflectionSymmetricMesh, Mesh]:
        if ((self.plane == 'xOz' and abs(normal[0]) < 1e-12)
            or(self.plane == 'yOz' and abs(normal[1]) < 1e-12)):
            clipped_half, indices = (
                    self.half
                    .with_metadata(index=np.arange(self.half.nb_faces))
                    .clipped(origin=origin, normal=normal)
                    .pop_metadata("index")
                    )
            all_indices = np.concatenate([indices, indices + self.half.nb_faces])
            metadata = {k: self.faces_metadata[k][all_indices] for k in self.faces_metadata.keys()}
            return ReflectionSymmetricMesh(
                    half=clipped_half,
                    plane=self.plane,
                    faces_metadata=metadata,
                    name=name)
        else:
            LOG.warning("Dropping mesh reflection symmetry with respect to"
                        f"{self.plane} when clipping with respect to plane"
                        f"with origin {origin} and normal {normal}")
            return self.merged().clipped(origin=origin, normal=normal, name=name)

    def show(self, *, backend=None, ghost_meshes=None, **kwargs):
        if ghost_meshes is None:
            ghost_meshes = []
        ghost_meshes = ghost_meshes + [self.other_half.merged()]
        return self.half.show(backend=backend, ghost_meshes=ghost_meshes, **kwargs)


# For some backward compatibility:
yOz_Plane = "yOz"
xOz_Plane = "xOz"
