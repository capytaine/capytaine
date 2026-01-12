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

import logging
from functools import cached_property
from typing import List, Union, Tuple, Dict, Optional, Literal

import xarray as xr
import numpy as np

from .abstract_meshes import AbstractMesh
from .geometry import (
    compute_faces_areas,
    compute_faces_centers,
    compute_faces_normals,
    compute_faces_radii,
    get_vertices_face,
)
from .clip import clip_faces
from .clean import clean_mesh
from .export import export_mesh
from .quality import _is_valid, check_mesh_quality
from .visualization import show_3d

LOG = logging.getLogger(__name__)


class Mesh(AbstractMesh):
    """Mesh class for representing and manipulating 3D surface meshes.

    Parameters
    ----------
    vertices : np.ndarray, optional
        Array of mesh vertices coordinates with shape (n_vertices, 3).
        Each row represents one vertex's (x, y, z) coordinates.
    faces : List[List[int]] or np.ndarray, optional
        Array of mesh connectivities for panels. Each row contains indices
        of vertices that form a face (triangles or quads).
    faces_metadata: Dict[str, np.ndarray]
        Some arrays with the same first dimension (should be the number of faces)
        storing some fields defined on all the faces of the mesh.
    name : str, optional
        Optional name for the mesh instance.
    auto_clean : bool, optional
        Whether to automatically clean the mesh upon initialization. Defaults to True.
    auto_check : bool, optional
        Whether to automatically check mesh quality upon initialization. Defaults to True.

    Attributes
    ----------
    vertices : np.ndarray
        Array of vertex coordinates with shape (n_vertices, 3).
    name : str or None
        Name of the mesh instance.
    """

    def __init__(
        self,
        vertices: np.ndarray = None,
        faces: Union[List[List[int]], np.ndarray] = None,
        *,
        faces_metadata: Optional[Dict[str, np.ndarray]] = None,
        name: Optional[str] = None,
        auto_clean: bool = True,
        auto_check: bool = True,
    ):
        # --- Vertices: always a NumPy array with shape (n,3) ---
        if vertices is None:
            self.vertices = np.empty((0, 3), dtype=np.float64)
        else:
            self.vertices = np.array(vertices, dtype=np.float64)

        # --- Faces: process using helper method ---
        self._faces: List[List[int]] = self._process_faces(faces)

        if faces_metadata is None:
            self.faces_metadata = {}
        else:
            self.faces_metadata = {k: np.asarray(faces_metadata[k]) for k in faces_metadata}

        for m in self.faces_metadata:
            assert self.faces_metadata[m].shape[0] == len(self._faces)

        # Optional name
        self.name = str(name) if name is not None else None

        # Cleaning/quality (unless mesh is completely empty)
        if not (len(self.vertices) == 0 and len(self._faces) == 0):
            if not _is_valid(vertices, faces):
                raise ValueError(
                    "Mesh is invalid: faces contain out-of-bounds or negative indices."
                )

            if np.any(np.isnan(vertices)):
                raise ValueError(
                    "Mesh is invalid: vertices coordinates contains NaN values."
                )

            if auto_clean:
                self.vertices, self._faces, self.faces_metadata = clean_mesh(
                    self.vertices, self._faces, self.faces_metadata, max_iter=5, tol=1e-8
                )

            if auto_check:
                check_mesh_quality(self)

    ## MAIN METRICS AND DISPLAY

    @cached_property
    def nb_vertices(self) -> int:
        """Number of vertices in the mesh."""
        return len(self.vertices)

    @cached_property
    def nb_faces(self) -> int:
        """Number of faces in the mesh."""
        return len(self._faces)

    @cached_property
    def nb_triangles(self) -> int:
        """Number of triangular faces (3-vertex) in the mesh."""
        return sum(1 for f in self._faces if len(f) == 3)

    @cached_property
    def nb_quads(self) -> int:
        """Number of quadrilateral faces (4-vertex) in the mesh."""
        return sum(1 for f in self._faces if len(f) == 4)

    def summary(self):
        """Print a summary of the mesh properties.

        Notes
        -----
        Displays the mesh name, vertex count, face count, and bounding box.
        """
        print("Mesh Summary")
        print(f"  Name           : {self.name}")
        print(f"  Vertices count : {self.nb_vertices}")
        print(f"  Faces count    : {self.nb_faces}")
        print(
            f"  Bounding box     : {self.vertices.min(axis=0)} to {self.vertices.max(axis=0)}"
        )
        print(f"  Metadata keys  : {self.faces_metadata.keys()}")

    def __str__(self) -> str:
        return (f"Mesh(vertices=[[... {self.nb_vertices} vertices ...]], "
                + f"faces=[[... {self.nb_faces} faces ...]]"
                + (f", name=\"{self.name}\")" if self.name is not None else ")"))

    def __short_str__(self) -> str:
        if self.name is not None:
            return f"Mesh(..., name={self.name})"
        else:
            return "Mesh(...)"

    def __repr__(self) -> str:
        return self.__str__()

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    def __rich_repr__(self):
        class CustomRepr:
            def __init__(self, n, kind):
                self.n = n
                self.kind = kind
            def __repr__(self):
                return "[[... {} {} ...]]".format(self.n, self.kind)
        yield "vertices", CustomRepr(self.nb_vertices, "vertices")
        yield "faces", CustomRepr(self.nb_faces, "faces")
        yield "name", self.name

    def show(self, *, backend=None, **kwargs):
        """Visualize the mesh using the specified backend.

        Parameters
        ----------
        backend : str, optional
            Visualization backend to use. Options are 'pyvista' or 'matplotlib'.
            By default, try several until an installed one is found.
        normal_vectors: bool, optional
            If True, display normal vectors on each face.
        **kwargs
            Additional keyword arguments passed to the visualization backend.
            See :mod:`~capytaine.new_meshes.visualization`

        Returns
        -------
        object
            Visualization object returned by the backend (e.g., matplotlib figure).

        Raises
        ------
        NotImplementedError
            If the specified backend is not supported.
        """
        return show_3d(self, backend=backend, **kwargs)


    ## INITIALISATION

    @staticmethod
    def _has_leading_count_column(arr: np.ndarray) -> bool:
        """Check if a 2D array has a leading column containing vertex counts.

        Parameters
        ----------
        arr : np.ndarray
            2D array of face data

        Returns
        -------
        bool
            True if the first column appears to be vertex counts
        """
        if arr.ndim != 2 or arr.shape[1] <= 3:
            return False

        expected_count = arr.shape[1] - 1
        for row in arr:
            # Check if first value could be a vertex count (3 or 4)
            # and if it matches the expected count (total cols - 1)
            if row[0] != expected_count and row[0] not in [3, 4]:
                return False
        return True

    def _process_faces(self, faces: List[List[int]] | np.ndarray) -> List[List[int]]:
        """Process the faces input for the Mesh class.

        Parameters
        ----------
        faces : np.ndarray or list
            The faces data to process.

        Returns
        -------
        list
            A list of faces, where each face is a list of vertex indices.

        Notes
        -----
        If the input is a 2D array with a leading column containing face vertex counts
        (e.g., [[3, v1, v2, v3], [4, v1, v2, v3, v4]]), the count column will be
        automatically stripped. This is checked per-row to support mixed triangle/quad meshes.
        """
        if faces is None:
            return []
        elif isinstance(faces, list):
            # assume it's already a list of lists of ints
            return [list(f) for f in faces]
        else:
            # fallback: convert array → nested list
            arr = np.asarray(faces, dtype=int)

            # Detect & strip a leading "count" column if present
            if self._has_leading_count_column(arr):
                arr = arr[:, 1:]

            return arr.tolist()

    @classmethod
    def from_list_of_faces(cls, list_faces, *, faces_metadata=None, name=None, auto_clean=True, auto_check=True) -> "Mesh":
        """
        Create a Mesh instance from a list of faces defined by vertex coordinates.

        Parameters
        ----------
        list_faces : list of list of list of float
            Each face is defined by a list of 3D coordinates. For example:
            [
                [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]],
                [[x4, y4, z4], [x5, y5, z5], [x6, y6, z6]]
            ]
        faces_metadata: Optional[Dict[str, np.ndarray]]
        name: str, optional
            A name for the new mesh.
        auto_clean : bool, optional
            Whether to automatically clean the mesh upon initialization. Defaults to True.
        auto_check : bool, optional
            Whether to automatically check mesh quality upon initialization. Defaults to True.

        Returns
        -------
        Mesh
            An instance of Mesh with:
            - unique vertices extracted from the input
            - faces defined as indices into the vertex array
        """
        unique_vertices = []
        vertices_map = {}
        indexed_faces = []

        for face in list_faces:
            indexed_face = []
            for coord in face:
                key = tuple(coord)
                if key not in vertices_map:
                    vertices_map[key] = len(unique_vertices)
                    unique_vertices.append(coord)
                indexed_face.append(vertices_map[key])
            indexed_faces.append(indexed_face)

        return cls(
            vertices=np.array(unique_vertices),
            faces=indexed_faces,
            faces_metadata=faces_metadata,
            name=name
        )

    def as_list_of_faces(self) -> List[List[List[float]]]:
        """
        Convert the Mesh instance to a list of faces defined by vertex coordinates.

        Returns
        -------
        list of list of list of float
            Each face is defined by a list of 3D coordinates. For example:
            [
                [[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]],
                [[x4, y4, z4], [x5, y5, z5], [x6, y6, z6]]
            ]
        """
        list_faces = []
        for face in self._faces:
            face_coords = [self.vertices[idx].tolist() for idx in face]
            list_faces.append(face_coords)
        if len(self.faces_metadata) > 0:
            LOG.info(f"Dropping metadata of {self} to export as list of faces.")
        return list_faces

    def as_array_of_faces(self) -> np.ndarray:
        """Similare to as_list_of_faces but returns an array of shape
        (nb_faces, 3, 3) if only triangles, or (nb_faces, 4, 3) otherwise.
        """
        array = self.vertices[self.faces[:, :], :]
        if self.nb_quads == 0:
            array = array[:, :3, :]
        return array

    def export(self, format, **kwargs):
        return export_mesh(self, format, **kwargs)

    ## INTERFACE FOR BEM SOLVER

    @cached_property
    def faces_vertices_centers(self) -> np.ndarray:
        """Calculate the center of vertices that form the faces.

        Returns
        -------
        np.ndarray
            Array of shape (n_faces, 3) containing the centroid of each face's vertices.
        """
        centers_vertices = []
        for face in self._faces:
            if face[3] != face[2]:
                a, b, c, d = get_vertices_face(face, self.vertices)
                mean = (a + b + c + d) / 4
                centers_vertices.append(mean)
            else:
                a, b, c = get_vertices_face(face, self.vertices)
                mean = (a + b + c) / 3
                centers_vertices.append(mean)
        return np.array(centers_vertices)

    @cached_property
    def faces_normals(self) -> np.ndarray:
        """Normal vectors for each face.

        Returns
        -------
        np.ndarray
            Array of shape (n_faces, 3) containing unit normal vectors.
        """
        return compute_faces_normals(self.vertices, self._faces)

    @cached_property
    def faces_areas(self) -> np.ndarray:
        """Surface area of each face.

        Returns
        -------
        np.ndarray
            Array of shape (n_faces,) containing the area of each face.
        """
        return compute_faces_areas(self.vertices, self._faces)

    @cached_property
    def faces_centers(self) -> np.ndarray:
        """Geometric centers of each face.

        Returns
        -------
        np.ndarray
            Array of shape (n_faces, 3) containing the center point of each face.
        """
        return compute_faces_centers(self.vertices, self._faces)

    @cached_property
    def faces_radiuses(self) -> np.ndarray:
        """Radii of each face (circumradius or characteristic size).

        Returns
        -------
        np.ndarray
            Array of shape (n_faces,) containing the radius of each face.
        """
        return compute_faces_radii(self.vertices, self._faces)

    @cached_property
    def faces(self) -> np.ndarray:
        """Face connectivity as quadrilateral array.

        Returns
        -------
        np.ndarray
            Array of shape (n_faces, 4) where triangular faces are padded
            by repeating the last vertex.

        Notes
        -----
        This property converts all faces to a uniform quad representation
        for compatibility with libraries expecting fixed-width face arrays.
        """
        faces_as_quad = [f if len(f) == 4 else f + [f[-1]] for f in self._faces]
        return np.array(faces_as_quad, dtype=int)

    @cached_property
    def quadrature_points(self) -> Tuple[np.ndarray, np.ndarray]:
        """Quadrature points and weights for numerical integration.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            (points, weights) where points has shape (n_faces, 1, 3) and
            weights has shape (n_faces, 1), using face centers and areas.
        """
        return (self.faces_centers.reshape((-1, 1, 3)), self.faces_areas.reshape(-1, 1))

    ## TRANSFORMATIONS

    def extract_faces(self, faces_id, *, name=None) -> "Mesh":
        """Extract a subset of faces by their indices and return a new Mesh instance.

        Parameters
        ----------
        faces_id : array_like
            Indices of faces to extract.
        name: str, optional
            A name for the new mesh

        Returns
        -------
        Mesh
            New mesh containing only the specified faces.
        """
        if isinstance(faces_id, np.ndarray):
            faces_id = faces_id.ravel()
        all_faces = self.as_list_of_faces()
        selected_faces = [all_faces[i] for i in faces_id]
        return Mesh.from_list_of_faces(
            selected_faces,
            faces_metadata={k: self.faces_metadata[k][selected_faces, ...] for k in self.faces_metadata},
            name=name,
            auto_clean=False,
            auto_check=False
        )

    def translated(self, shift, *, name=None) -> "Mesh":
        """Return a new Mesh translated along vector-like `shift`."""
        return Mesh(
            vertices=self.vertices + np.asarray(shift),
            faces=self._faces,
            faces_metadata=self.faces_metadata,
            name=name,
            auto_clean=False,
            auto_check=False,
        )

    def rotated_with_matrix(self, R, *, name=None) -> "Mesh":
        """Return a new Mesh rotated using the provided 3×3 rotation matrix."""
        new_vertices = self.vertices @ R.T
        return Mesh(
            vertices=new_vertices,
            faces=self._faces,
            name=name,
            faces_metadata=self.faces_metadata,
            auto_clean=False,
            auto_check=False,
        )

    def mirrored(self, plane: Literal['xOz', 'yOz'], *, name=None):
        new_vertices = self.vertices.copy()
        if plane == "xOz":
            new_vertices[:, 1] *= -1
        elif plane == "yOz":
            new_vertices[:, 0] *= -1
        else:
            raise ValueError(f"Unsupported value for plane: {plane}")
        new_faces = [f[::-1] for f in self._faces]  # Invert normals
        if name is None and self.name is not None:
            name = f"mirrored_{self.name}"
        return Mesh(
            new_vertices,
            new_faces,
            faces_metadata=self.faces_metadata,
            name=name,
            auto_clean=False,
            auto_check=False
        )

    def join_meshes(*meshes: List["Mesh"], return_masks=False, name=None) -> "Mesh":
        """Join two meshes and return a new Mesh instance.

        Parameters
        ----------
        meshes: List[Mesh]
            Meshes to be joined
        return_masks: bool, optional
            If True, additionally return a list of numpy masks establishing the
            origin of each face in the new mesh.
            (Default: False)
        name: str, optional
            A name for the new object

        Returns
        -------
        Mesh
            New mesh containing vertices and faces from all meshes.

        See Also
        --------
        __add__ : Implements the + operator for mesh joining.
        """
        if not all(isinstance(m, Mesh) for m in meshes):
            raise TypeError("Only Mesh instances can be added together.")

        if return_masks:
            meshes = [m.with_metadata(origin_mesh_index=np.array([i]*m.nb_faces))
                      for i, m in enumerate(meshes)]

        faces = sum((m.as_list_of_faces() for m in meshes), [])

        faces_metadata = {k: np.concatenate([m.faces_metadata[k] for m in meshes], axis=0)
                             for k in AbstractMesh._common_metadata_keys(*meshes)}

        joined_mesh = Mesh.from_list_of_faces(faces, faces_metadata=faces_metadata, name=name)
        # If list of faces is trimmed for some reason, metadata will be updated accordingly

        if return_masks:
            masks = [joined_mesh.faces_metadata['origin_mesh_index'] == i for i in range(len(meshes))]
            return joined_mesh.without_metadata('origin_mesh_index'), masks
        else:
            return joined_mesh

    def with_normal_vector_going_down(self, **kwargs) -> "Mesh":
        # Kwargs are for backward compatibility with former inplace implementation of this.
        # It could be removed in the final release.
        """Ensure normal vectors point downward (negative z-direction).

        Returns
        -------
        Mesh
            Self if normals already point down, otherwise modifies face orientation.

        Notes
        -----
        Used for lid meshes to avoid irregular frequency issues by ensuring
        consistent normal vector direction.
        """
        # For lid meshes for irregular frequencies removal
        if np.allclose(self.faces_normals[:, 2], np.ones((self.nb_faces,))):
            # The mesh is horizontal with normal vectors going up
            LOG.warning(
                f"Inverting the direction of the normal vectors of {self} to be downward."
            )
            return Mesh(
                vertices=self.vertices,
                faces=self.faces[:, ::-1],
                faces_metadata=self.faces_metadata,
                name=self.name
            )
        else:
            return self

    def copy(self, *, faces_metadata=None, name=None) -> Mesh:
        # No-op for backward compatibility
        if faces_metadata is None:
            faces_metadata = self.faces_metadata.copy()
        if name is None:
            name = self.name
        return Mesh(
            vertices=self.vertices,
            faces=self._faces,
            faces_metadata=faces_metadata,
            name=name,
            auto_clean=False,
            auto_check=False
        )

    def merged(self, *, name=None) -> Mesh:
        # No-op to be extended to symmetries
        return self.copy(name=name)

    def clipped(self, *, origin, normal, name=None) -> "Mesh":
        """
        Clip the mesh by a plane defined by `origin` and `normal`.

        Parameters
        ----------
        origin : np.ndarray
            The point in space where the clipping plane intersects (3D point).
        normal : np.ndarray
            The normal vector defining the orientation of the clipping plane.
        name: Optional[str]
            A name for the newly created mesh

        Returns
        -------
        Mesh
            A new Mesh instance that has been clipped.
        """
        new_vertices, new_faces, face_parent = \
                clip_faces(self.vertices, self._faces, normal, origin)
        new_metadata = {k: self.faces_metadata[k][face_parent] for k in self.faces_metadata}
        if name is None and self.name is not None:
            name = f"{self.name}_clipped"
        return Mesh(vertices=new_vertices, faces=new_faces, faces_metadata=new_metadata, name=name)


def to_new_mesh(old_mesh):
    # Temporary method for testing new method while the former implementation
    # is still the default
    return Mesh(old_mesh.vertices, old_mesh.faces, name=old_mesh.name)
