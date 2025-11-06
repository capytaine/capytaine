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
from typing import List

import numpy as np

from .geometry import (
    compute_faces_areas,
    compute_faces_centers,
    compute_faces_normals,
    compute_faces_radii,
    get_vertices_face,
)
from .clean import clean_mesh
from .quality import check_mesh_quality

LOG = logging.getLogger(__name__)


class Mesh:
    """Mesh class for representing and manipulating 3D surface meshes.

    Parameters
    ----------
    vertices : np.ndarray, optional
        Array of mesh vertices coordinates with shape (n_vertices, 3).
        Each row represents one vertex's (x, y, z) coordinates.
    faces : List[List[int]] or np.ndarray, optional
        Array of mesh connectivities for panels. Each row contains indices
        of vertices that form a face (triangles or quads).
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
        name: str = None,
        *,
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

        # Optional name
        self.name = str(name) if name is not None else None

        # Skip cleaning/quality if mesh is completely empty
        n_faces = len(self._faces)
        n_vertices = len(self.vertices)
        if not (n_vertices == 0 and n_faces == 0):
            if auto_clean:
                self._clean()
            if auto_check:
                self.check_quality()

    ## MAIN METRICS AND DISPLAY

    @property
    def nb_vertices(self) -> int:
        """Number of vertices in the mesh."""
        return len(self.vertices)

    @property
    def nb_faces(self) -> int:
        """Number of faces in the mesh."""
        return len(self._faces)

    @property
    def nb_triangles(self) -> int:
        """Number of triangular faces (3-vertex) in the mesh."""
        return sum(1 for f in self._faces if len(f) == 3)

    @property
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

    def __str__(self) -> str:
        return (f"Mesh(vertices=[[... {self.nb_vertices} vertices ...]], "
                + f"faces=[[... {self.nb_faces} faces ...]]"
                + ", name=\"{self.name}\")" if self.name is not None else ")")

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
    def from_list_of_faces(cls, list_faces, *, name=None) -> "Mesh":
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
        name: str, optional
            A name for the new mesh.

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

        return cls(vertices=np.array(unique_vertices), faces=indexed_faces, name=name)

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
        return list_faces

    ## INTERFACE FOR BEM SOLVER

    @property
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

    @property
    def faces_normals(self) -> np.ndarray:
        """Normal vectors for each face.

        Returns
        -------
        np.ndarray
            Array of shape (n_faces, 3) containing unit normal vectors.
        """
        return compute_faces_normals(self.vertices, self._faces)

    @property
    def faces_areas(self) -> np.ndarray:
        """Surface area of each face.

        Returns
        -------
        np.ndarray
            Array of shape (n_faces,) containing the area of each face.
        """
        return compute_faces_areas(self.vertices, self._faces)

    @property
    def faces_centers(self) -> np.ndarray:
        """Geometric centers of each face.

        Returns
        -------
        np.ndarray
            Array of shape (n_faces, 3) containing the center point of each face.
        """
        return compute_faces_centers(self.vertices, self._faces)

    @property
    def faces_radiuses(self) -> np.ndarray:
        """Radii of each face (circumradius or characteristic size).

        Returns
        -------
        np.ndarray
            Array of shape (n_faces,) containing the radius of each face.
        """
        return compute_faces_radii(self.vertices, self._faces)

    @property
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

    @property
    def quadrature_points(self) -> Tuple[np.ndarray, np.ndarray]:
        """Quadrature points and weights for numerical integration.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            (points, weights) where points has shape (n_faces, 1, 3) and
            weights has shape (n_faces, 1), using face centers and areas.
        """
        return (self.faces_centers.reshape((-1, 1, 3)), self.faces_areas.reshape(-1, 1))

    ## QUALITY CHECK

    def _clean(self, max_iter=5, tol=1e-8):
        """Clean the mesh by applying geometric and topological simplifications iteratively.

        Parameters
        ----------
        max_iter : int, optional
            Maximum number of iterations to perform. Defaults to 5.
        tol : float, optional
            Tolerance for merging vertices and removing small faces. Defaults to 1e-8.
        """
        self.vertices, self._faces = clean_mesh(
            self.vertices, self._faces, max_iter, tol
        )

    def check_quality(self):
        """Run geometric and metric quality checks on the mesh instance.

        Notes
        -----
        Logs warnings for potential mesh quality issues such as degenerate faces,
        non-manifold edges, or improper face orientations.
        """
        check_mesh_quality(self.vertices, self._faces)

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
        return Mesh.from_list_of_faces(selected_faces, name=name)

    def translated(self, dx=0.0, dy=0.0, dz=0.0, *, name=None) -> "Mesh":
        """Return a new Mesh translated along x, y, z axes.

        Parameters
        ----------
        dx : float, optional
            Translation along the x-axis. Defaults to 0.0.
        dy : float, optional
            Translation along the y-axis. Defaults to 0.0.
        dz : float, optional
            Translation along the z-axis. Defaults to 0.0.
        name: str, optional
            A name for the new translated object

        Returns
        -------
        Mesh
            New translated mesh instance.
        """
        shift = np.array([dx, dy, dz])
        return Mesh(vertices=self.vertices + shift, faces=self._faces, name=name)

    def rotated(
        self,
        angle: float = 0.0,
        axis: str = "z",
        angle_type: str = "rad",
        rotation_matrix: np.ndarray = None,
        *,
        name=None
    ) -> "Mesh":
        """Return a new Mesh rotated either by an angle around an axis or by a provided rotation matrix.

        Parameters
        ----------
        angle : float, optional
            Angle of rotation. Interpreted in radians by default, unless angle_type is 'deg'.
            Ignored if rotation_matrix is provided.
        axis : str, optional
            Axis of rotation: 'x', 'y', or 'z'. Ignored if rotation_matrix is provided.
        angle_type : str, optional
            Unit of the angle: 'rad' (default) or 'deg'.
        rotation_matrix : np.ndarray, optional
            Directly provide a 3x3 rotation matrix. If given, angle and axis are ignored.
        name: str, optional
            A name for the new translated object

        Returns
        -------
        Mesh
            New rotated mesh instance.
        """
        if rotation_matrix is not None:
            if rotation_matrix.shape != (3, 3):
                raise ValueError("rotation_matrix must be of shape (3, 3)")
            R = rotation_matrix
        else:
            # Convert angle if needed
            if angle_type == "deg":
                angle = np.deg2rad(angle)
            elif angle_type != "rad":
                raise ValueError("angle_type must be 'rad' or 'deg'.")

            c, s = np.cos(angle), np.sin(angle)
            if axis == "x":
                R = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
            elif axis == "y":
                R = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
            elif axis == "z":
                R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
            else:
                raise ValueError("Axis must be 'x', 'y', or 'z'.")

        # Apply rotation
        new_vertices = self.vertices @ R.T
        return Mesh(vertices=new_vertices, faces=self._faces, name=name)


    def join_meshes(self, other: "Mesh", *, name=None) -> "Mesh":
        """Join two meshes and return a new Mesh instance.

        Parameters
        ----------
        other : Mesh
            Another mesh to combine with this one.
        name: str, optional
            A name for the new object

        Returns
        -------
        Mesh
            New mesh containing vertices and faces from both meshes.

        See Also
        --------
        __add__ : Implements the + operator for mesh joining.
        """
        if not isinstance(other, Mesh):
            raise TypeError("Only Mesh instances can be added together.")

        # 1) Stack the vertices
        offset = len(self.vertices)
        new_vertices = np.vstack([self.vertices, other.vertices])

        # 2) Offset each face index in `other.faces` by `offset`
        #    and concatenate the two face‐lists
        new_faces = []
        # copy self._faces
        for face in self._faces:
            new_faces.append(list(face))
        # offset other.faces
        for face in other._faces:
            new_faces.append([idx + offset for idx in face])

        # 3) Return a new Mesh
        return Mesh(vertices=new_vertices, faces=new_faces, name=name)

    def __add__(self, other: "Mesh") -> "Mesh":
        """Combine two meshes using the + operator.

        Parameters
        ----------
        other : Mesh
            Another mesh to combine with this one.

        Returns
        -------
        Mesh
            New mesh containing vertices and faces from both meshes.

        Raises
        ------
        TypeError
            If other is not a Mesh instance.

        Notes
        -----
        Vertex indices in the second mesh are automatically offset to account
        for vertices from the first mesh.
        """
        if self.name is not None or other.name is not None:
            name = f"{self.name}+{other.name}"
        else:
            name = None
        return self.join_meshes(other, name=name)

    def with_normal_vector_going_down(self) -> "Mesh":
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
            return Mesh(vertices=self.vertices, faces=self.faces[:, ::-1], name=self.name)
        else:
            return self
