"""A set of meshes that can be used as a Mesh."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
from itertools import chain, accumulate
from functools import lru_cache
from typing import Iterable, Union

import numpy as np

from capytaine.meshes.geometry import Abstract3DObject, ClippableMixin, inplace_transformation
from capytaine.meshes.surface_integrals import SurfaceIntegralsMixin
from capytaine.meshes.meshes import Mesh

LOG = logging.getLogger(__name__)


class CollectionOfMeshes(ClippableMixin, SurfaceIntegralsMixin, Abstract3DObject):
    """A tuple of meshes.
    It gives access to all the vertices of all the sub-meshes as if it were a mesh itself.
    Collections can be nested to store meshes in a tree structure.

    Parameters
    ----------
    meshes: Iterable of Mesh or CollectionOfMeshes
        meshes in the collection
    name : str, optional
        a name for the collection
    """

    def __init__(self, meshes: Iterable[Union[Mesh, 'CollectionOfMeshes']], name=None):

        self._meshes = tuple(meshes)

        for mesh in self._meshes:
            assert isinstance(mesh, Mesh) or isinstance(mesh, CollectionOfMeshes)

        if name is None:
            self.name = f'collection_of_meshes_{next(Mesh._ids)}'
        else:
            self.name = str(name)

        LOG.debug(f"New collection of meshes: {repr(self)}")

    def __short_str__(self):
        return (f"{self.__class__.__name__}(..., name=\"{self.name}\")")

    def __str__(self):
        if len(self._meshes) < 3:
            meshes_str = ', '.join(m.__short_str__() for m in self._meshes)
        else:
            meshes_str = self._meshes[0].__short_str__() + ", ..., " + self._meshes[-1].__short_str__()
        return f"{self.__class__.__name__}([{meshes_str}], name=\"{self.name}\")"

    def __repr__(self):
        return f"{self.__class__.__name__}({', '.join(str(m) for m in self._meshes)}, name=\"{self.name}\")"

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    def __rich_repr__(self):
        class WrappedString:
            def __init__(self, s):
                self.s = s
            def __repr__(self):
                return self.s
        for m in self._meshes:
            yield m
        yield "name", self.name

    def __iter__(self):
        return iter(self._meshes)

    def __len__(self):
        return len(self._meshes)

    def __getitem__(self, item):
        return self._meshes.__getitem__(item)

    def __eq__(self, other):
        if isinstance(other, CollectionOfMeshes):
            return self._meshes == other._meshes
        else:
            return NotImplemented

    def __hash__(self):
        return hash(self._meshes)

    def tree_view(self, **kwargs):
        body_tree_views = []
        for i, mesh in enumerate(self):
            tree_view = mesh.tree_view(**kwargs)
            if i == len(self)-1:
                prefix = ' └─'
                shift  = '   '
            else:
                prefix = ' ├─'
                shift  = ' │ '
            body_tree_views.append(prefix + tree_view.replace('\n', '\n' + shift))

        return self.__short_str__() + '\n' + '\n'.join(body_tree_views)

    def path_to_leaf(self):
        """
        Builds a list of lists of paths from the collection corresponding to the
        root of the tree to the submeshes corresponding to the leaves
        """
        ptl = []
        for i, mesh in enumerate(self):
           for path in mesh.path_to_leaf():
               ptl.append([i] + path)
        return ptl

    def copy(self, name=None):
        from copy import deepcopy
        new_mesh = deepcopy(self)
        if name is not None:
            new_mesh.name = name
        return new_mesh

    @inplace_transformation
    def heal_mesh(self, closed_mesh=False):
        for mesh in self:
            mesh.heal_mesh(closed_mesh=closed_mesh)

    @inplace_transformation
    def with_normal_vector_going_down(self):
        for mesh in self:
            mesh.with_normal_vector_going_down()

    ##############
    # Properties #
    ##############

    @property
    def nb_submeshes(self):
        return len(self)

    @property
    def nb_vertices(self):
        return sum(mesh.nb_vertices for mesh in self)

    @property
    def nb_faces(self):
        return sum(mesh.nb_faces for mesh in self)

    @property
    def volume(self):
        return sum(mesh.volume for mesh in self)

    @property
    def vertices(self):
        return np.concatenate([mesh.vertices for mesh in self])

    @property
    def faces(self):
        """Return the indices of the vertices forming each of the faces. For the
        later submeshes, the indices of the vertices has to be shifted to
        correspond to their index in the concatenated array self.vertices.
        """
        nb_vertices = accumulate(chain([0], (mesh.nb_vertices for mesh in self[:-1])))
        return np.concatenate([mesh.faces + nbv for mesh, nbv in zip(self, nb_vertices)])

    @property
    def faces_normals(self):
        return np.concatenate([mesh.faces_normals for mesh in self])

    @property
    def faces_areas(self):
        return np.concatenate([mesh.faces_areas for mesh in self])

    @property
    def faces_centers(self):
        return np.concatenate([mesh.faces_centers for mesh in self])

    @property
    def faces_radiuses(self):
        return np.concatenate([mesh.faces_radiuses for mesh in self])

    @property
    def quadrature_points(self):
        quad_submeshes = [mesh.quadrature_points for mesh in self]
        return (
            np.concatenate([quad[0] for quad in quad_submeshes]),  # Points
            np.concatenate([quad[1] for quad in quad_submeshes])   # Weights
        )

    @property
    def quadrature_method(self):
        methods_submeshes = [mesh.quadrature_method for mesh in self]
        if len(set(methods_submeshes)) == 1:
            return methods_submeshes[0]  # All the same methods
        else:
            return "Mixed quadrature method"

    def compute_quadrature(self, method):
        for mesh in self:
            mesh.compute_quadrature(method)

    @property
    def center_of_mass_of_nodes(self):
        return sum([mesh.nb_vertices*mesh.center_of_mass_of_nodes for mesh in self])/self.nb_vertices

    @property
    @lru_cache(maxsize=1024)
    def diameter_of_nodes(self):
        return self.merged().diameter_of_nodes  # TODO: improve implementation

    def indices_of_mesh(self, mesh_index: int) -> slice:
        """Return the indices of the faces for the sub-mesh given as argument."""
        start = sum((mesh.nb_faces for mesh in self[:mesh_index]))  # Number of faces in previous meshes
        return slice(start, start + self[mesh_index].nb_faces)

    def submesh_containing_face(self, id_face):
        total_faces = 0
        for id_mesh in range(self.nb_submeshes):
            total_faces += self[id_mesh].nb_faces
            if id_face < total_faces:
                return id_mesh, id_face - (total_faces - self[id_mesh].nb_faces)

    ##################
    # Transformation #
    ##################

    def join_meshes(*meshes, name=None):
        return CollectionOfMeshes(meshes, name=name)

    def __add__(self, mesh_to_add):
        return self.join_meshes(mesh_to_add)

    def merged(self, name=None) -> Mesh:
        """Merge the sub-meshes and return a full mesh.
        If the collection contains other collections, they are merged recursively.
        Optionally, a new name can be given to the resulting mesh."""
        if name is None:
            name = self.name
        merged = Mesh(self.vertices, self.faces, name=name)
        merged.merge_duplicates()
        merged.heal_triangles()
        return merged

    def extract_one_face(self, id_face):
        id_mesh, relative_id_face = self.submesh_containing_face(id_face)
        mesh = self[id_mesh]

        extracted_mesh = mesh.extract_one_face(relative_id_face)

        if hasattr(mesh, '__internals__'):
            for prop in mesh.__internals__:
                if prop[:4] == "face":
                    extracted_mesh.__internals__[prop] = mesh.__internals__[prop][[relative_id_face]]

        return extracted_mesh

    def extract_faces(self, *args, **kwargs):
        return self.merged().extract_faces(*args, **kwargs)

    def sliced_by_plane(self, plane):
        return CollectionOfMeshes([mesh.sliced_by_plane(plane) for mesh in self], name=self.name)

    @inplace_transformation
    def translate(self, vector):
        for mesh in self:
            mesh.translate(vector)

    @inplace_transformation
    def rotate(self, axis, angle):
        for mesh in self:
            mesh.rotate(axis, angle)

    @inplace_transformation
    def mirror(self, plane):
        for mesh in self:
            mesh.mirror(plane)

    @inplace_transformation
    def clip(self, plane):
        self._clipping_data = {'faces_ids': []}
        faces_shifts = list(accumulate(chain([0], (mesh.nb_faces for mesh in self[:-1]))))
        for mesh, faces_shift in zip(self, faces_shifts):
            mesh.clip(plane)
            self._clipping_data['faces_ids'].extend([i + faces_shift for i in mesh._clipping_data['faces_ids']])
        self._clipping_data['faces_ids'] = np.asarray(self._clipping_data['faces_ids'])
        self.prune_empty_meshes()

    def symmetrized(self, plane):
        from capytaine.meshes.symmetric import ReflectionSymmetricMesh
        half = self.clipped(plane, name=f"{self.name}_half")
        return ReflectionSymmetricMesh(half, plane=plane, name=f"symmetrized_of_{self.name}")

    @inplace_transformation
    def prune_empty_meshes(self):
        """Remove empty meshes from the collection."""
        remaining_meshes = tuple(mesh for mesh in self if mesh.nb_faces > 0 and mesh.nb_vertices > 0)
        if len(remaining_meshes) == 0:
            self._meshes = (Mesh(name="empty_mesh"),)
        else:
            self._meshes = remaining_meshes

    @property
    def axis_aligned_bbox(self):
        """Get the axis aligned bounding box of the mesh.

        Returns
        -------
        tuple
            (xmin, xmax, ymin, ymax, zmin, zmax)
        """
        if self.nb_vertices > 0:
            x, y, z = self.vertices.T
            return (x.min(), x.max(),
                    y.min(), y.max(),
                    z.min(), z.max())
        else:
            return tuple(np.zeros(6))

    def show(self, **kwargs):
        from capytaine.ui.vtk.mesh_viewer import MeshViewer

        viewer = MeshViewer()
        for mesh in self:
            viewer.add_mesh(mesh.merged(), **kwargs)
        viewer.show()
        viewer.finalize()

    def show_matplotlib(self, *args, **kwargs):
        self.merged().show_matplotlib(*args, **kwargs)

    def lowest_lid_position(self, *args, **kwargs):
        return self.merged().lowest_lid_position(*args, **kwargs)

    def generate_lid(self, *args, **kwargs):
        return self.merged().generate_lid(*args, **kwargs)
