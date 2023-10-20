"""Special meshes with symmetries, useful to speed up the computations."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging
import reprlib
from typing import Union, Callable, Iterable

import numpy as np

from capytaine.meshes.meshes import Mesh
from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.meshes.geometry import Axis, Plane, Oz_axis, inplace_transformation

LOG = logging.getLogger(__name__)


class SymmetricMesh(CollectionOfMeshes):
    def __repr__(self):
        reprer = reprlib.Repr()
        reprer.maxstring = 90
        reprer.maxother = 90
        slice_name = reprer.repr(self._meshes[0])
        if self.name is not None:
            return f"{self.__class__.__name__}({slice_name}, name={self.name})"
        else:
            return f"{self.__class__.__name__}({slice_name})"


class ReflectionSymmetricMesh(SymmetricMesh):
    """A mesh with one vertical symmetry plane.

    Parameters
    ----------
    half : Mesh or CollectionOfMeshes
        a mesh describing half of the body
    plane : Plane
        the symmetry plane across which the half body is mirrored
    name :str, optional
        a name for the mesh
    """

    def __init__(self, half: Union[Mesh, CollectionOfMeshes], plane: Plane, name=None):
        assert isinstance(half, Mesh) or isinstance(half, CollectionOfMeshes)
        assert isinstance(plane, Plane)
        assert plane.normal[2] == 0, "Only vertical reflection planes are supported in ReflectionSymmetry classes."

        other_half = half.mirrored(plane, name=f"mirrored_of_{half.name}")

        if name is None:
            name = f"reflection_of_{half.name}"

        self.plane = plane.copy()

        super().__init__((half, other_half), name=name)

        if self.name is not None:
            LOG.debug(f"New mirror symmetric mesh: {self.name}.")
        else:
            LOG.debug(f"New mirror symmetric mesh.")

    def __str__(self):
        return f"{self.__class__.__name__}({self.half.__short_str__()}, plane={self.plane}, name=\"{self.name}\")"

    def __repr__(self):
        return f"{self.__class__.__name__}({self.half}, plane={self.plane}, name=\"{self.name}\")"

    def __rich_repr__(self):
        yield self.half
        yield "plane", self.plane
        yield "name", self.name

    @property
    def half(self):
        return self[0]

    def tree_view(self, fold_symmetry=True, **kwargs):
        if fold_symmetry:
            return (self.__short_str__() + '\n' + ' ├─' + self.half.tree_view().replace('\n', '\n │ ') + '\n'
                    + f" └─mirrored copy of the above {self.half.__short_str__()}")
        else:
            return CollectionOfMeshes.tree_view(self, **kwargs)

    def __deepcopy__(self, *args):
        return ReflectionSymmetricMesh(self.half.copy(), self.plane, name=self.name)

    def join_meshes(*meshes, name=None):
        assert all(isinstance(mesh, ReflectionSymmetricMesh) for mesh in meshes), \
            "Only meshes with the same symmetry can be joined together."
        assert all(meshes[0].plane == mesh.plane for mesh in meshes), \
            "Only reflection symmetric meshes with the same reflection plane can be joined together."
        half_mesh = CollectionOfMeshes([mesh.half for mesh in meshes], name=f"half_of_{name}" if name is not None else None)
        return ReflectionSymmetricMesh(half_mesh, plane=meshes[0].plane, name=name)

    @inplace_transformation
    def translate(self, vector):
        self.plane.translate(vector)
        CollectionOfMeshes.translate(self, vector)
        return self

    @inplace_transformation
    def rotate(self, axis: Axis, angle: float):
        self.plane.rotate(axis, angle)
        CollectionOfMeshes.rotate(self, axis, angle)
        return self

    @inplace_transformation
    def mirror(self, plane: Plane):
        self.plane.mirror(plane)
        CollectionOfMeshes.mirror(self, plane)
        return self


class TranslationalSymmetricMesh(SymmetricMesh):
    """A mesh with a repeating pattern by translation.

    Parameters
    ----------
    mesh_slice : Mesh or CollectionOfMeshes
        the pattern that will be repeated to form the whole body
    translation : array(3)
        the vector of the translation
    nb_repetitions : int, optional
        the number of repetitions of the pattern (excluding the original one, default: 1)
    name : str, optional
        a name for the mesh
    """

    def __init__(self, mesh_slice: Union[Mesh, CollectionOfMeshes], translation, nb_repetitions=1, name=None):
        assert isinstance(mesh_slice, Mesh) or isinstance(mesh_slice, CollectionOfMeshes)
        assert isinstance(nb_repetitions, int)
        assert nb_repetitions >= 1

        translation = np.asarray(translation).copy()
        assert translation.shape == (3,)
        assert translation[2] == 0  # Only horizontal translation are supported.

        slices = [mesh_slice]
        for i in range(1, nb_repetitions+1):
            slices.append(mesh_slice.translated(vector=i*translation, name=f"repetition_{i}_of_{mesh_slice.name}"))

        if name is None:
            name = f"translation_of_{mesh_slice.name}"

        self.translation = translation

        super().__init__(slices, name=name)

        if self.name is not None:
            LOG.debug(f"New translation symmetric mesh: {self.name}.")
        else:
            LOG.debug(f"New translation symmetric mesh.")

    @property
    def first_slice(self):
        return self[0]

    def __str__(self):
        return f"{self.__class__.__name__}({self.first_slice.__short_str__()}, translation={self.translation}, nb_repetitions={len(self)-1}, name=\"{self.name}\")"

    def __repr__(self):
        return f"{self.__class__.__name__}({self.first_slice}, translation={self.translation}, nb_repetitions={len(self)-1}, name=\"{self.name}\")"

    def __rich_repr__(self):
        yield self.first_slice
        yield "translation", self.translation
        yield "nb_repetitions", len(self)-1
        yield "name", self.name

    def tree_view(self, fold_symmetry=True, **kwargs):
        if fold_symmetry:
            return (self.__short_str__() + '\n' + ' ├─' + self.first_slice.tree_view().replace('\n', '\n │ ') + '\n'
                    + f" └─{len(self)-1} translated copies of the above {self.first_slice.__short_str__()}")
        else:
            return CollectionOfMeshes.tree_view(self, **kwargs)

    def __deepcopy__(self, *args):
        return TranslationalSymmetricMesh(self.first_slice.copy(), self.translation, nb_repetitions=len(self) - 1, name=self.name)

    @inplace_transformation
    def translate(self, vector):
        CollectionOfMeshes.translate(self, vector)
        return self

    @inplace_transformation
    def rotate(self, axis: Axis, angle: float):
        self.translation = axis.rotate_vectors([self.translation], angle)[0, :]
        CollectionOfMeshes.rotate(self, axis, angle)
        return self

    @inplace_transformation
    def mirror(self, plane: Plane):
        self.translation -= 2 * (self.translation @ plane.normal) * plane.normal
        CollectionOfMeshes.mirror(self, plane)
        return self

    def join_meshes(*meshes, name=None):
        assert all(isinstance(mesh, TranslationalSymmetricMesh) for mesh in meshes), \
            "Only meshes with the same symmetry can be joined together."
        assert all(np.allclose(meshes[0].translation, mesh.translation) for mesh in meshes), \
            "Only translation symmetric meshes with the same translation vector can be joined together."
        assert all(len(meshes[0]) == len(mesh) for mesh in meshes), \
            "Only symmetric meshes with the same number of elements can be joined together."
        mesh_strip = CollectionOfMeshes([mesh.first_slice for mesh in meshes], name=f"strip_of_{name}" if name is not None else None)
        return TranslationalSymmetricMesh(mesh_strip, translation=meshes[0].translation, nb_repetitions=len(meshes[0]) - 1, name=name)


def build_regular_array_of_meshes(base_mesh, distance, nb_bodies):
    """Create an array of objects using TranslationalSymmetries.

    Parameters
    ----------
    base_mesh : Mesh or CollectionOfMeshes or SymmetricMesh
        The mesh to duplicate to create the array
    distance : float
        Center-to-center distance between objects in the array
    nb_bodies : couple of ints
        Number of objects in the x and y directions.

    Returns
    -------
    TranslationalSymmetricMesh
    """
    if nb_bodies[0] == 1:
        line = base_mesh
    else:
        line = TranslationalSymmetricMesh(base_mesh, translation=(distance, 0.0, 0.0), nb_repetitions=nb_bodies[0] - 1,
                                          name=f'line_of_{base_mesh.name}')
    if nb_bodies[1] == 1:
        array = line
    else:
        array = TranslationalSymmetricMesh(line, translation=(0.0, distance, 0.0), nb_repetitions=nb_bodies[1] - 1,
                                           name=f'array_of_{base_mesh.name}')
    return array


class AxialSymmetricMesh(SymmetricMesh):
    """A mesh with a repeating pattern by rotation.

    Parameters
    ----------
    mesh_slice : Mesh or CollectionOfMeshes
        the pattern that will be repeated to form the whole body
    axis : Axis, optional
        symmetry axis
    nb_repetitions : int, optional
        the number of repetitions of the pattern (excluding the original one, default: 1)
    name : str, optional
        a name for the mesh
    """
    def __init__(self, mesh_slice: Union[Mesh, CollectionOfMeshes], axis: Axis=Oz_axis, nb_repetitions: int=1, name=None):
        assert isinstance(mesh_slice, Mesh) or isinstance(mesh_slice, CollectionOfMeshes)
        assert isinstance(nb_repetitions, int)
        assert nb_repetitions >= 1
        assert isinstance(axis, Axis)

        slices = [mesh_slice]
        for i in range(1, nb_repetitions+1):
            slices.append(mesh_slice.rotated(axis, angle=2*i*np.pi/(nb_repetitions+1),
                                             name=f"rotation_{i}_of_{mesh_slice.name}"))

        if name is None:
            name = f"rotation_of_{mesh_slice.name}"

        self.axis = axis.copy()

        super().__init__(slices, name=name)

        if not axis.is_parallel_to(Oz_axis):
            LOG.warning(f"{self.name} is an axi-symmetric mesh along a non vertical axis.")

        if self.name is not None:
            LOG.debug(f"New rotation symmetric mesh: {self.name}.")
        else:
            LOG.debug(f"New rotation symmetric mesh.")

    @staticmethod
    def from_profile(profile: Union[Callable, Iterable[float]],
                     z_range: Iterable[float]=np.linspace(-5, 0, 20),
                     axis: Axis=Oz_axis,
                     nphi: int=20,
                     name=None):
        """Return a floating body using the axial symmetry.
        The shape of the body can be defined either with a function defining the profile as [f(z), 0, z] for z in z_range.
        Alternatively, the profile can be defined as a list of points.
        The number of vertices along the vertical direction is len(z_range) in the first case and profile.shape[0] in the second case.

        Parameters
        ----------
        profile : function(float → float)  or  array(N, 3)
            define the shape of the body either as a function or a list of points.
        z_range: array(N), optional
            used only if the profile is defined as a function.
        axis : Axis
            symmetry axis
        nphi : int, optional
            number of vertical slices forming the body
        name : str, optional
            name of the generated body (optional)

        Returns
        -------
        AxialSymmetricMesh
            the generated mesh
        """

        if name is None:
            name = "axisymmetric_mesh"

        if callable(profile):
            z_range = np.asarray(z_range)
            x_values = [profile(z) for z in z_range]
            profile_array = np.stack([x_values, np.zeros(len(z_range)), z_range]).T
        else:
            profile_array = np.asarray(profile)

        assert len(profile_array.shape) == 2
        assert profile_array.shape[1] == 3

        n = profile_array.shape[0]
        angle = 2 * np.pi / nphi

        nodes_slice = np.concatenate([profile_array, axis.rotate_points(profile_array, angle)])
        faces_slice = np.array([[i, i+n, i+n+1, i+1] for i in range(n-1)])
        body_slice = Mesh(nodes_slice, faces_slice, name=f"slice_of_{name}")
        body_slice.merge_duplicates()
        body_slice.heal_triangles()

        return AxialSymmetricMesh(body_slice, axis=axis, nb_repetitions=nphi - 1, name=name)

    @property
    def first_slice(self):
        return self[0]

    def __str__(self):
        return f"{self.__class__.__name__}({self.first_slice.__short_str__()}, axis={self.axis}, nb_repetitions={len(self)-1}, name=\"{self.name}\")"

    def __repr__(self):
        return f"{self.__class__.__name__}({self.first_slice}, axis={self.axis}, nb_repetitions={len(self)-1}, name=\"{self.name}\")"

    def __rich_repr__(self):
        yield self.first_slice
        yield "axis", self.axis
        yield "nb_repetitions", len(self)-1
        yield "name", self.name

    def tree_view(self, fold_symmetry=True, **kwargs):
        if fold_symmetry:
            return (self.__short_str__() + '\n' + ' ├─' + self.first_slice.tree_view().replace('\n', '\n │ ') + '\n'
                    + f" └─{len(self)-1} rotated copies of the above {self.first_slice.__short_str__()}")
        else:
            return CollectionOfMeshes.tree_view(self, **kwargs)

    def __deepcopy__(self, *args):
        return AxialSymmetricMesh(self.first_slice.copy(), axis=self.axis.copy(), nb_repetitions=len(self) - 1, name=self.name)

    def join_meshes(*meshes, name=None):
        assert all(isinstance(mesh, AxialSymmetricMesh) for mesh in meshes), \
            "Only meshes with the same symmetry can be joined together."
        assert all(meshes[0].axis == mesh.axis for mesh in meshes), \
            "Only axisymmetric meshes with the same symmetry axis can be joined together."
        assert all(len(meshes[0]) == len(mesh) for mesh in meshes), \
            "Only axisymmetric meshes with the same number of elements can be joined together."
        mesh_slice = CollectionOfMeshes([mesh.first_slice for mesh in meshes], name=f"slice_of_{name}" if name is not None else None)
        return AxialSymmetricMesh(mesh_slice, axis=meshes[0].axis, nb_repetitions=len(meshes[0]) - 1, name=name)

    @inplace_transformation
    def translate(self, vector):
        self.axis.translate(vector)
        CollectionOfMeshes.translate(self, vector)
        return self

    @inplace_transformation
    def rotate(self, other_axis: Axis, angle: float):
        self.axis.rotate(other_axis, angle)
        CollectionOfMeshes.rotate(self, other_axis, angle)
        return self

    @inplace_transformation
    def mirror(self, plane: Plane):
        self.axis.mirror(plane)
        CollectionOfMeshes.mirror(self, plane)
        return self
