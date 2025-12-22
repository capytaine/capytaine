"""Floating bodies to be used in radiation-diffraction problems."""
# Copyright (C) 2017-2024 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import logging
import copy
from itertools import chain, accumulate
from functools import lru_cache

import numpy as np
import xarray as xr

from capytaine.meshes.mesh_like_protocol import MeshLike
from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.meshes.geometry import Abstract3DObject, ClippableMixin, inplace_transformation
from capytaine.meshes.properties import connected_components, connected_components_of_waterline
from capytaine.meshes.meshes import Mesh
from capytaine.meshes.symmetric import build_regular_array_of_meshes
from capytaine.new_meshes.abstract_meshes import AbstractMesh
from capytaine.bodies.dofs import RigidBodyDofsPlaceholder, TRANSLATION_DOFS_DIRECTIONS, ROTATION_DOFS_AXIS
from capytaine.bodies.hydrostatics import _HydrostaticsMixin

from capytaine.tools.optional_imports import silently_import_optional_dependency
meshio = silently_import_optional_dependency("meshio")

LOG = logging.getLogger(__name__)


class FloatingBody(_HydrostaticsMixin, ClippableMixin, Abstract3DObject):
    """A floating body described as a mesh and some degrees of freedom.

    The mesh structure is stored as a Mesh from capytaine.mesh.mesh or a
    CollectionOfMeshes from capytaine.mesh.meshes_collection.

    The degrees of freedom (dofs) are stored as a dict associating a name to
    a complex-valued array of shape (nb_faces, 3). To each face of the body
    (as indexed in the mesh) corresponds a complex-valued 3d vector, which
    defines the displacement of the center of the face in frequency domain.

    Parameters
    ----------
    mesh : MeshLike, optional
        the mesh describing the geometry of the hull of the floating body.
        If none is given, a empty one is created.
    lid_mesh : MeshLike or None, optional
        a mesh of an internal lid for irregular frequencies removal.
        Unlike the mesh of the hull, no dof is defined on the lid_mesh.
        If none is given, none is used when solving the Boundary Integral Equation.
    dofs : dict, optional
        the degrees of freedom of the body.
        If none is given, a empty dictionary is initialized.
    mass : float or None, optional
        the mass of the body in kilograms.
        Required only for some hydrostatics computation.
        If None, the mass is implicitly assumed to be the mass of displaced water.
    center_of_mass: 3-element array, optional
        the position of the center of mass.
        Required only for some hydrostatics computation.
    name : str, optional
        a name for the body.
        If none is given, the one of the mesh is used.
    """

    def __init__(self, mesh=None, dofs=None, mass=None, center_of_mass=None, name=None, *, lid_mesh=None):
        if mesh is None:
            self.mesh = Mesh(name="dummy_mesh")

        elif meshio is not None and isinstance(mesh, meshio._mesh.Mesh):
            from capytaine.io.meshio import load_from_meshio
            self.mesh = load_from_meshio(mesh)

        elif isinstance(mesh, MeshLike):
            self.mesh = mesh

        else:
            raise TypeError("Unrecognized `mesh` object passed to the FloatingBody constructor.")

        if lid_mesh is not None:
            if lid_mesh.nb_faces == 0:
                LOG.warning("Lid mesh %s provided for body initialization is empty. The lid mesh is ignored.", lid_mesh)
                self.lid_mesh = None
            else:
                self.lid_mesh = lid_mesh.with_normal_vector_going_down(inplace=False)
        else:
            self.lid_mesh = None

        if name is None and mesh is None:
            self.name = "dummy_body"
        elif name is None:
            if hasattr(self.mesh, "name") and self.mesh.name is not None:
                self.name = self.mesh.name
            else:
                self.name = "anonymous_body"
        else:
            self.name = name

        self.mass = mass
        if center_of_mass is not None:
            self.center_of_mass = np.asarray(center_of_mass, dtype=float)
        else:
            self.center_of_mass = None

        if hasattr(self.mesh, "heal_mesh") and self.mesh.nb_vertices > 0 and self.mesh.nb_faces > 0:
            self.mesh.heal_mesh()

        if dofs is None:
            self.dofs = {}
        elif isinstance(dofs, RigidBodyDofsPlaceholder):
            if dofs.rotation_center is not None:
                self.rotation_center = np.asarray(dofs.rotation_center, dtype=float)
            self.dofs = {}
            self.add_all_rigid_body_dofs()
        else:
            self.dofs = dofs

        self._evaluate_full_mesh()

        LOG.info(f"New floating body: {self.__str__()}.")

        self._check_dofs_shape_consistency()

    def _evaluate_full_mesh(self):
        """Merge the mesh and lid_mesh, while keeping track of where each panel came from."""
        if self.lid_mesh is None:
            self.mesh_including_lid = self.mesh
            self.hull_mask = np.full((self.mesh.nb_faces,), True)
        else:
            self.mesh_including_lid, masks = self.mesh.join_meshes(self.lid_mesh, return_masks=True)
            self.hull_mask = masks[0]

    @staticmethod
    def from_meshio(mesh, name=None) -> 'FloatingBody':
        """Create a FloatingBody from a meshio mesh object.
        Kinda deprecated, use cpt.load_mesh instead."""
        from capytaine.io.meshio import load_from_meshio
        return FloatingBody(mesh=load_from_meshio(mesh, name), name=name)

    @staticmethod
    def from_file(filename: str, file_format=None, name=None) -> 'FloatingBody':
        """Create a FloatingBody from a mesh file using meshmagick.
        Kinda deprecated, use cpt.load_mesh instead."""
        from capytaine.io.mesh_loaders import load_mesh
        if name is None: name = filename
        mesh = load_mesh(filename, file_format, name=f"{name}_mesh")
        return FloatingBody(mesh, name=name)

    def __lt__(self, other: 'FloatingBody') -> bool:
        """Arbitrary order. The point is to sort together the problems involving the same body."""
        return self.name < other.name

    ##########
    #  Dofs  #
    ##########

    @property
    def nb_dofs(self) -> int:
        """Number of degrees of freedom."""
        return len(self.dofs)

    def add_translation_dof(self, direction=None, name=None, amplitude=1.0) -> None:
        """Add a new translation dof (in place).
        If no direction is given, the code tries to infer it from the name.

        Parameters
        ----------
        direction : array of shape (3,), optional
            the direction of the translation
        name : str, optional
            a name for the degree of freedom
        amplitude : float, optional
            amplitude of the dof (default: 1.0 m/s)
        """
        if direction is None:
            if name is not None and name.lower() in TRANSLATION_DOFS_DIRECTIONS:
                direction = TRANSLATION_DOFS_DIRECTIONS[name.lower()]
            else:
                raise ValueError("A direction needs to be specified for the dof.")

        if name is None:
            name = f"dof_{self.nb_dofs}_translation"

        direction = np.asarray(direction)
        assert direction.shape == (3,)

        motion = np.empty((self.mesh.nb_faces, 3))
        motion[:, :] = direction
        self.dofs[name] = amplitude * motion

    def add_rotation_dof(self, axis=None, name=None, amplitude=1.0) -> None:
        """Add a new rotation dof (in place).
        If no axis is given, the code tries to infer it from the name.

        Parameters
        ----------
        axis: Axis, optional
            the axis of the rotation
        name : str, optional
            a name for the degree of freedom
        amplitude : float, optional
            amplitude of the dof (default: 1.0)
        """
        if axis is None:
            if name is not None and name.lower() in ROTATION_DOFS_AXIS:
                axis_direction = ROTATION_DOFS_AXIS[name.lower()]
                for point_attr in ('rotation_center', 'center_of_mass', 'geometric_center'):
                    if hasattr(self, point_attr) and getattr(self, point_attr) is not None:
                        axis_point = getattr(self, point_attr)
                        LOG.info(f"The rotation dof {name} has been initialized around the point: "
                                 f"{self.__short_str__()}.{point_attr} = {getattr(self, point_attr)}")
                        break
                else:
                    axis_point = np.array([0, 0, 0])
                    LOG.warning(f"The rotation dof {name} has been initialized "
                                f"around the origin of the domain (0, 0, 0).")
            else:
                raise ValueError("A direction needs to be specified for the dof.")
        else:
            axis_point = axis.point
            axis_direction = axis.vector

        if name is None:
            name = f"dof_{self.nb_dofs}_rotation"

        if self.mesh.nb_faces == 0:
            self.dofs[name] = np.empty((self.mesh.nb_faces, 3))
        else:
            motion = np.cross(axis_point - self.mesh.faces_centers, axis_direction)
            self.dofs[name] = amplitude * motion

    def add_all_rigid_body_dofs(self) -> None:
        """Add the six degrees of freedom of rigid bodies (in place)."""
        self.add_translation_dof(name="Surge")
        self.add_translation_dof(name="Sway")
        self.add_translation_dof(name="Heave")
        self.add_rotation_dof(name="Roll")
        self.add_rotation_dof(name="Pitch")
        self.add_rotation_dof(name="Yaw")

    def integrate_pressure(self, pressure):
        forces = {}
        for dof_name in self.dofs:
            # Scalar product on each face:
            normal_dof_amplitude_on_face = - np.sum(self.dofs[dof_name] * self.mesh.faces_normals, axis=1)
            # The minus sign in the above line is because we want the force of the fluid on the body and not the force of the body on the fluid.
            # Sum over all faces:
            forces[dof_name] = np.sum(pressure * normal_dof_amplitude_on_face * self.mesh.faces_areas)
        return forces

    @inplace_transformation
    def keep_only_dofs(self, dofs):
        for dof in list(self.dofs.keys()):
            if dof not in dofs:
                del self.dofs[dof]

        if hasattr(self, 'inertia_matrix'):
            self.inertia_matrix = self.inertia_matrix.sel(radiating_dof=dofs, influenced_dof=dofs)
        if hasattr(self, 'hydrostatic_stiffness'):
            self.hydrostatic_stiffness = self.hydrostatic_stiffness.sel(radiating_dof=dofs, influenced_dof=dofs)

        return self

    def add_dofs_labels_to_vector(self, vector):
        """Helper function turning a bare vector into a vector labelled by the name of the dofs of the body,
        to be used for instance for the computation of RAO."""
        return xr.DataArray(data=np.asarray(vector), dims=['influenced_dof'],
                            coords={'influenced_dof': list(self.dofs)},
                            )

    def add_dofs_labels_to_matrix(self, matrix):
        """Helper function turning a bare matrix into a matrix labelled by the name of the dofs of the body,
        to be used for instance for the computation of RAO."""
        return xr.DataArray(data=np.asarray(matrix), dims=['influenced_dof', 'radiating_dof'],
                            coords={'influenced_dof': list(self.dofs), 'radiating_dof': list(self.dofs)},
                            )

    def _check_dofs_shape_consistency(self):
        for dof_name, dof in self.dofs.items():
            if np.array(dof).shape != (self.mesh.nb_faces, 3):
                raise ValueError(f"The array defining the dof {dof_name} of body {self.name} does not have the expected shape.\n"
                                 f"Expected shape: ({self.mesh.nb_faces}, 3)\n"
                                 f"  Actual shape: {dof.shape}")



    ###################
    # Transformations #
    ###################

    def __add__(self, body_to_add: 'FloatingBody') -> 'FloatingBody':
        return self.join_bodies(body_to_add)

    def join_bodies(*bodies, name=None) -> 'FloatingBody':
        if name is None:
            name = "+".join(body.name for body in bodies)
        meshes = [body.mesh.copy() for body in bodies if body.mesh.nb_faces > 0]
        if len(meshes) > 0:
            joined_mesh = meshes[0].join_meshes(*meshes[1:], name=f"{name}_mesh")
        else:
            joined_mesh = Mesh()

        if all(body.lid_mesh is None for body in bodies):
            joined_lid = None
        else:
            lid_meshes = [body.lid_mesh.copy() for body in bodies if body.lid_mesh is not None]
            joined_lid = lid_meshes[0].join_meshes(*lid_meshes[1:], name=f"{name}_lid_mesh")

        dofs = FloatingBody.combine_dofs(bodies)

        if all(body.mass is not None for body in bodies):
            new_mass = sum(body.mass for body in bodies)
        else:
            new_mass = None

        if (all(body.mass is not None for body in bodies)
                and all(body.center_of_mass is not None for body in bodies)):
            new_cog = sum(body.mass*np.asarray(body.center_of_mass) for body in bodies)/new_mass
        else:
            new_cog = None

        joined_bodies = FloatingBody(
            mesh=joined_mesh, lid_mesh=joined_lid, dofs=dofs,
            mass=new_mass, center_of_mass=new_cog, name=name
            )

        for matrix_name in ["inertia_matrix", "hydrostatic_stiffness"]:
            if all(hasattr(body, matrix_name) for body in bodies):
                from scipy.linalg import block_diag
                setattr(joined_bodies, matrix_name, joined_bodies.add_dofs_labels_to_matrix(
                        block_diag(*[getattr(body, matrix_name) for body in bodies])
                        ))

        return joined_bodies

    @staticmethod
    def combine_dofs(bodies) -> dict:
        """Combine the degrees of freedom of several bodies."""
        for body in bodies:
            body._check_dofs_shape_consistency()
        dofs = {}
        cum_nb_faces = accumulate(chain([0], (body.mesh.nb_faces for body in bodies)))
        total_nb_faces = sum(body.mesh.nb_faces for body in bodies)
        for body, nbf in zip(bodies, cum_nb_faces):
            # nbf is the cumulative number of faces of the previous subbodies,
            # that is the offset of the indices of the faces of the current body.
            for name, dof in body.dofs.items():
                new_dof = np.zeros((total_nb_faces, 3))
                new_dof[nbf:nbf+len(dof), :] = dof
                if '__' not in name:
                    new_dof_name = '__'.join([body.name, name])
                else:
                    # The body is probably a combination of bodies already.
                    # So for the associativity of the + operation,
                    # it is better to keep the same name.
                    new_dof_name = name
                dofs[new_dof_name] = new_dof
        return dofs

    def copy(self, name=None) -> 'FloatingBody':
        """Return a deep copy of the body.

        Parameters
        ----------
        name : str, optional
            a name for the new copy
        """
        self._check_dofs_shape_consistency()

        new_body = copy.deepcopy(self)
        if name is None:
            new_body.name = f"copy_of_{self.name}"
            LOG.debug(f"Copy {self.name}.")
        else:
            new_body.name = name
            LOG.debug(f"Copy {self.name} under the name {name}.")
        return new_body

    def assemble_regular_array(self, distance, nb_bodies):
        """Create an regular array of identical bodies.

        Parameters
        ----------
        distance : float
            Center-to-center distance between objects in the array
        nb_bodies : couple of ints
            Number of objects in the x and y directions.

        Returns
        -------
        FloatingBody
        """
        bodies = (self.translated((i*distance, j*distance, 0), name=f"{i}_{j}") for j in range(nb_bodies[1]) for i in range(nb_bodies[0]))
        array = FloatingBody.join_bodies(*bodies)
        array.mesh = build_regular_array_of_meshes(self.mesh, distance, nb_bodies)
        array.name = f"array_of_{self.name}"
        return array

    def assemble_arbitrary_array(self, locations:np.ndarray):

        if not isinstance(locations, np.ndarray):
            raise TypeError('locations must be of type np.ndarray')
        assert locations.shape[1] == 2, 'locations must be of shape nx2, received {:}'.format(locations.shape)

        fb_list = []
        for idx, li in enumerate(locations):
            fb1 = self.copy()
            fb1.translate(np.append(li,0))
            fb1.name = 'arbitrary_array_body{:02d}'.format(idx)
            fb_list.append(fb1)

        arbitrary_array = fb_list[0].join_bodies(*fb_list[1:])

        return arbitrary_array

    def extract_faces(self, id_faces_to_extract, return_index=False):
        """Create a new FloatingBody by extracting some faces from the mesh.
        The dofs evolve accordingly.
        The lid_mesh, center_of_mass, mass and hydrostatics data are discarded.
        """
        if isinstance(self.mesh, CollectionOfMeshes):
            raise NotImplementedError  # TODO

        if return_index:
            new_mesh, id_v = self.mesh.extract_faces(id_faces_to_extract, return_index)
        else:
            new_mesh = self.mesh.extract_faces(id_faces_to_extract, return_index)
        new_body = FloatingBody(new_mesh)
        LOG.info(f"Extract floating body from {self.name}.")

        new_body.dofs = {}
        for name, dof in self.dofs.items():
            new_body.dofs[name] = dof[id_faces_to_extract, :]

        if return_index:
            return new_body, id_v
        else:
            return new_body

    @inplace_transformation
    def mirror(self, plane):
        self.mesh.mirror(plane)
        if self.lid_mesh is not None:
            self.lid_mesh.mirror(plane)
        self._evaluate_full_mesh()
        for dof in self.dofs:
            self.dofs[dof] -= 2 * np.outer(np.dot(self.dofs[dof], plane.normal), plane.normal)
        for point_attr in ('geometric_center', 'rotation_center', 'center_of_mass'):
            if point_attr in self.__dict__ and self.__dict__[point_attr] is not None:
                point = np.array(self.__dict__[point_attr])
                shift = - 2 * (np.dot(point, plane.normal) - plane.c) * plane.normal
                self.__dict__[point_attr] = point + shift
        return self

    @inplace_transformation
    def translate(self, vector, *args, **kwargs):
        self.mesh.translate(vector, *args, **kwargs)
        if self.lid_mesh is not None:
            self.lid_mesh.translate(vector, *args, **kwargs)
        self._evaluate_full_mesh()
        for point_attr in ('geometric_center', 'rotation_center', 'center_of_mass'):
            if point_attr in self.__dict__ and self.__dict__[point_attr] is not None:
                self.__dict__[point_attr] = np.array(self.__dict__[point_attr]) + vector
        return self

    @inplace_transformation
    def rotate(self, axis, angle):
        self.mesh.rotate(axis, angle)
        if self.lid_mesh is not None:
            self.lid_mesh.rotate(axis, angle)
        self._evaluate_full_mesh()
        for point_attr in ('geometric_center', 'rotation_center', 'center_of_mass'):
            if point_attr in self.__dict__ and self.__dict__[point_attr] is not None:
                self.__dict__[point_attr] = axis.rotate_points([self.__dict__[point_attr]], angle)[0, :]
        for dof in self.dofs:
            self.dofs[dof] = axis.rotate_vectors(self.dofs[dof], angle)
        return self

    @inplace_transformation
    def clip(self, plane):
        self._check_dofs_shape_consistency()

        # Clip mesh
        LOG.info(f"Clipping {self.name} with respect to {plane}")
        self.mesh.clip(plane)
        if self.lid_mesh is not None:
            self.lid_mesh.clip(plane)
            if self.lid_mesh.nb_faces == 0:
                LOG.warning("Lid mesh %s is empty after clipping. The lid mesh is removed.", self.lid_mesh)
                self.lid_mesh = None
        self._evaluate_full_mesh()

        # Clip dofs
        ids = self.mesh._clipping_data['faces_ids']
        for dof in self.dofs:
            if len(ids) > 0:
                self.dofs[dof] = np.array(self.dofs[dof])[ids]
            else:
                self.dofs[dof] = np.empty((0, 3))
        return self

    def _apply_on_mesh(self, func, args, kwargs):
        mesh_with_dofs = self.mesh.with_metadata(**self.dofs)
        transformed_mesh = func(mesh_with_dofs, *args, **kwargs)
        new_dofs = {k: transformed_mesh.faces_metadata[k] for k in self.dofs}
        transformed_mesh = transformed_mesh.without_metadata(*self.dofs.keys())

        if self.lid_mesh is not None:
            transformed_lid_mesh = func(self.lid_mesh, *args, **kwargs)
            if transformed_lid_mesh.nb_faces == 0:
                LOG.warning("Empty lid mesh %s has been removed.", self.lid_mesh)
                transformed_lid_mesh = None
        else:
            transformed_lid_mesh = None

        return transformed_mesh, transformed_lid_mesh, new_dofs

    def clipped(self, *, origin, normal, name=None) -> "FloatingBody":
        if not isinstance(self.mesh, NewMesh):
            # Legacy path, to be removed
            return self.clip(
                    cpt.Plane(point=origin, normal=normal),
                    inplace=False,
                    name=name
                    )

        clipped_mesh, clipped_lid_mesh, updated_dofs = self._apply_on_mesh(
            self.mesh.__class__.clipped,
            (),
            {'origin': origin, 'normal': normal}
        )
        if name is None:
            name = "clipped_" + self.name
        return FloatingBody(
            mesh=clipped_mesh,
            lid_mesh=clipped_lid_mesh,
            dofs=updated_dofs,
            name=name
        )

    def immersed_part(self, free_surface=0.0, *, sea_bottom=None, water_depth=None, name=None) -> "FloatingBody":
        if not isinstance(self.mesh, AbstractMesh):
            # Legacy path, to be removed
            return self.keep_immersed_part(
                    free_surface,
                    inplace=False,
                    sea_bottom=sea_bottom,
                    water_depth=water_depth,
                    name=self.name if name is None else name,
                    )

        clipped_mesh, clipped_lid_mesh, updated_dofs = self._apply_on_mesh(
            self.mesh.__class__.immersed_part,
            (free_surface,),
            {'sea_bottom': sea_bottom, 'water_depth': water_depth}
        )
        if name is None:
            name = "immersed_part_of_" + self.name
        return FloatingBody(
            mesh=clipped_mesh,
            lid_mesh=clipped_lid_mesh,
            dofs=updated_dofs,
            name=name
        )

    #############
    #  Display  #
    #############

    def __short_str__(self):
        return (f"{self.__class__.__name__}(..., name=\"{self.name}\")")

    def _optional_params_str(self):
        items = []
        if self.mass is not None: items.append(f"mass={self.mass}, ")
        if self.center_of_mass is not None: items.append(f"center_of_mass={self.center_of_mass}, ")
        return ''.join(items)

    def __str__(self):
        short_dofs = '{' + ', '.join('"{}": ...'.format(d) for d in self.dofs) + '}'

        if self.lid_mesh is not None:
            lid_mesh_str = self.lid_mesh.__short_str__()
        else:
            lid_mesh_str = str(None)

        return (f"{self.__class__.__name__}(mesh={self.mesh.__short_str__()}, lid_mesh={lid_mesh_str}, "
                f"dofs={short_dofs}, {self._optional_params_str()}name=\"{self.name}\")")

    def __repr__(self):
        short_dofs = '{' + ', '.join('"{}": ...'.format(d) for d in self.dofs) + '}'

        if self.lid_mesh is not None:
            lid_mesh_str = str(self.lid_mesh)
        else:
            lid_mesh_str = str(None)

        return (f"{self.__class__.__name__}(mesh={str(self.mesh)}, lid_mesh={lid_mesh_str}, "
                f"dofs={short_dofs}, {self._optional_params_str()}name=\"{self.name}\")")

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__())

    def __rich_repr__(self):
        class DofWithShortRepr:
            def __repr__(self):
                return '...'
        yield "mesh", self.mesh
        yield "lid_mesh", self.lid_mesh
        yield "dofs", {d: DofWithShortRepr() for d in self.dofs}
        if self.mass is not None:
            yield "mass", self.mass, None
        if self.center_of_mass is not None:
            yield "center_of_mass", tuple(self.center_of_mass)
        yield "name", self.name

    def show(self, **kwargs):
        from capytaine.ui.vtk.body_viewer import FloatingBodyViewer
        viewer = FloatingBodyViewer()
        viewer.add_body(self, **kwargs)
        viewer.show()
        viewer.finalize()

    def show_matplotlib(self, *args, **kwargs):
        return self.mesh.show_matplotlib(*args, **kwargs)

    def animate(self, motion, *args, **kwargs):
        """Display a motion as a 3D animation.

        Parameters
        ==========
        motion: dict or pd.Series or str
            A dict or series mapping the name of the dofs to its amplitude.
            If a single string is passed, it is assumed to be the name of a dof
            and this dof with a unit amplitude will be displayed.
        """
        from capytaine.ui.vtk.animation import Animation
        if isinstance(motion, str):
            motion = {motion: 1.0}
        elif isinstance(motion, xr.DataArray):
            motion = {k: motion.sel(radiating_dof=k).data for k in motion.coords["radiating_dof"].data}

        if any(dof not in self.dofs for dof in motion):
            missing_dofs = set(motion.keys()) - set(self.dofs.keys())
            raise ValueError(f"Trying to animate the body {self.name} using dof(s) {missing_dofs}, but no dof of this name is defined for {self.name}.")

        animation = Animation(*args, **kwargs)
        animation._add_actor(self.mesh.merged(), faces_motion=sum(motion[dof_name] * dof for dof_name, dof in self.dofs.items() if dof_name in motion))
        return animation

    @property
    def minimal_computable_wavelength(self):
        """For accuracy of the resolution, wavelength should not be smaller than this value."""
        if self.lid_mesh is not None:
            return max(8*self.mesh.faces_radiuses.max(), 8*self.lid_mesh.faces_radiuses.max())
        else:
            return 8*self.mesh.faces_radiuses.max()

    @lru_cache
    def first_irregular_frequency_estimate(self, *, g=9.81):
        r"""Estimates the angular frequency of the lowest irregular
        frequency.
        This is based on the formula for the lowest irregular frequency of a
        parallelepiped of size :math:`L \times B` and draft :math:`H`:

        .. math::
            \omega = \sqrt{
                        \frac{\pi g \sqrt{\frac{1}{B^2} + \frac{1}{L^2}}}
                             {\tanh\left(\pi H \sqrt{\frac{1}{B^2} + \frac{1}{L^2}} \right)}
                     }

        The formula is applied to all shapes to get an estimate that is usually
        conservative.
        The definition of a lid (supposed to be fully covering and horizontal)
        is taken into account.
        """
        if self.lid_mesh is None:
            draft = abs(self.mesh.vertices[:, 2].min())
        else:
            draft = abs(self.lid_mesh.vertices[:, 2].min())
            if draft < 1e-6:
                return np.inf

        # Look for the x and y span of each components (e.g. for multibody) and
        # keep the one causing the lowest irregular frequency.
        # The draft is supposed to be same for all components.
        omega = np.inf
        for comp in connected_components(self.mesh):
            for ccomp in connected_components_of_waterline(comp):
                x_span = ccomp.vertices[:, 0].max() - ccomp.vertices[:, 0].min()
                y_span = ccomp.vertices[:, 1].max() - ccomp.vertices[:, 1].min()
                p = np.hypot(1/x_span, 1/y_span)
                omega_comp = np.sqrt(np.pi*g*p/(np.tanh(np.pi*draft*p)))
                omega = min(omega, omega_comp)
        return omega
