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

import importlib
from typing import Optional, List

import numpy as np

from capytaine import __version__
from capytaine.tools.optional_imports import import_optional_dependency


def show_3d(mesh, *, backend=None, **kwargs):
    """Dispatch the 3D viewing to one of the available backends below."""
    backends_functions = {
            "pyvista": show_pyvista,
            "matplotlib": show_matplotlib,
            }
    if backend is not None:
        if backend in backends_functions:
            return backends_functions[backend](mesh, **kwargs)
        else:
            raise NotImplementedError(f"Backend '{backend}' is not implemented.")
    else:
        for backend in backends_functions:
            try:
                return backends_functions[backend](mesh, **kwargs)
            except (NotImplementedError, ImportError):
                pass
        raise NotImplementedError(f"No compatible backend found to show the mesh {mesh}.\n"
                                  "Consider installing `matplotlib` or `pyvista`.")


def show_pyvista(mesh, **kwargs) -> Optional["pv.Plotter"]:  # noqa: F821
    """Visualize the mesh using PyVista. See :class:`PyVistaViewer` for parameters."""
    viewer = PyVistaViewer(mesh, **kwargs)
    return viewer.show()


class PyVistaViewer:
    """Interactive PyVista mesh viewer with toggleable features.

    PyVista default keyboard controls: https://docs.pyvista.org/api/plotting/plotting

    Parameters
    ----------
    mesh : Mesh
        The mesh object to visualize.
    ghost_meshes: List[Mesh], optional
        Additional meshes, shown in transparency, used for symmetric meshes
    plotter: pv.Plotter, optional
        If provided, use this PyVista plotter and return it at the end.
        Otherwise a new one is created and the 3D view is displayed at the end.
    normal_vectors: bool, optional
        If True, display normal vector (default: False)
    display_free_surface: bool, optional
        If True, display free surface and if ``water_depth`` is finite display the sea bottom.
        (default: True)
    water_depth: float, optional
        Where to display the sea bottom if ``display_free_surface`` is True
    color_field: array of shape (nb_faces, ), optional
        Scalar field to be plot on the mesh.
    cmap: matplotlib colormap, optional
        Colormap to use for scalar field plotting.
    cbar_label: string, optional
        Label for colorbar showing color field scale
    mesh_kwargs : additional optional arguments
        Additional arguments passed to PyVista's add_mesh methods for customization (e.g. mesh color).
    """

    help_text = [
        f"Capytaine version {__version__}\n",
        "Keyboard controls:",
        "  b: toggle scale and bounding box",
        "  h: toggle free surface (and sea bottom if water depth was given)",
        "  n: toggle normal vectors",
        "  T,B,P,S,F,R: view [T]op, [B]ottom, [P]ort, [S]tarboard, [F]ront, [R]ear",
        "  X, Y: toggle displaying clipped mesh in x or y direction",
        "  q: exit"
    ]

    def __init__(
        self,
        mesh,
        *,
        ghost_meshes=None,
        plotter=None,
        normal_vectors=False,
        display_free_surface=True,
        water_depth=np.inf,
        color_field=None,
        cbar_label=None,
        **mesh_kwargs,
    ):
        pv = import_optional_dependency("pyvista")

        all_meshes = [mesh] if ghost_meshes is None else [mesh, *ghost_meshes]

        self.mesh = mesh
        self.water_depth = water_depth
        self.mesh_kwargs = mesh_kwargs
        self.pv_meshes = [m.export_to_pyvista() for m in all_meshes]

        # Distribute color field across sub-meshes
        if color_field is not None and isinstance(color_field, np.ndarray):
            acc = 0
            for m in self.pv_meshes:
                m.cell_data["color_field"] = color_field[acc:acc + m.n_cells]
                acc += m.n_cells

        # Plotter
        self._owns_plotter = plotter is None
        self.plotter = pv.Plotter() if plotter is None else plotter

        # Scene extents
        self._scene_min = np.min([m.vertices.min(axis=0) for m in all_meshes], axis=0)
        self._scene_max = np.max([m.vertices.max(axis=0) for m in all_meshes], axis=0)

        # Build scene
        if color_field is not None:
            self.mesh_kwargs.setdefault("scalars", "color_field")
        if cbar_label is not None:
            self.mesh_kwargs.setdefault("scalar_bar_args", {"title": cbar_label})

        self.plotter.add_mesh(
            self.pv_meshes[0], name="hull", show_edges=True, **self.mesh_kwargs,
        )
        for i, g_mesh in enumerate(self.pv_meshes[1:]):
            self.plotter.add_mesh(
                g_mesh, name=f"symmetric_hull_{i}",
                opacity=0.4, show_edges=False, **self.mesh_kwargs,
            )

        # Toggle states
        self._normals_visible = normal_vectors
        if self._normals_visible:
            self._show_normals()
        self._free_surface_visible = display_free_surface
        if self._free_surface_visible:
            self._show_free_surface()
        self._bounds_visible = True
        if self._bounds_visible:
            self._show_bounds()
        self._clipping = {'x': 0, 'y': 0}

        p = self.plotter
        p.add_key_event("n", self.toggle_normals)
        p.add_key_event("h", self.toggle_free_surface)
        p.add_key_event("b", self.toggle_bounds)
        p.add_key_event("T", p.view_xy)
        p.add_key_event("B", lambda: p.view_xy(negative=True))
        p.add_key_event("S", p.view_xz)
        p.add_key_event("P", lambda: p.view_xz(negative=True))
        p.add_key_event("F", p.view_yz)
        p.add_key_event("R", lambda: p.view_yz(negative=True))
        p.add_key_event("X", lambda: self.toggle_clipping("x"))
        p.add_key_event("Y", lambda: self.toggle_clipping("y"))

        self.plotter.add_text("\n".join(self.help_text), position="upper_left", font_size=10)

        self.plotter.show_axes()

    # -- Normal vectors --

    def _show_normals(self):
        mini = self.mesh.vertices.min()
        maxi = self.mesh.vertices.max()
        self.plotter.add_arrows(
            self.mesh.faces_centers, self.mesh.faces_normals,
            name="normals", mag=0.04 * (maxi - mini), show_scalar_bar=False,
        )

    def toggle_normals(self):
        if self._normals_visible:
            self.plotter.remove_actor("normals")
            self._normals_visible = False
        else:
            self._show_normals()
            self._normals_visible = True

    # -- Free surface --

    def _show_free_surface(self):
        pv = import_optional_dependency("pyvista")
        center = (self._scene_min[:2] + self._scene_max[:2]) / 2
        diam = 1.1 * (self._scene_max[:2] - self._scene_min[:2])
        plane = pv.Plane(
            center=(*center, 0), direction=(0, 0, 1),
            i_size=diam[0], j_size=diam[1],
        )
        self.plotter.add_mesh(plane, color="blue", opacity=0.5, name="display_free_surface")
        if self.water_depth != np.inf:
            plane = pv.Plane(
                center=(*center, -self.water_depth), direction=(0, 0, 1),
                i_size=diam[0], j_size=diam[1],
            )
            self.plotter.add_mesh(plane, color="brown", opacity=0.5, name="display_sea_bottom")

    def toggle_free_surface(self):
        if self._free_surface_visible:
            self.plotter.remove_actor("display_free_surface")
            if self.water_depth != np.inf:
                self.plotter.remove_actor("display_sea_bottom")
            self._free_surface_visible = False
        else:
            self._show_free_surface()
            self._free_surface_visible = True

    # -- Bounds --

    def _show_bounds(self):
        self.plotter.show_bounds(
            grid="back", location="outer", n_xlabels=2, n_ylabels=2, n_zlabels=2,
        )

    def toggle_bounds(self):
        if self._bounds_visible:
            self.plotter.remove_bounds_axes()
            self._bounds_visible = False
        else:
            self._show_bounds()
            self.plotter.update()
            self._bounds_visible = True

    # -- Clipping --

    def _clipped_mesh(self):
        result = self.pv_meshes[0]
        for axis in ["x", "y"]:
            if self._clipping[axis] == 1:
                result = result.clip(axis)
            elif self._clipping[axis] == -1:
                result = result.clip("-" + axis)
        return result

    def toggle_clipping(self, direction):
        """Cycle clipping for *direction* (``'x'`` or ``'y'``): none → positive → negative → none."""
        c = self._clipping[direction]
        self._clipping[direction] = {0: 1, 1: -1, -1: 0}[c]
        self.plotter.add_mesh(
            self._clipped_mesh(), name="hull", show_edges=True, **self.mesh_kwargs,
        )

    # ── Show ────────────────────────────────────────────────────

    def show(self):
        """Display the viewer window (if plotter was created internally) and return plotter."""
        if self._owns_plotter:
            self.plotter.show()
        return self.plotter


###################################################################################################

def show_matplotlib(
    mesh,
    *,
    ghost_meshes=None,
    ax=None,
    bounding_box=None,
    normal_vectors=False,
    scale_normal_vector=None,
    color_field=None,
    cmap=None,
    cbar_label=None,
    **kwargs
):
    """
    Visualize the mesh using Matplotlib.

    Parameters
    ----------
    mesh : Mesh
        The mesh object to visualize.
    ghost_meshes: List[Mesh], optional
        Additional meshes. In the matplotlib viewer, they are just merged with the main mesh.
    ax: matplotlib axis
        The 3d axis in which to plot the mesh. If not provided, create a new one.
    bounding_box: tuple[tuple[int]], optional
        Min and max coordinates values to display in each three dimensions.
    normal_vectors: bool, optional
        If True, display normal vector.
    scale_normal_vector: array of shape (nb_faces, ), optional
        Scale separately each of the normal vectors.
    color_field: array of shape (nb_faces, ), optional
        Scalar field to be plot on the mesh (optional).
    cmap: matplotlib colormap, optional
        Colormap to use for scalar field plotting.
    cbar_label: string, optional
        Label for colorbar show color field scale

    Other parameters are passed to Poly3DCollection.
    """
    matplotlib = import_optional_dependency("matplotlib")
    plt = importlib.import_module("matplotlib.pyplot")
    cm = importlib.import_module("matplotlib.cm")

    mpl_toolkits = import_optional_dependency("mpl_toolkits", package_name="matplotlib")
    Poly3DCollection = mpl_toolkits.mplot3d.art3d.Poly3DCollection

    default_axis = ax is None
    if default_axis:
        fig = plt.figure(layout="constrained")
        ax = fig.add_subplot(111, projection="3d")
        ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio

    from capytaine.meshes import Mesh  # Just for the typehint below
    all_meshes_in_scene: List[Mesh] = [mesh] if ghost_meshes is None else [mesh, *ghost_meshes]

    faces = []
    for m in all_meshes_in_scene:
        for face in m.faces:
            vertices = [m.vertices[int(index_vertex), :] for index_vertex in face]
            faces.append(vertices)

    if color_field is None:
        if 'facecolors' not in kwargs:
            kwargs['facecolors'] = "yellow"
    else:
        if cmap is None:
            cmap = matplotlib.colormaps['coolwarm']
        m = cm.ScalarMappable(cmap=cmap)
        m.set_array([min(color_field), max(color_field)])
        m.set_clim(vmin=min(color_field), vmax=max(color_field))
        colors = m.to_rgba(color_field)
        kwargs['facecolors'] = colors
    kwargs.setdefault("edgecolor", "k")

    ax.add_collection3d(Poly3DCollection(faces, **kwargs))

    if color_field is not None:
        cbar = plt.colorbar(m, ax=ax)
        if cbar_label is not None:
            cbar.set_label(cbar_label)

    # Plot normal vectors.
    if normal_vectors:
        if scale_normal_vector is not None:
            vectors = (scale_normal_vector * mesh.faces_normals.T).T
        else:
            vectors = mesh.faces_normals
        ax.quiver(*zip(*mesh.faces_centers), *zip(*vectors), length=0.2)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    if bounding_box is None:
        # auto cube around mesh
        mini = mesh.vertices.min(axis=0)
        maxi = mesh.vertices.max(axis=0)
        center = (mini + maxi) / 2
        radius = (maxi - mini).max() / 2
        ax.set_xlim(center[0] - radius, center[0] + radius)
        ax.set_ylim(center[1] - radius, center[1] + radius)
        ax.set_zlim(center[2] - radius, center[2] + radius)
    else:
        (xmin, xmax), (ymin, ymax), (zmin, zmax) = bounding_box
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_zlim(zmin, zmax)

    if default_axis:
        plt.show()
    return ax
