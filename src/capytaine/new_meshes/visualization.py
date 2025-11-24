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

from capytaine.tools.optional_imports import import_optional_dependency
from capytaine.new_meshes.export import mesh_to_pyvista


def show_pyvista(mesh, **kwargs):
    """
    Visualize the mesh using PyVista.

    Parameters
    ----------
    mesh : Mesh
        The mesh object to visualize.
    kwargs : additional optional arguments
        Additional arguments passed to PyVista's add_mesh and plotter methods for customization.
    """
    pv = import_optional_dependency("pyvista")
    pv_mesh = mesh_to_pyvista(mesh.vertices, mesh._faces)
    plotter = pv.Plotter()
    kwargs.setdefault("show_edges", True)
    plotter.add_mesh(pv_mesh, **kwargs)
    plotter.show()


def show_matplotlib(mesh, ax=None, bounding_box=None,
                    normal_vectors=False, scale_normal_vector=None,
                    color_field=None, cmap=None,
                    cbar_label=None,
                    **kwargs):
    """
    Visualize the mesh using Matplotlib.

    Parameters
    ----------
    ax: matplotlib axis
        The 3d axis in which to plot the mesh. If not provided, create a new one.
    bounding_box: tuple[tuple[int]], optional
        Min and max coordinates values to display in each three dimensions.
    normal_vectors: bool, optional
        If True, print normal vector.
    scale_normal_vector: array of shape (nb_faces, ), optional
        Scale separately each of the normal vectors.
    color_field: array of shape (nb_faces, ), optional
        Scalar field to be plot on the mesh (optional).
    cmap: matplotlib colormap, optional
        Colormap to use for field plotting.
    cbar_label: string, optional
        Label for colormap

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

    faces = []
    for face in mesh.faces:
        vertices = [mesh.vertices[int(index_vertex), :] for index_vertex in face]
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
