# Copyright 2026 Capytaine developers
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
from time import sleep
from dataclasses import dataclass

import numpy as np

from capytaine.meshes.abstract_meshes import AbstractMesh
from capytaine.tools.optional_imports import import_optional_dependency


_default_camera = dict(
    pos=(-20.0, +20.0, -20.0),
    viewup=(0.0, 0.0, 1.0),
    focal_point=(0.0, 0.0, 0.0),
    )


@dataclass
class AnimationComponent:
    mean_mesh: AbstractMesh  # The mean position of the mesh, used as a reference for the motion
    vedo_mesh: "vedo.Mesh"  # The mesh to be animated
    vertices_motion: Optional[np.ndarray] = None  # Complex array of shape (n_vertices, 3) representing the motion of each vertex as a function of time
    vertices_colors: Optional[np.ndarray] = None  # Complex array of shape (n_vertices,) representing the color of each vertex as a function of time
    cmap: Optional[str] = None  # Name of a colormap, a matplotlib colormap or a list of colors, as accepted by vedo.Mesh.cmap
    vmin: Optional[float] = None
    vmax: Optional[float] = None


class Animation:
    vedo = import_optional_dependency("vedo")

    def __init__(self, *, loop_duration, fps=24):
        self.loop_duration = loop_duration
        self.omega = 2 * np.pi / self.loop_duration
        self.fps = fps
        self.components = []


    def _time_steps(self):
        # The end of the loop is excluded, since it would be the same frame as the first one.
        return np.linspace(0.0, self.loop_duration, int(self.loop_duration * self.fps), endpoint=False)

    def add_body(self, body, vertices_pressure=None, vertices_motion=None, *, color=None, edge_color=None, cmap="RdBu_r"):
        """Add a body to the animation.

        Parameters
        ----------
        body: FloatingBody
            The body to be displayed
        vertices_pressure: array of complex of shape (n_vertices,), optional
            Complex amplitude of a scalar field displayed with a colormap on the body.
        vertices_motion: array of complex of shape (n_vertices, 3), optional
            Complex amplitude of the motion of each vertex of the mesh.
        color: optional
            Uniform color of the body (any color accepted by vedo, e.g. ``"#21355d"``).
            Ignored if ``vertices_pressure`` is given.
        edge_color: optional
            If given, the edges of the faces of the mesh are drawn with this color.
        cmap: optional
            Colormap used for ``vertices_pressure``: name of a colormap, a matplotlib colormap object or a list of colors.
        """
        if vertices_pressure is not None:
            pressure_amplitude = np.abs(vertices_pressure).max()
            cmap_bounds = {"vmin": -pressure_amplitude, "vmax": pressure_amplitude}
        else:
            cmap_bounds = {}
        self.components.append(
                AnimationComponent(
                    mean_mesh=body.mesh,
                    vedo_mesh=self._vedo_mesh(body.mesh, color=color, edge_color=edge_color),
                    vertices_colors=vertices_pressure,
                    vertices_motion=vertices_motion,
                    cmap=cmap,
                    **cmap_bounds
                    )
                )

    def add_free_surface(self, free_surface_mesh, vertices_elevation, *, color=None, edge_color=None, cmap="Blues_r"):
        """Add a free surface to the animation.

        Parameters
        ----------
        free_surface_mesh: Mesh
            The mesh of the free surface
        vertices_elevation: array of complex of shape (n_vertices,) or None
            Complex amplitude of the free surface elevation at each vertex.
            It is used both for the vertical motion and for the colors.
        color: optional
            Uniform color of the free surface (any color accepted by vedo).
            Ignored if ``vertices_elevation`` is given.
        edge_color: optional
            If given, the edges of the faces of the mesh are drawn with this color.
        cmap: optional
            Colormap used for the elevation: name of a colormap, a matplotlib colormap object or a list of colors.
        """
        if vertices_elevation is not None:
            elevation_amplitude = np.abs(vertices_elevation).max()
            cmap_bounds = {"vmin": -elevation_amplitude, "vmax": elevation_amplitude}
            elevation_as_motion = np.vstack(
                [np.zeros_like(vertices_elevation),
                np.zeros_like(vertices_elevation),
                vertices_elevation]
            ).T
        else:
            cmap_bounds = {}
            elevation_amplitude = None
            elevation_as_motion = None
        self.components.append(
                AnimationComponent(
                    mean_mesh=free_surface_mesh,
                    vedo_mesh=self._vedo_mesh(free_surface_mesh, color=color, edge_color=edge_color),
                    vertices_motion=elevation_as_motion,
                    vertices_colors=vertices_elevation,
                    cmap=cmap,
                    **cmap_bounds
                    )
                )

    @staticmethod
    def _vedo_mesh(mesh, *, color=None, edge_color=None):
        vedo_mesh = mesh.export("vedo")
        if color is not None:
            vedo_mesh.color(color)
        if edge_color is not None:
            vedo_mesh.linecolor(edge_color).linewidth(1)
        return vedo_mesh

    def update(self, t):
        for component in self.components:
            if component.vertices_motion is not None:
                motion = np.real(component.vertices_motion * np.exp(-1j * self.omega * t))
                new_pts = component.vedo_mesh.vertices
                new_pts[:, :] = component.mean_mesh.vertices + motion[:, :]
                component.vedo_mesh.vertices = new_pts
            if component.vertices_colors is not None:
                colors = np.real(component.vertices_colors * np.exp(-1j * self.omega * t)).ravel()
                component.vedo_mesh.cmap(component.cmap, colors, vmin=component.vmin, vmax=component.vmax)


    def run(self, camera=None, lights=[], resolution=(800, 600), **kwargs):
        if camera is None:
            camera = {}
        for k, v in _default_camera.items():
            camera.setdefault(k, v)
        kwargs.setdefault("size", resolution)
        kwargs.setdefault("axes", 1)
        plt = self.vedo.Plotter(interactive=False, **kwargs)
        plt.show(
            *[comp.vedo_mesh for comp in self.components],
            *lights,
            camera=camera
        )
        t_range = self._time_steps()
        for t in t_range:
            self.update(t)
            plt.render()
            sleep(1.0 / self.fps)
        plt.interactive()


    def save(self, filename, camera=None, lights=[], resolution=(800, 600), **kwargs):
        if camera is None:
            camera = {}
        for k, v in _default_camera.items():
            camera.setdefault(k, v)
        kwargs.setdefault("size", resolution)
        video = self.vedo.Video(filename, duration=self.loop_duration, fps=self.fps)
        kwargs.setdefault("axes", 1)
        plt = self.vedo.Plotter(interactive=False, **kwargs)
        plt.show(
            *[comp.vedo_mesh for comp in self.components],
            *lights,
            camera=camera
        )
        t_range = self._time_steps()
        for t in t_range:
            self.update(t)
            video.add_frame()
        video.close()
        plt.close()
