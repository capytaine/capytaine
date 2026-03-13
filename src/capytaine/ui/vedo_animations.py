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
    cmap: Optional[str] = None
    vmin: Optional[float] = None
    vmax: Optional[float] = None


class Animation:
    vedo = import_optional_dependency("vedo")

    def __init__(self, *, loop_duration, fps=24):
        self.loop_duration = loop_duration
        self.omega = 2 * np.pi / self.loop_duration
        self.fps = fps
        self.components = []


    def add_body(self, body, vertices_pressure=None, vertices_motion=None):
        if vertices_pressure is not None:
            pressure_amplitude = np.abs(vertices_pressure).max()
            cmap_bounds = {"vmin": -pressure_amplitude, "vmax": pressure_amplitude}
        else:
            cmap_bounds = {}
        self.components.append(
                AnimationComponent(
                    mean_mesh=body.mesh,
                    vedo_mesh=body.mesh.export("vedo"),
                    vertices_colors=vertices_pressure,
                    vertices_motion=vertices_motion,
                    cmap="RdBu_r",
                    **cmap_bounds
                    )
                )

    def add_free_surface(self, free_surface_mesh, vertices_elevation):
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
                    vedo_mesh=free_surface_mesh.export("vedo"),
                    vertices_motion=elevation_as_motion,
                    vertices_colors=vertices_elevation,
                    cmap="Blues_r",
                    **cmap_bounds
                    )
                )

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


    def run(self, camera=None, resolution=(800, 600), **kwargs):
        if camera is None:
            camera = {}
        for k, v in _default_camera.items():
            camera.setdefault(k, v)
        kwargs.setdefault("size", resolution)
        plt = self.vedo.Plotter(axes=1, interactive=False, **kwargs)
        plt.show(
            *[comp.vedo_mesh for comp in self.components],
            camera=camera
        )
        t_range = np.linspace(0.0, self.loop_duration, int(self.loop_duration * self.fps))
        for t in t_range:
            self.update(t)
            plt.render()
            sleep(1.0 / self.fps)
        plt.interactive()


    def save(self, filename, camera=None, resolution=(800, 600), **kwargs):
        if camera is None:
            camera = {}
        for k, v in _default_camera.items():
            camera.setdefault(k, v)
        kwargs.setdefault("size", resolution)
        video = self.vedo.Video(filename, duration=self.loop_duration, fps=self.fps)
        plt = self.vedo.Plotter(axes=1, interactive=False, **kwargs)
        plt.show(
            *[comp.vedo_mesh for comp in self.components],
            camera=camera
        )
        t_range = np.linspace(0.0, self.loop_duration, int(self.loop_duration * self.fps))
        for t in t_range:
            self.update(t)
            video.add_frame()
        video.close()
        plt.close()
