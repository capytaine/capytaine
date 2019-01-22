#!/usr/bin/env python
# coding: utf-8
"""VTK animation for the free surface elevation."""
# This file is part of "capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import logging

import numpy as np
from numpy import pi
import vtk

from capytaine.ui.vtk.helpers import compute_node_data, compute_vtk_polydata

LOG = logging.getLogger(__name__)


class Animation:
    """Class to generate an animation of a result of Capytaine,
    including the elevation of the free surface.

    The animation is made of a short loop of a single period of the solution in frequency domain.

    Parameters
    ----------
    loop_duration: float
        Duration in the loop. For real time animation, the period of the motion.
    fps: int, optional
        Number of frames per second in the animation (default: 24).

    Attributes
    ----------
    frames_per_loop: int
        Number of frames in one loop
    actors: list of vtk actor objects
        The objects in the scene.
    """

    def __init__(self, loop_duration, fps=24):
        self.fps = fps
        self.frames_per_loop = int(fps * loop_duration)
        self.actors = []

        self._precomputed_polydatas = {}
        self._current_frame = 0

    def _add_actor(self, mesh, faces_motion=None, lut=False, edges=False):
        """Add an animated object to the scene."""
        if faces_motion is not None:
            nodes_motion = compute_node_data(mesh.merged(), faces_motion)
        else:
            nodes_motion = None

        base_polydata = compute_vtk_polydata(mesh.merged())
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(base_polydata)

        actor = vtk.vtkActor()
        if edges:
            actor.GetProperty().EdgeVisibilityOn()
        actor.GetProperty().SetInterpolationToGouraud()
        actor.SetMapper(mapper)

        if nodes_motion is not None:
            LOG.info(f"Precompute motions of {mesh.name} before animation.")
            self._precomputed_polydatas[actor] = []

            for i_frame in range(self.frames_per_loop):
                current_deformation = np.abs(nodes_motion) * np.cos(np.angle(nodes_motion)
                                                                    - 2*pi*i_frame/self.frames_per_loop)

                new_polydata = vtk.vtkPolyData()
                new_polydata.DeepCopy(base_polydata)
                points = new_polydata.GetPoints()
                for j in range(mesh.nb_vertices):
                    point = points.GetPoint(j)
                    point = np.asarray(point) + current_deformation[j]
                    points.SetPoint(j, tuple(point))

                # Store elevation for the LUT of the free surface
                if lut:
                    ef = vtk.vtkElevationFilter()
                    ef.SetLowPoint(0, 0, min(abs(nodes_motion[2])))
                    ef.SetHighPoint(0, 0, max(abs(nodes_motion[2])))
                    ef.SetInputData(new_polydata)
                    ef.Update()
                    new_polydata = ef.GetPolyDataOutput()

                self._precomputed_polydatas[actor].append(new_polydata)
        else:
            self._precomputed_polydatas[actor] = None

        self.actors.append(actor)

        return actor

    def add_body(self, body, faces_motion=None, edges=False):
        """Add an floating body to the scene.

        Parameters
        ----------
        body: FloatingBody
            The object to include in the scene.
        faces_motion: dof, optional
            The motion of the body defined at the center of the faces.
        edges: bool, optional
            Draw the edges of the mesh in the scene.

        Returns
        -------
        vtk actor object
        """
        actor = self._add_actor(body.mesh.merged(), faces_motion=faces_motion,
                                edges=edges)
        actor.GetProperty().SetColor((1, 1, 0))
        return actor

    def add_free_surface(self, free_surface, faces_elevation):
        """Add the free surface to the scene.

        Parameters
        ----------
        free_surface: FreeSurface
            The free surface object
        faces_elevation: array of complex numbers
            The elevation of each face of the meshed free surface given as a complex number.

        Returns
        -------
        vtk actor object
        """
        faces_motion = np.array([(0, 0, elevation) for elevation in faces_elevation])
        actor = self._add_actor(free_surface.mesh, faces_motion=faces_motion, lut=True)

        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(50)
        lut.SetHueRange(0.58, 0.58)
        lut.SetSaturationRange(0.5, 0.5)
        lut.SetValueRange(0.5, 0.6)
        lut.Build()
        actor.GetMapper().SetLookupTable(lut)

        return actor

    # def add_result(self, result):
    #     from capytaine.bodies import TRANSLATION_DOFS_DIRECTIONS, ROTATION_DOFS_AXIS
    #
    #         if display_dof.lower() in TRANSLATION_DOFS_DIRECTIONS:
    #             direction = np.asarray(TRANSLATION_DOFS_DIRECTIONS[display_dof.lower()])
    #             def translation_motion(self, frame):
    #                 nonlocal direction
    #                 pos = np.asarray(self.body_actor.GetPosition())
    #                 pos = (1 - direction) * pos + \
    #                       direction * np.cos(2*np.pi*(frame % self.frames_per_loop)/self.frames_per_loop)
    #                 self.body_actor.SetPosition(*pos)
    #             self.update_body_position = translation_motion
    #
    #         elif display_dof.lower() in ROTATION_DOFS_AXIS:
    #             direction = np.asarray(ROTATION_DOFS_AXIS[display_dof.lower()])
    #             def rotation_motion(self, frame):
    #                 nonlocal direction
    #                 pos = np.asarray(self.body_actor.GetOrientation())
    #                 pos = (1 - direction) * pos + \
    #                       direction * np.cos(2*np.pi*(frame % self.frames_per_loop)/self.frames_per_loop)
    #                 self.body_actor.SetOrientation(*pos)
    #             self.update_body_position = rotation_motion

    def _callback(self, renderer, event):
        for actor in self.actors:
            if self._precomputed_polydatas[actor] is not None:
                actor.GetMapper().SetInputData(self._precomputed_polydatas[actor][self._current_frame % self.frames_per_loop])

        if hasattr(self, 'writer') and self._current_frame < self.frames_per_loop:
            self.imageFilter.Modified()
            self.writer.Write()
            if self._current_frame == self.frames_per_loop - 1:
                self.writer.End()
                render_window = renderer.GetRenderWindow()
                render_window.Finalize()
                renderer.TerminateApp()

        renderer.GetRenderWindow().Render()

        self._current_frame += 1

    def run(self, camera_position=(10.0, 10.0, 10.0), resolution=(1280, 720), out_file_path=None):
        """Run the animation.

        Parameters
        ----------
        camera_position: 3-ple of floats, optional
            The starting position of the camera in the scene.
        resolution: 2-ple of ints, optional
            Resolution of the video in pixels.
        out_file_path: string, optional
            File in which to save the animation
        """
        # Setup a renderer, render window, and interactor
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(1, 1, 1)  # Background color white
        for actor in self.actors:
            renderer.AddActor(actor)

        camera = vtk.vtkCamera()
        camera.SetPosition(*camera_position)
        camera.SetFocalPoint(0, 0, 0)
        camera.SetViewUp(0, 0, 1)
        renderer.SetActiveCamera(camera)

        render_window = vtk.vtkRenderWindow()
        render_window.SetSize(*resolution)
        render_window.SetWindowName("Capytaine animation")
        render_window.AddRenderer(renderer)

        render_window_interactor = vtk.vtkRenderWindowInteractor()
        render_window_interactor.SetRenderWindow(render_window)
        render_window_interactor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

        if out_file_path is not None:
            self.imageFilter = vtk.vtkWindowToImageFilter()
            self.imageFilter.SetInput(render_window)
            self.imageFilter.SetInputBufferTypeToRGB()
            self.imageFilter.ReadFrontBufferOff()

            self.writer = vtk.vtkOggTheoraWriter()
            self.writer.SetInputConnection(self.imageFilter.GetOutputPort())
            self.writer.SetFileName(out_file_path)
            self.writer.SetRate(self.fps)

        # Render and interact
        render_window.Render()
        render_window_interactor.Initialize()  # Initialize must be called prior to creating timer events.

        render_window_interactor.AddObserver('TimerEvent', self._callback)
        render_window_interactor.CreateRepeatingTimer(int(1000 / self.fps))

        if out_file_path is not None:
            self.writer.Start()

        # Start the interaction and timer
        render_window_interactor.Start()

        if out_file_path is not None:
            del self.imageFilter
            del self.writer

        del render_window_interactor
        del render_window

