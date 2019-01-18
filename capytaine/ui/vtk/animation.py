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

    def __init__(self, loop_duration, fps=24):
        self.fps = fps
        self.frames_per_loop = int(fps * loop_duration)

        self.actors = []
        self._precomputed_polydatas = {}

        self.current_frame = 0

    def add_actor(self, mesh, faces_motion=None, lut=False):
        if faces_motion is not None:
            nodes_motion = compute_node_data(mesh.merged(), faces_motion)
        else:
            nodes_motion = None

        base_polydata = compute_vtk_polydata(mesh.merged())
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(base_polydata)

        actor = vtk.vtkActor()
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
                    ef = vtk.vtkSimpleElevationFilter()
                    ef.SetInputData(new_polydata)
                    ef.Update()
                    new_polydata = ef.GetPolyDataOutput()

                self._precomputed_polydatas[actor].append(new_polydata)
        else:
            self._precomputed_polydatas[actor] = None

        self.actors.append(actor)

        return actor

    def add_body(self, body, faces_motion=None):
        actor = self.add_actor(body.mesh.merged(), faces_motion=faces_motion)
        actor.GetProperty().SetColor((1, 1, 0))
        return actor

    def add_free_surface(self, free_surface, faces_elevation):
        faces_motion = np.array([(0, 0, elevation) for elevation in faces_elevation])
        actor = self.add_actor(free_surface.mesh, faces_motion=faces_motion, lut=True)

        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(20)
        lut.SetHueRange(0.58, 0.58)
        lut.SetSaturationRange(0.5, 0.5)
        lut.SetValueRange(0.5, 1.0)
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

    def callback(self, renderer, event):
        for actor in self.actors:
            if self._precomputed_polydatas[actor] is not None:
                actor.GetMapper().SetInputData(self._precomputed_polydatas[actor][self.current_frame % self.frames_per_loop])

        if hasattr(self, 'writer') and self.current_frame < self.frames_per_loop:
            self.imageFilter.Modified()
            self.writer.Write()
            if self.current_frame == self.frames_per_loop - 1:
                self.writer.End()
                render_window = renderer.GetRenderWindow()
                render_window.Finalize()
                renderer.TerminateApp()

        renderer.GetRenderWindow().Render()

        self.current_frame += 1

    def run(self, out_file_path=None):
        # Setup a renderer, render window, and interactor
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(1, 1, 1)  # Background color white
        for actor in self.actors:
            renderer.AddActor(actor)

        camera = vtk.vtkCamera()
        camera.SetPosition(80, 80, 90)
        camera.SetFocalPoint(0, 0, 0)
        camera.SetViewUp(0, 0, 1)
        renderer.SetActiveCamera(camera)

        render_window = vtk.vtkRenderWindow()
        render_window.SetSize(1024, 768)
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

        render_window_interactor.AddObserver('TimerEvent', self.callback)
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

