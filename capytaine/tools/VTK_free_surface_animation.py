#!/usr/bin/env python
# coding: utf-8
"""
VTK animation for the free surface elevation.

This file is part of "capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import numpy as np
import vtk

from capytaine.results import LinearPotentialFlowResult
from capytaine.geometric_bodies.free_surface import FreeSurface


class Animation:
    def __init__(self,
                 result: LinearPotentialFlowResult,
                 free_surface: FreeSurface,
                 complex_node_elevation=None,
                 fps: int = 24,
                 ):

        if complex_node_elevation is None:
            complex_node_elevation = free_surface.elevation_at_nodes(result.fs_elevation[free_surface])
        assert len(complex_node_elevation) == free_surface.mesh.nb_vertices

        # Create an object for the floating body
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(result.body.mesh._vtk_polydata())

        self.body_actor = vtk.vtkActor()
        self.body_actor.GetProperty().SetColor((1, 1, 0))
        self.body_actor.GetProperty().SetInterpolationToGouraud()
        self.body_actor.SetMapper(mapper)

        # Create an object for the free surface
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(20)
        lut.SetHueRange(0.58, 0.58)
        lut.SetSaturationRange(0.5, 0.5)
        lut.SetValueRange(0.5, 1.0)
        lut.Build()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(free_surface.mesh._vtk_polydata())
        mapper.SetLookupTable(lut)

        self.fs_actor = vtk.vtkActor()
        self.fs_actor.GetProperty().SetInterpolationToGouraud()
        self.fs_actor.SetMapper(mapper)

        # Precompute data before loop
        self.fps = fps
        self.frames_per_period = int(fps * result.period)
        base_polydata = self.fs_actor.GetMapper().GetInput()
        self.polydatas = []

        for i_frame in range(self.frames_per_period):
            current_node_elevation = np.abs(complex_node_elevation) * np.cos(np.angle(complex_node_elevation) - 2 * np.pi * i_frame / self.frames_per_period)

            new_polydata = vtk.vtkPolyData()
            new_polydata.DeepCopy(base_polydata)
            points = new_polydata.GetPoints()
            for j in range(len(complex_node_elevation)):
                point = points.GetPoint(j)
                points.SetPoint(j, (point[0], point[1], current_node_elevation[j]))

            # Store elevation for the LUT
            ef = vtk.vtkSimpleElevationFilter()
            ef.SetInputData(new_polydata)
            ef.Update()
            new_polydata = ef.GetPolyDataOutput()

            self.polydatas.append(new_polydata)

        self.current_frame = 0

    def callback(self, renderer, event):
        # Update body position and free surface elevation
        # update_position(self.body,
        #                 2 * np.pi * (self.current_frame % self.frames_per_period) / self.frames_per_period)

        self.fs_actor.GetMapper().SetInputData(self.polydatas[self.current_frame % self.frames_per_period])

        if hasattr(self, 'writer') and self.current_frame < self.frames_per_period:
            self.imageFilter.Modified()
            self.writer.Write()
            if self.current_frame == self.frames_per_period - 1:
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
        renderer.AddActor(self.body_actor)
        renderer.AddActor(self.fs_actor)

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
