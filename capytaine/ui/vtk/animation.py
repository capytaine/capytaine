#!/usr/bin/env python
# coding: utf-8
"""VTK animation for the free surface elevation."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np
from numpy import pi

from capytaine.ui.vtk.helpers import compute_node_data, compute_vtk_polydata
from capytaine.tools.optional_imports import import_optional_dependency

vtk = import_optional_dependency("vtk")

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

    # @classmethod
    # def from_result(self, result):
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


    def _add_actor(self, mesh, faces_motion=None, faces_colors=None, edges=False):
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

        if nodes_motion is not None or faces_colors is not None:
            LOG.info(f"Precompute motions of {mesh.name} before animation.")
            self._precomputed_polydatas[actor] = []

            for i_frame in range(self.frames_per_loop):
                new_polydata = vtk.vtkPolyData()
                new_polydata.DeepCopy(base_polydata)

                if nodes_motion is not None:
                    # Change points positions at frame i
                    current_deformation = (
                            np.abs(nodes_motion)*np.cos(np.angle(nodes_motion)-2*pi*i_frame/self.frames_per_loop)
                    )

                    points = new_polydata.GetPoints()
                    for j in range(mesh.nb_vertices):
                        point = points.GetPoint(j)
                        point = np.asarray(point) + current_deformation[j]
                        points.SetPoint(j, tuple(point))

                if faces_colors is not None:
                    # Evaluate scalar field at frame i
                    current_colors = (
                        np.abs(faces_colors)*np.cos(np.angle(faces_colors)-2*pi*i_frame/self.frames_per_loop)
                    )
                    max_val = max(abs(faces_colors))
                    vtk_faces_colors = vtk.vtkFloatArray()
                    for i, color in enumerate(current_colors):
                        vtk_faces_colors.InsertValue(i, (color+max_val)/(2*max_val))
                    new_polydata.GetCellData().SetScalars(vtk_faces_colors)

                self._precomputed_polydatas[actor].append(new_polydata)
        else:
            self._precomputed_polydatas[actor] = None

        self.actors.append(actor)

        return actor

    def add_body(self, body, faces_motion=None, faces_colors=None, edges=False):
        """Add an floating body to the scene.

        Parameters
        ----------
        body: FloatingBody
            The object to include in the scene.
        faces_motion: dof, optional
            The motion of the body defined at the center of the faces.
        faces_colors: iterable of complex numbers, optional
            Scalar field over the surface of the body that should be displayed with colors.
        edges: bool, optional
            Draw the edges of the mesh in the scene.

        Returns
        -------
        vtk actor object
        """
        actor = self._add_actor(body.mesh.merged(), faces_motion=faces_motion,
                                faces_colors=faces_colors, edges=edges)
        if faces_colors is None:
            actor.GetProperty().SetColor((1, 1, 0))
        else:
            lut = vtk.vtkLookupTable()
            lut.SetNumberOfColors(50)
            lut.SetHueRange(0, 0.6)
            lut.SetSaturationRange(0.5, 0.5)
            lut.SetValueRange(0.8, 0.8)
            lut.Build()
            actor.GetMapper().SetLookupTable(lut)

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
        actor = self._add_actor(free_surface.mesh, faces_motion=faces_motion, faces_colors=faces_motion[:, 2])

        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(50)
        lut.SetHueRange(0.58, 0.58)
        lut.SetSaturationRange(0.5, 0.5)
        lut.SetValueRange(0.4, 0.6)
        lut.Build()
        actor.GetMapper().SetLookupTable(lut)

        return actor

    def _callback(self, renderer, event):
        for actor in self.actors:
            if self._precomputed_polydatas[actor] is not None:
                actor.GetMapper().SetInputData(self._precomputed_polydatas[actor][self._current_frame % self.frames_per_loop])
        renderer.GetRenderWindow().Render()
        self._current_frame += 1

    def run(self, camera_position=(-10.0, -10.0, 10.0), resolution=(1280, 720)):
        """Run the animation.

        Parameters
        ----------
        camera_position: 3-ple of floats, optional
            The starting position of the camera in the scene.
        resolution: 2-ple of ints, optional
            Resolution of the video in pixels.
        """
        # Setup a renderer, render window, and interactor
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(1, 1, 1)  # Background color white
        for actor in self.actors:
            renderer.AddActor(actor)
        renderer.Modified()

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

        render_window.Render()

        render_window_interactor.Initialize()  # Initialize must be called prior to creating timer events.

        render_window_interactor.AddObserver('TimerEvent', self._callback)
        render_window_interactor.CreateRepeatingTimer(int(1000 / self.fps))

        render_window_interactor.Start()

        # Run until stopped by user.

        del render_window_interactor
        del render_window

    def save(self, filepath, nb_loops=1, camera_position=(-10.0, -10.0, 10.0), resolution=(1280, 720)):
        """Save the animation in a video file.

        Parameters
        ----------
        filepath: string
            Path of the output file.
        nb_loop: int, optional
            Number of periods to save in the file.
        camera_position: 3-ple of floats, optional
            The starting position of the camera in the scene.
        resolution: 2-ple of ints, optional
            Resolution of the video in pixels.
        """
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(1, 1, 1)  # Background color white
        for actor in self.actors:
            renderer.AddActor(actor)
        renderer.Modified()

        camera = vtk.vtkCamera()
        camera.SetPosition(*camera_position)
        camera.SetFocalPoint(0, 0, 0)
        camera.SetViewUp(0, 0, 1)
        renderer.SetActiveCamera(camera)

        render_window = vtk.vtkRenderWindow()
        render_window.SetSize(*resolution)
        render_window.OffScreenRenderingOn()
        render_window.AddRenderer(renderer)

        image_filter = vtk.vtkWindowToImageFilter()
        image_filter.SetInput(render_window)
        image_filter.SetInputBufferTypeToRGB()
        image_filter.ReadFrontBufferOff()

        writer = vtk.vtkOggTheoraWriter()
        writer.SetInputConnection(image_filter.GetOutputPort())
        writer.SetFileName(filepath)
        writer.SetRate(self.fps)

        writer.Start()

        for i_frame in range(nb_loops*self.frames_per_loop):
            self._callback(renderer, None)
            image_filter.Modified()
            writer.Write()

        writer.End()
        render_window.Finalize()

        del image_filter
        del writer
        del render_window

    def embed_in_notebook(self, resolution=(640, 360), **kwargs):
        from tempfile import mkstemp
        from IPython.core.display import Video

        class CustomIPythonVideo(Video):
            def __init__(self, *args, html_attributes="controls", **kwargs):
                self.html_attributes = html_attributes
                super(CustomIPythonVideo, self).__init__(*args, **kwargs)

            def _repr_html_(self):
                import mimetypes
                from binascii import b2a_hex, b2a_base64, hexlify

                width = height = ''
                if self.width:
                    width = ' width="%d"' % self.width
                if self.height:
                    height = ' height="%d"' % self.height

                # External URLs and potentially local files are not embedded into the
                # notebook output.
                if not self.embed:
                    url = self.url if self.url is not None else self.filename
                    output = """<video src="{0}" {1} {2} {3}>
Your browser does not support the <code>video</code> element.
</video>""".format(url, self.html_attributes, width, height)
                    return output

                # Embedded videos are base64-encoded.
                mimetype = self.mimetype
                if self.filename is not None:
                    if not mimetype:
                        mimetype, _ = mimetypes.guess_type(self.filename)

                    with open(self.filename, 'rb') as f:
                        video = f.read()
                else:
                    video = self.data
                if isinstance(video, str):
                    # unicode input is already b64-encoded
                    b64_video = video
                else:
                    b64_video = b2a_base64(video).decode('ascii').rstrip()

                output = """<video {0} {1} {2}>
 <source src="data:{3};base64,{4}" type="{3}">
 Your browser does not support the video tag.
 </video>""".format(self.html_attributes, width, height, mimetype, b64_video)
                return output

        filepath = mkstemp(suffix=".ogv")[1]
        self.save(filepath, nb_loops=1, resolution=resolution, **kwargs)
        return CustomIPythonVideo(filepath, embed=True, width=resolution[0], html_attributes="controls loop autoplay")




