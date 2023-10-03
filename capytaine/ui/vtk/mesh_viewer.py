"""3D display of a mesh with VTK.
Based on meshmagick <https://github.com/LHEEA/meshmagick> by François Rongère.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin, based on the work of François Rongère
# See LICENSE file at <https://github.com/mancellin/capytaine>

import datetime
from itertools import cycle
from os import getcwd

from capytaine.ui.vtk.helpers import compute_vtk_polydata
from capytaine.tools.optional_imports import import_optional_dependency

vtk = import_optional_dependency("vtk")

__year__ = datetime.datetime.now().year

COLORS = cycle([(1, 1, 0), (1, 0, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1)])


class MeshViewer:
    """This class implements a viewer based on VTK"""
    def __init__(self):

        # Building renderer
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(0.7706, 0.8165, 1.0)

        # Building render window
        self.render_window = vtk.vtkRenderWindow()
        self.render_window.SetSize(1024, 768)
        self.render_window.SetWindowName("Mesh viewer")
        self.render_window.AddRenderer(self.renderer)

        # Building interactor
        self.render_window_interactor = vtk.vtkRenderWindowInteractor()
        self.render_window_interactor.SetRenderWindow(self.render_window)
        self.render_window_interactor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()
        self.render_window_interactor.AddObserver('KeyPressEvent', self.on_key_press, 0.0)

        # Building axes view
        axes = vtk.vtkAxesActor()
        widget = vtk.vtkOrientationMarkerWidget()
        widget.SetOrientationMarker(axes)
        self.widget = widget

        self.widget.SetInteractor(self.render_window_interactor)
        self.widget.SetEnabled(1)
        self.widget.InteractiveOn()

        # Building command annotations
        command_text = ("left mouse : rotate\n"
                        "right mouse : zoom\n"
                        "middle mouse : pan\n"
                        "ctrl+left mouse : spin\n"
                        "n : (un)show normals\n"
                        "b : (un)show axes box\n"
                        "f : focus on the mouse cursor\n"
                        "r : reset view\n"
                        "s : surface representation\n"
                        "w : wire representation\n"
                        "h : (un)show Oxy plane\n"
                        "x : save\n"
                        "c : screenshot\n"
                        "q : quit")

        corner_annotation = vtk.vtkCornerAnnotation()
        corner_annotation.SetLinearFontScaleFactor(2)
        corner_annotation.SetNonlinearFontScaleFactor(1)
        corner_annotation.SetMaximumFontSize(20)
        corner_annotation.SetText(3, command_text)
        corner_annotation.GetTextProperty().SetColor(0., 0., 0.)
        self.renderer.AddViewProp(corner_annotation)

        copyright_text = (f"Capytaine — Copyright 2017-{__year__} University College Dublin\n"
                          f"based on Meshmagick Viewer — Copyright 2014-{__year__}, École Centrale de Nantes")

        copyright_annotation = vtk.vtkCornerAnnotation()
        copyright_annotation.SetLinearFontScaleFactor(0.5)
        copyright_annotation.SetNonlinearFontScaleFactor(1)
        copyright_annotation.SetMaximumFontSize(12)
        copyright_annotation.SetText(1, copyright_text)
        copyright_annotation.GetTextProperty().SetColor(0., 0., 0.)
        self.renderer.AddViewProp(copyright_annotation)

        self.polydatas = []
        self.normals = []
        self.axes = []
        self.oxy_plane = None

    # DISPLAY OF A MESH

    def add_polydata(self, polydata, color=(1, 1, 0), representation='surface'):
        """Add a polydata object to the viewer

        Parameters
        ----------
        polydata : vtkPolyData
            the object to be added
        color : array_like, optional
            the color of the object. Default is yellow (1, 1, 0)
        representation : str
            the representation mode of the object ('surface' or 'wireframe'). Default is 'surface'.
        """

        assert isinstance(polydata, vtk.vtkPolyData)
        assert representation in ('surface', 'wireframe')

        self.polydatas.append(polydata)

        # Building mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polydata)

        # Building actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # Properties setting
        actor.GetProperty().SetColor(color)
        actor.GetProperty().EdgeVisibilityOn()
        actor.GetProperty().SetEdgeColor(0, 0, 0)
        actor.GetProperty().SetLineWidth(1)
        actor.GetProperty().SetPointSize(10)
        if representation == 'wireframe':
            actor.GetProperty().SetRepresentationToWireframe()

        self.renderer.AddActor(actor)
        self.renderer.Modified()

    def add_mesh(self, mesh, color=None, representation='surface'):
        vtk_polydata = compute_vtk_polydata(mesh)
        if color is None:
            color = next(COLORS)
        self.add_polydata(vtk_polydata, color=color, representation=representation)

    # ADD MORE DETAILS

    def add_oxy_plane(self):
        """Displays the Oxy plane"""

        one_mesh_in_the_viewer = self.polydatas[0]

        plane = vtk.vtkPlaneSource()
        (xmin, xmax, ymin, ymax, _, _) = one_mesh_in_the_viewer.GetBounds()

        dx = max(0.1 * (xmax - xmin), 0.1)
        dy = max(0.1 * (ymax - ymin), 0.1)

        plane.SetOrigin(xmin - dx, ymax + dy, 0)
        plane.SetPoint1(xmin - dx, ymin - dy, 0)
        plane.SetPoint2(xmax + dx, ymax + dy, 0)
        plane.Update()
        polydata = plane.GetOutput()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polydata)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        color = [0., 102. / 255, 204. / 255]
        actor.GetProperty().SetColor(color)
        actor.GetProperty().SetEdgeColor(0, 0, 0)
        actor.GetProperty().SetLineWidth(1)

        self.renderer.AddActor(actor)
        self.renderer.Modified()
        self.oxy_plane = actor

    def show_normals(self):
        """Shows the normals of the current objects"""

        for polydata in self.polydatas:
            normals = vtk.vtkPolyDataNormals()
            normals.ConsistencyOff()
            normals.ComputeCellNormalsOn()
            normals.SetInputData(polydata)
            normals.Update()

            normals_at_centers = vtk.vtkCellCenters()
            normals_at_centers.SetInputConnection(normals.GetOutputPort())

            arrows = vtk.vtkArrowSource()
            arrows.SetTipResolution(16)
            arrows.SetTipLength(0.5)
            arrows.SetTipRadius(0.1)

            glyph = vtk.vtkGlyph3D()
            glyph.SetSourceConnection(arrows.GetOutputPort())
            glyph.SetInputConnection(normals_at_centers.GetOutputPort())
            glyph.SetVectorModeToUseNormal()
            glyph.SetScaleModeToScaleByVector()
            glyph.SetScaleFactor(1)  # FIXME: may be too big ...
            # glyph.SetVectorModeToUseNormal()
            # glyph.SetVectorModeToUseVector()
            # glyph.SetScaleModeToDataScalingOff()
            glyph.Update()

            glyph_mapper = vtk.vtkPolyDataMapper()
            glyph_mapper.SetInputConnection(glyph.GetOutputPort())

            glyph_actor = vtk.vtkActor()
            glyph_actor.SetMapper(glyph_mapper)

            self.renderer.AddActor(glyph_actor)
            self.normals.append(glyph_actor)

    def show_axes(self):
        """Shows the axes around the main object"""

        tprop = vtk.vtkTextProperty()
        tprop.SetColor(0., 0., 0.)
        tprop.ShadowOn()

        axes = vtk.vtkCubeAxesActor2D()
        axes.SetInputData(self.polydatas[0])

        axes.SetCamera(self.renderer.GetActiveCamera())
        axes.SetLabelFormat("%6.4g")
        axes.SetFlyModeToOuterEdges()
        axes.SetFontFactor(0.8)
        axes.SetAxisTitleTextProperty(tprop)
        axes.SetAxisLabelTextProperty(tprop)
        # axes.DrawGridLinesOn()

        self.renderer.AddViewProp(axes)
        self.axes.append(axes)

    def save(self):
        """Saves the main object in a 'mmviewer_save.vtp' vtp file is the current folder"""

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetDataModeToAscii()
        writer.SetFileName('mmviewer_save.vtp')

        for polydata in self.polydatas:
            writer.SetInputData(polydata)
        writer.Write()

        print("File 'mmviewer_save.vtp' written in %s" % getcwd())
        return

    def screenshot(self):
        """Saves a screenshot of the current window in a file screenshot.png"""
        w2if = vtk.vtkWindowToImageFilter()
        w2if.SetInput(self.render_window)
        w2if.Update()

        writer = vtk.vtkPNGWriter()
        writer.SetFileName("screenshot.png")
        writer.SetInputData(w2if.GetOutput())
        writer.Write()

        print("File 'screenshot.png' written in %s" % getcwd())
        return

    # INTERACTION

    def show(self):
        """Show the viewer"""
        self.renderer.ResetCamera()
        self.render_window.Render()
        self.render_window_interactor.Start()

    def on_key_press(self, obj, event):
        """Event trig at keystroke"""
        key = obj.GetKeySym()

        if key == 'n':
            if self.normals:
                # self.normals = False
                for actor in self.normals:
                    self.renderer.RemoveActor(actor)
                self.renderer.Render()
                self.normals = []
            else:
                self.show_normals()
                self.renderer.Render()

        elif key == 'b':
            if self.axes:
                for axis in self.axes:
                    self.renderer.RemoveActor(axis)
                self.axes = []
            else:
                self.show_axes()

        elif key == 'x':
            self.save()

        elif key == 'c':
            self.screenshot()

        elif key == 'h':
            if self.oxy_plane:
                self.renderer.RemoveActor(self.oxy_plane)
                self.oxy_plane = None
            else:
                self.add_oxy_plane()

        elif key == 'e' or key == 'q':
            self.render_window_interactor.GetRenderWindow().Finalize()
            self.render_window_interactor.TerminateApp()

    def finalize(self):
        """Cleanly close the viewer"""
        del self.render_window
        del self.render_window_interactor

    # =======================================
    # =======================================
    # =======================================

    # OTHER METHODS THAT ARE CURRENTLY UNUSED

    def add_point(self, pos, color=(0, 0, 0)):
        """Add a point to the viewer

        Parameters
        ----------
        pos : array_like
            The point's position
        color : array_like, optional
            The RGB color required for the point. Default is (0, 0, 0) corresponding to black.

        Returns
        -------
        vtkPolyData
        """

        assert len(pos) == 3

        p = vtk.vtkPoints()
        v = vtk.vtkCellArray()

        i = p.InsertNextPoint(pos)
        v.InsertNextCell(1)
        v.InsertCellPoint(i)

        pd = vtk.vtkPolyData()
        pd.SetPoints(p)
        pd.SetVerts(v)

        self.add_polydata(pd, color=color)

        return pd

    def add_line(self, p0, p1, color=(0, 0, 0)):
        """Add a line to the viewer

        Parameters
        ----------
        p0 : array_like
            position of one end point of the line
        p1 : array_like
            position of a second end point of the line
        color : array_like, optional
            RGB color of the line. Default is black (0, 0, 0)

        Returns
        -------
        vtkPolyData
        """

        assert len(p0) == 3 and len(p1) == 3

        points = vtk.vtkPoints()
        points.InsertNextPoint(p0)
        points.InsertNextPoint(p1)

        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, 0)
        line.GetPointIds().SetId(1, 1)

        lines = vtk.vtkCellArray()
        lines.InsertNextCell(line)

        lines_pd = vtk.vtkPolyData()
        lines_pd.SetPoints(points)
        lines_pd.SetLines(lines)

        self.add_polydata(lines_pd, color=color)

        return lines_pd

    def add_vector(self, point, value, scale=1, color=(0, 0, 0)):
        """Add a vector to the viewer

        Parameters
        ----------
        point : array_like
            starting point position of the vector
        value : float
            the magnitude of the vector
        scale : float, optional
            the scaling to apply to the vector for better visualization. Default is 1.
        color : array_like
            The color of the vector. Default is black (0, 0, 0)
        """

        points = vtk.vtkPoints()
        idx = points.InsertNextPoint(point)

        vert = vtk.vtkCellArray()
        vert.InsertNextCell(1)
        vert.InsertCellPoint(idx)
        pd_point = vtk.vtkPolyData()
        pd_point.SetPoints(points)
        pd_point.SetVerts(vert)

        arrow = vtk.vtkArrowSource()
        arrow.SetTipResolution(16)
        arrow.SetTipLength(0.1)
        arrow.SetTipRadius(0.02)
        arrow.SetShaftRadius(0.005)

        vec = vtk.vtkFloatArray()
        vec.SetNumberOfComponents(3)
        v0, v1, v2 = value / scale
        vec.InsertTuple3(idx, v0, v1, v2)
        pd_point.GetPointData().SetVectors(vec)

        g_glyph = vtk.vtkGlyph3D()
        # g_glyph.SetScaleModeToDataScalingOff()
        g_glyph.SetVectorModeToUseVector()
        g_glyph.SetInputData(pd_point)
        g_glyph.SetSourceConnection(arrow.GetOutputPort())
        g_glyph.SetScaleModeToScaleByVector()
        # g_glyph.SetScaleFactor(10)
        g_glyph.ScalingOn()
        g_glyph.Update()

        g_glyph_mapper = vtk.vtkPolyDataMapper()
        g_glyph_mapper.SetInputConnection(g_glyph.GetOutputPort())

        g_glyph_actor = vtk.vtkActor()
        g_glyph_actor.SetMapper(g_glyph_mapper)
        g_glyph_actor.GetProperty().SetColor(color)

        self.renderer.AddActor(g_glyph_actor)

    def add_plane(self, center, normal):
        """Add a plane to the viewer

        Parameters
        ----------
        center : array_like
            The origin of the plane
        normal : array_like
            The normal of the plane
        """

        plane = vtk.vtkPlaneSource()
        plane.SetCenter(center)
        plane.SetNormal(normal)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(plane.GetOutput())

        # FIXME: terminer l'implementation et l'utiliser pour le plan de la surface libre
