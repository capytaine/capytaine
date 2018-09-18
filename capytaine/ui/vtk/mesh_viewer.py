#!/usr/bin/env python
# coding: utf-8
"""VTK display of the floating body with its degrees of freedom."""
# This file is part of "capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import vtk

from capytaine.ui.vtk.MMviewer import MMViewer


class FloatingBodyViewer(MMViewer):
    def show_normals(self):
        """Shows the normals of the current objects"""

        for polydata in self.polydatas:
            normals = vtk.vtkPolyDataNormals()
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
            glyph.SetVectorModeToUseVector()
            glyph.SetScaleModeToScaleByVector()
            glyph.SetScaleFactor(1)  # FIXME: may be too big ...
            glyph.Update()

            glyph_mapper = vtk.vtkPolyDataMapper()
            glyph_mapper.SetInputConnection(glyph.GetOutputPort())

            glyph_actor = vtk.vtkActor()
            glyph_actor.SetMapper(glyph_mapper)

            self.renderer.AddActor(glyph_actor)
            self.normals.append(glyph_actor)
