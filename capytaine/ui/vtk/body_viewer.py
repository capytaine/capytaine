"""3D display of floating body with VTK."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

from capytaine.tools.optional_imports import import_optional_dependency
from capytaine.ui.vtk.mesh_viewer import MeshViewer

vtk = import_optional_dependency("vtk")

class FloatingBodyViewer(MeshViewer):

    def __init__(self):
        super().__init__()
        self.dofs_data = {}

    def add_body(self, body, **kwargs):
        self.add_mesh(body.mesh, **kwargs)
        if body.lid_mesh is not None:
            self.add_mesh(body.lid_mesh, **kwargs)

        for dof in body.dofs:
            vtk_data_array = vtk.vtkFloatArray()
            vtk_data_array.SetNumberOfComponents(3)
            vtk_data_array.SetNumberOfTuples(body.mesh.nb_faces)
            for i, vector in enumerate(body.dofs[dof]):
                vtk_data_array.SetTuple3(i, *vector)
            self.dofs_data[dof] = vtk_data_array
            # vtk_polydata.GetCellData().SetVectors(vtk_data_array)
