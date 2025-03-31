import numpy as np
import xarray as xr
import capytaine as cpt

from capytaine.meshes.interface import MeshLike


def test_existing_classes():
    mesh = cpt.mesh_sphere()
    assert isinstance(mesh, MeshLike)
    assert isinstance(cpt.CollectionOfMeshes([mesh, mesh.translated_x(0.1)]), MeshLike)


def test_implementing_interface():
    class MyMesh:
        def __init__(self):
            self.vertices = np.array([[0.0, 0.0, -1.0], [1.0, 0.0, -1.0], [1.0, 1.0, -1.0], [0.0, 1.0, -1.0]])
            self.faces= np.array([[0, 1, 2, 3]])
            self.faces_centers = np.array([[0.5, 0.5, -1.0]])
            self.faces_areas = np.array([1.0])
            self.faces_radiuses = np.array([np.sqrt(2)/2])
            self.faces_normals = np.array([[0.0, 0.0, 1.0]])
            self.quadrature_points = (self.faces_centers.reshape(1, 1, 3), self.faces_areas)
            self.nb_vertices = 3
            self.nb_faces = 1

        def __short_str__(self):
            return str(self)

        def extract_faces(self, *args, **kwargs):
            return self

    assert isinstance(MyMesh(), MeshLike)

    mesh = MyMesh()
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    test_matrix = xr.Dataset(coords={"wavenumber":1.0, "radiating_dof": ["Heave"], "wave_direction": 1.0})
    solver = cpt.BEMSolver()
    solver.fill_dataset(test_matrix, body)
