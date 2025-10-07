import numpy as np
import xarray as xr
import capytaine as cpt

from capytaine.meshes.mesh_like_protocol import MeshLike


def test_existing_classes_are_meshlike():
    mesh = cpt.mesh_sphere()
    assert isinstance(mesh, MeshLike)
    assert isinstance(cpt.CollectionOfMeshes([mesh, mesh.translated_x(0.1)]), MeshLike)

class MyMesh:
    """A set of horizontal panels at a list of depths"""
    def __init__(self, zs=[-1.0]):
        self.zs = zs
        self.vertices = np.array(
            sum(
                ([[0.0, 0.0, z], [0.0, 1.0, z], [1.0, 1.0, z], [1.0, 0.0, z]]
                for z in zs),

                []
            )
        )
        self.faces = np.array([[0, 1, 2, 3] for z in zs])
        self.faces_centers = np.array([[0.5, 0.5, z] for z in zs])
        self.faces_areas = np.array([1.0 for z in zs])
        self.faces_radiuses = np.array([np.sqrt(2) / 2 for z in zs])
        self.faces_normals = np.array([[0.0, 0.0, -1.0] for z in zs])
        self.quadrature_points = (
            self.faces_centers.reshape(len(zs), 1, 3),
            self.faces_areas,
        )
        self.nb_vertices = 4*len(zs)
        self.nb_faces = len(zs)

    def __short_str__(self):
        return str(self)

    def extract_faces(self, *args, **kwargs):
        return self

    def join_meshes(*meshes, return_masks=False):
        joined = MyMesh(sum([m.zs for m in meshes], []))
        # Concatenate the lists of vertical positions
        if not return_masks:
            return joined
        else:
            masks = [np.full((joined.nb_faces,), False) for _ in meshes]
            accumulate_shifts = 0
            for i, m in enumerate(meshes):
                masks[i][accumulate_shifts:accumulate_shifts+m.nb_faces] = True
                accumulate_shifts = accumulate_shifts + m.nb_faces
            return joined, masks

    def with_normal_vector_going_down(self, **kwargs):
        return self

def test_minimal_implementation_is_meshlike():
    assert isinstance(MyMesh(), MeshLike)

def test_minimal_implementation_solve_pb():
    mesh = MyMesh(zs=[-1.0])
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb = cpt.RadiationProblem(body=body, omega=1.0, radiating_dof="Heave")
    solver.solve(pb, method="indirect")

def test_minimal_implementation_compute_velocity():
    mesh = MyMesh(zs=[-1.0])
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb = cpt.RadiationProblem(body=body, omega=1.0, radiating_dof="Heave")
    res = solver.solve(pb, method="indirect")
    solver.compute_velocity(mesh, res)

def test_minimal_implementation_solve_pb_with_lid():
    mesh = MyMesh(zs=[-1.0])
    lid_mesh = MyMesh(zs=[0.0])
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=lid_mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb = cpt.RadiationProblem(body=body, omega=1.0, radiating_dof="Heave")
    solver.solve(pb, method="indirect")

def test_minimal_implementation_fill_dataset():
    mesh = MyMesh(zs=[-1.0])
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    test_matrix = xr.Dataset(
        coords={"wavenumber": 1.0, "radiating_dof": ["Heave"], "wave_direction": 1.0}
    )
    solver.fill_dataset(test_matrix, body, method="direct")
