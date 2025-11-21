import numpy as np
import capytaine as cpt
from fakeblocks.bodies import join_bodies, cluster_bodies, minced

def test_cluster_bodies():
    centers = np.array([[0.0, 0.0, 0.0],
                        [10.0, 0.0, 0.0],
                        [11.0, 0.0, 0.0],
                        [10.0, 1.0, 0.0]])
    meshes = [cpt.mesh_sphere(center=c, name=str(c)) for c in centers]
    bodies = [cpt.FloatingBody(mesh=m) for m in meshes]
    joined_bodies = join_bodies(*bodies)
    clustered_bodies = cluster_bodies(*bodies)
    print(joined_bodies.mesh.tree_view())
    print(clustered_bodies.mesh.tree_view())
    assert isinstance(clustered_bodies, cpt.FloatingBody)
    assert isinstance(clustered_bodies.mesh, cpt.CollectionOfMeshes)
    assert clustered_bodies.mesh.merged() == joined_bodies.mesh.merged()
    assert meshes[0] in clustered_bodies.mesh  # The first body is at the top level independently from the other three


def test_mincing():
    body = cpt.FloatingBody(mesh=cpt.mesh_horizontal_cylinder(length=10, radius=0.5))
    body = body.minced((4, 1, 1))
    assert len(body.mesh) == 2
    assert np.all(body.mesh[0].faces_centers[:, 0] < 0)
    assert isinstance(body.mesh[0][0], cpt.Mesh)
    body = body.minced((1, 2, 2))
    assert isinstance(body.mesh[0][0][0][0], cpt.Mesh)
