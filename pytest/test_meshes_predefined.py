
import pytest
import capytaine as cpt


def test_mesh_disk():
    from capytaine.meshes.predefined.cylinders import mesh_disk
    d = mesh_disk(resolution=(3, 5))
    assert isinstance(d, cpt.Mesh)
    assert d.nb_faces == 15

