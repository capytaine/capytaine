import logging

import numpy as np
import pytest
import capytaine as cpt

from capytaine.meshes.meshes import Mesh
from capytaine.meshes.symmetric_meshes import ReflectionSymmetricMesh, RotationSymmetricMesh

logging.basicConfig(level=logging.DEBUG)

# DISK

def test_mesh_disk():
    from capytaine.meshes.predefined.cylinders import mesh_disk
    d = mesh_disk(resolution=(3, 6), name="foo")
    assert isinstance(d, Mesh)
    assert d.nb_faces == 18
    assert d.name == "foo"

def test_mesh_disk_position():
    from capytaine.meshes.predefined.cylinders import mesh_disk
    d = mesh_disk(center=(1.0, 0.5, -0.5), normal=(0, 1, 0))
    assert np.allclose(d.vertices[:, 1], 0.5)

def test_mesh_disk_reflection_symmetry():
    from capytaine.meshes.predefined.cylinders import mesh_disk
    d = mesh_disk(resolution=(3, 6), reflection_symmetry=True,)
    assert isinstance(d, ReflectionSymmetricMesh)
    assert np.all(d.half.vertices[:, 1] >= 0.0)

def test_mesh_disk_axial_symmetry():
    from capytaine.meshes.predefined.cylinders import mesh_disk
    d = mesh_disk(resolution=(3, 6), axial_symmetry=True)
    assert isinstance(d, RotationSymmetricMesh)

def test_mesh_disk_both_symmetries():
    from capytaine.meshes.predefined.cylinders import mesh_disk
    with pytest.raises(NotImplementedError):
        mesh_disk(axial_symmetry=True, reflection_symmetry=True)

@pytest.mark.parametrize("max_rad", [0.5, 1.0, 5.0])
def test_mesh_disk_max_size(max_rad):
    from capytaine.meshes.predefined.cylinders import mesh_disk
    d = mesh_disk(radius=10.0, faces_max_radius=max_rad)
    assert d.faces_radiuses.max() <= max_rad

# VERTICAL CYLINDER

def test_mesh_vertical_cylinder():
    from capytaine.meshes.predefined.cylinders import mesh_vertical_cylinder
    v = mesh_vertical_cylinder(resolution=(2, 8, 10), name="foo")
    assert isinstance(v, Mesh)
    assert v.nb_faces == (10+2*2)*8
    assert v.name == "foo"

def test_mesh_vertical_cylinder_position():
    from capytaine.meshes.predefined.cylinders import mesh_vertical_cylinder
    v = mesh_vertical_cylinder(
        center=(-10.0, 5.0, 1.0),
        length=6.0,
        radius=4.0,
        resolution=(0, 100, 5)
    )
    assert np.isclose(np.min(v.vertices[:, 0]), -14.0)
    assert np.isclose(np.min(v.vertices[:, 1]), 1.0)
    assert np.isclose(np.min(v.vertices[:, 2]), -2.0)

def test_mesh_vertical_cylinder_reflection_symmetry():
    from capytaine.meshes.predefined.cylinders import mesh_vertical_cylinder
    v = mesh_vertical_cylinder(reflection_symmetry=True, center=(0.0, 0.0, 6.0))
    assert isinstance(v, ReflectionSymmetricMesh)
    assert v.merged().nb_faces == mesh_vertical_cylinder(reflection_symmetry=False).nb_faces
    assert np.all(v.half.vertices[:, 0] >= 0.0)
    assert not np.all(v.vertices[:, 0] >= 0.0)
    assert np.all(v.vertices[:, 2] > 0.0)

def test_mesh_vertical_cylinder_axial_symmetry():
    from capytaine.meshes.predefined.cylinders import mesh_vertical_cylinder
    v = mesh_vertical_cylinder(axial_symmetry=True)
    assert isinstance(v, RotationSymmetricMesh)

def test_mesh_vertical_cylinder_both_symmetries():
    from capytaine.meshes.predefined.cylinders import mesh_vertical_cylinder
    with pytest.raises(NotImplementedError):
        mesh_vertical_cylinder(axial_symmetry=True, reflection_symmetry=True)

@pytest.mark.parametrize("max_rad", [0.5, 1.0, 5.0])
def test_mesh_vertical_max_size(max_rad):
    from capytaine.meshes.predefined.cylinders import mesh_vertical_cylinder
    d = mesh_vertical_cylinder(radius=10.0, faces_max_radius=max_rad)
    assert d.faces_radiuses.max() <= max_rad

# HORIZONTAL CYLINDER

def test_mesh_horizontal_cylinder():
    from capytaine.meshes.predefined.cylinders import mesh_horizontal_cylinder
    h = mesh_horizontal_cylinder(resolution=(2, 8, 10), name="foo")
    assert isinstance(h, Mesh)
    assert h.nb_faces == (10+2*2)*8
    assert h.name == "foo"

def test_mesh_horizontal_cylinder_position():
    from capytaine.meshes.predefined.cylinders import mesh_horizontal_cylinder
    h = mesh_horizontal_cylinder(length=10.0, radius=1.0, resolution=(0, 100, 3), center=(10.0, 1.0, 0.0))
    assert np.isclose(np.min(h.vertices[:, 0]), 5.0)
    assert np.isclose(np.min(h.vertices[:, 1]), 0.0)

def test_mesh_horizontal_cylinder_reflection_symmetry():
    from capytaine.meshes.predefined.cylinders import mesh_horizontal_cylinder
    h = mesh_horizontal_cylinder(reflection_symmetry=True, center=(0.0, 0.0, 5.0))
    assert isinstance(h, ReflectionSymmetricMesh)
    assert h.merged().nb_faces == mesh_horizontal_cylinder(reflection_symmetry=False).nb_faces
    assert not np.all(h.vertices[:, 1] >= 0.0)  # All panels are on the same side
    assert np.all(h.half.vertices[:, 1] >= 0.0)  # All panels are on the same side
    assert np.all(h.vertices[:, 2] > 0.0)

@pytest.mark.parametrize("max_rad", [0.5, 1.0, 5.0])
def test_mesh_horizontal_max_size(max_rad):
    from capytaine.meshes.predefined.cylinders import mesh_horizontal_cylinder
    d = mesh_horizontal_cylinder(radius=10.0, faces_max_radius=max_rad)
    assert d.faces_radiuses.max() <= max_rad

# SPHERE

def test_mesh_sphere():
    from capytaine.meshes.predefined.spheres import mesh_sphere
    d = mesh_sphere(resolution=(3, 6), name="foo")
    assert isinstance(d, Mesh)
    assert d.nb_faces == 18
    assert d.name == "foo"

def test_mesh_sphere_position():
    from capytaine.meshes.predefined.spheres import mesh_sphere
    d = mesh_sphere(radius=0.1, center=(1.0, 0.5, -0.5))
    assert np.all(d.vertices[:, 0] <= 1.1)
    assert np.all(0.9 <= d.vertices[:, 0])
    assert np.all(d.vertices[:, 2] <= -0.4)
    assert np.all(-0.6 <= d.vertices[:, 2])

def test_mesh_sphere_axial_symmetry():
    from capytaine.meshes.predefined.spheres import mesh_sphere
    d = mesh_sphere(resolution=(3, 6), axial_symmetry=True)
    assert isinstance(d, RotationSymmetricMesh)
    assert d.axis == 'z+'

@pytest.mark.parametrize("max_rad", [0.5, 1.0, 5.0])
def test_mesh_sphere_max_size(max_rad):
    from capytaine.meshes.predefined.spheres import mesh_sphere
    d = mesh_sphere(radius=10.0, faces_max_radius=max_rad)
    assert d.faces_radiuses.max() <= max_rad
    d = mesh_sphere(
            radius=10.0, resolution=(10, 10),
            faces_max_radius=max_rad)
    assert d.faces_radiuses.max() <= max_rad

# RECTANGLES

def test_invalid_input():
    from capytaine.meshes.predefined.rectangles import mesh_rectangle
    with pytest.raises(AssertionError):
        mesh_rectangle(size=(-4, 5))
    with pytest.raises(AssertionError):
        mesh_rectangle(resolution=(4.1, 5))
    with pytest.raises(ValueError):
        mesh_rectangle(resolution=(3, 3), reflection_symmetry=True)

def test_mesh_rectangle():
    from capytaine.meshes.predefined.rectangles import mesh_rectangle
    d = mesh_rectangle(resolution=(3, 6), name="foo")
    assert isinstance(d, Mesh)
    assert d.nb_faces == 18
    assert d.name == "foo"

def test_mesh_rectangle_position():
    from capytaine.meshes.predefined.rectangles import mesh_rectangle
    d = mesh_rectangle(center=(1.0, 0.5, -0.5), normal=(0, 1, 0))
    assert np.allclose(d.vertices[:, 1], 0.5)

def test_mesh_rectangle_reflection_symmetry():
    from capytaine.meshes.predefined.rectangles import mesh_rectangle
    d = mesh_rectangle(resolution=(4, 6), reflection_symmetry=True, normal=(0, 0, -1))
    assert isinstance(d, ReflectionSymmetricMesh)
    assert np.allclose(d.faces_normals, (0, 0, -1))

def test_mesh_rectangle_invalid_reflection_symmetry():
    from capytaine.meshes.predefined.rectangles import mesh_rectangle
    with pytest.raises(ValueError):
        mesh_rectangle(resolution=(4, 6), reflection_symmetry=True, normal=(1, 1, 1))

@pytest.mark.parametrize("max_rad", [0.5, 1.0, 5.0])
def test_mesh_rectangle_max_size(max_rad):
    from capytaine.meshes.predefined.rectangles import mesh_rectangle
    d = mesh_rectangle(size=(10.0, 10.0), faces_max_radius=max_rad)
    assert d.faces_radiuses.max() <= max_rad

# PARALLELEPIPED

def test_mesh_parallelepiped():
    from capytaine.meshes.predefined.rectangles import mesh_parallelepiped
    p = mesh_parallelepiped(resolution=(3, 6, 3), name="foo")
    assert isinstance(p, Mesh)
    assert p.nb_faces == 90
    assert p.name == "foo"

def test_mesh_parallelepiped_sides():
    from capytaine.meshes.predefined.rectangles import mesh_parallelepiped
    p = mesh_parallelepiped(size=(5.0, 5.0, 5.0), missing_sides=["top", "bottom"])
    assert p.nb_faces == 64
    assert not np.any((abs(p.faces_centers[:, 0]) < 2.0) & (abs(p.faces_centers[:, 1]) < 2.0))
    # That is: no face with center (x, y, z) such that -2 < x < 2 and -2 < y < 2

def test_mesh_parallelepiped_position():
    from capytaine.meshes.predefined.rectangles import mesh_parallelepiped
    p = mesh_parallelepiped(center=(2.0, -2.0, 1.0))
    assert np.all(p.faces_centers[:, 0] <= 3.0)
    assert np.all(1.0 <= p.faces_centers[:, 0])

def test_mesh_parallelepiped_reflection_symmetry():
    from capytaine.meshes.predefined.rectangles import mesh_parallelepiped
    p = mesh_parallelepiped(resolution=(4, 4, 4), reflection_symmetry=True)
    assert isinstance(p, ReflectionSymmetricMesh)
    assert isinstance(p.half, ReflectionSymmetricMesh)
    assert p.nb_faces == 6*16

@pytest.mark.parametrize("max_rad", [0.5, 1.0, 5.0])
def test_mesh_parallelepiped_max_size(max_rad):
    from capytaine.meshes.predefined.rectangles import mesh_parallelepiped
    d = mesh_parallelepiped(size=(10.0, 10.0, 10.0), faces_max_radius=max_rad)
    assert d.faces_radiuses.max() <= max_rad
    d = mesh_parallelepiped(size=(1.0, 10.0, 1.0), faces_max_radius=max_rad)
    assert d.faces_radiuses.max() <= max_rad
