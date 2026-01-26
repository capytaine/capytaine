"""Test related to the definition and use of geometric objects (planes, axes, ...)."""

import numpy as np

from capytaine.meshes.geometry import (
    e_x, e_y, e_z,
    Axis, Ox_axis, Oy_axis, Oz_axis,
    Plane, xOz_Plane, yOz_Plane, xOy_Plane,
    orthogonal_vectors, parallel_vectors,
    parallel_vectors_with_same_direction
)


def test_helper_functions():
    assert orthogonal_vectors(e_x, e_y)
    assert not orthogonal_vectors(e_x, 2*e_x)
    assert parallel_vectors(e_x, -e_x)
    assert not parallel_vectors(e_x, e_y)
    assert parallel_vectors_with_same_direction(e_y, 2*e_y)
    assert not parallel_vectors_with_same_direction(e_y, -e_y)


def test_axis():
    assert np.allclose(Axis(vector=(1, 1, 1)).vector, np.sqrt(3)/3 * np.ones((3,)))

    assert Axis(vector=(1, 1, 1), point=(0, 0, 0)) == Axis(vector=(1, 1, 1), point=(2, 2, 2))
    assert Axis(vector=(0, 0, 1), point=(0, 0, 0)) == Axis(vector=(2e-16, 3e-16, 1), point=(0,  0,  1.5))
    assert Axis(vector=(1, 1, 1), point=(0, 0, 0)) != Axis(vector=(1, 1, 1), point=(2, 2, 0))

    assert (2, 0, 0) in Ox_axis
    assert (0, 1, 0) not in Ox_axis

    assert Ox_axis.is_orthogonal_to(yOz_Plane)
    assert not Ox_axis.is_orthogonal_to((1, 1, 1))

    assert Ox_axis.angle_with_respect_to(Oy_axis) == np.pi/2
    assert Oy_axis.angle_with_respect_to(Ox_axis) == np.pi/2


def test_axis_transformation():
    assert Ox_axis.translated_x(10) == Ox_axis
    assert Ox_axis.translated_y(10) == Axis(vector=(1, 0, 0), point=(0, 10, 0))

    assert Ox_axis.rotated(Ox_axis, angle=np.pi/2) == Ox_axis
    assert Ox_axis.rotated(Oy_axis, angle=-np.pi/2) == Oz_axis

    assert Ox_axis.mirrored(plane=yOz_Plane) == Ox_axis
    assert Ox_axis.mirrored(plane=xOz_Plane.translated_y(2)) == Axis(vector=(1, 0, 0), point=(0, 4, 0))

    axis1 = Axis(vector=(1, 1, 1), point=(0, 0, 0))
    axis2 = Axis(vector=(1, 2, 3), point=(0, 0, 0))
    assert axis1.rotated_around_center_to_align_vectors(axis1.point, axis1.vector, axis2.vector) == axis2
    axis1.rotate(axis2, np.pi)
    assert axis1.rotated_around_center_to_align_vectors(axis1.point, axis1.vector, axis2.vector) == axis2
    axis1.vector *= -1
    assert axis1.rotated_around_center_to_align_vectors(axis1.point, axis1.vector, axis2.vector) == axis2

    axis1 = Axis(vector=(1, 1, 1), point=(1, 2, 0))
    axis2 = Axis(vector=(2, 2, 2), point=(0, 0, 0))
    assert axis1.translated_point_to_point(axis1.point, axis2.point) == axis2


def test_rotation_non_reference_axis():
    p = np.random.rand(3)
    axis = Axis(vector=(0, 0, 1), point=p)
    rotated_point = axis.rotate_points(p, np.pi)
    assert np.allclose(rotated_point, p)


def test_invariance_of_rotation_center():
    p = np.random.rand(3)
    axis1 = Axis(vector=(0, 0, 1), point=p)
    axis2 = Axis(vector=(0, 1, 0), point=p)
    assert np.allclose(axis1.rotated(axis2, np.pi/3).point, p)


def test_rotation_around_center_to_align_vectors_commutes_with_translation():
    n = np.random.rand(3)
    n /= np.linalg.norm(n)
    axis = Axis(vector=n, point=(0, 0, 0))
    a = axis.rotated_around_center_to_align_vectors((0, 0, 0), n, (0, 1, 0)).translated((1, 1, 1))
    b = axis.translated((1, 1, 1)).rotated_around_center_to_align_vectors((1, 1, 1), n, (0, 1, 0))
    assert np.allclose(a.vector, b.vector, (0, 1, 0))
    assert np.allclose(a.point, b.point)


def test_plane():
    assert (0, 1, 1) in yOz_Plane
    assert Oy_axis in yOz_Plane

    assert np.allclose(Plane(normal=(1, 1, 0)).normal, (np.sqrt(2)/2, np.sqrt(2)/2, 0))
    assert yOz_Plane == Plane(point=(0, 1, 1), normal=(2, 0, 0))

    assert xOy_Plane.is_orthogonal_to(Oz_axis)

    points_in_xplus = np.random.rand(10, 3) + np.array([1.0, -0.5, -0.5])
    assert np.all(yOz_Plane.distance_to_point(points_in_xplus) > 0)
    assert np.all(yOz_Plane.translated_x(-5.0).distance_to_point(points_in_xplus) > 0)
    assert not np.any(yOz_Plane.translated_x(5.0).distance_to_point(points_in_xplus) > 0)

    points_in_xminus = np.random.rand(10, 3) + np.array([-2.0, -0.5, -0.5])
    assert np.all(yOz_Plane.distance_to_point(points_in_xminus) < 0)
    assert not np.any(yOz_Plane.translated_x(-5.0).distance_to_point(points_in_xminus) < 0)
    assert np.all(yOz_Plane.translated_x(5.0).distance_to_point(points_in_xminus) < 0)


def test_plane_transformations():
    # TRANSLATIONS
    translated_plane = xOz_Plane.translate(vector=(1, 0, 0), inplace=False)
    assert xOz_Plane is not translated_plane
    assert xOz_Plane == translated_plane

    assert yOz_Plane.translated_x(10).rotated_y(np.pi/8).c == 10

    translated_plane = xOz_Plane.translate(vector=(0, 1, 0), inplace=False)
    assert translated_plane.c == 1
    assert np.all(translated_plane.normal == xOz_Plane.normal)

    # ROTATIONS
    rotated_plane = xOz_Plane.rotate(Oy_axis, angle=np.pi/12, inplace=False)
    assert rotated_plane == xOz_Plane.rotated(Oy_axis, angle=np.pi/12)
    assert xOz_Plane is not rotated_plane
    assert xOz_Plane == rotated_plane

    rotated_plane = xOz_Plane.rotate(Ox_axis, angle=np.pi/2, inplace=False)
    assert rotated_plane == xOy_Plane

    # MIRRORED BY ITSELF
    plane = Plane(normal=(1, 0, 0), point=(0.3, 0.2, 0.6))
    assert plane.mirrored(plane) != plane
    assert plane.mirrored(plane) == Plane(normal=(-1, 0, 0), point=(0.3, 0.2, 0.6))

    flipped_plane = plane.rotate(Axis(point=plane.point, vector=(0, 1, 0)), np.pi)
    assert flipped_plane == plane.mirror(plane)
