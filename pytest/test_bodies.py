#!/usr/bin/env python
# coding: utf-8

import pytest

import numpy as np
import xarray as xr

from capytaine import BEMSolver
from capytaine.bodies import FloatingBody
from capytaine.meshes.meshes import Mesh
from capytaine.meshes.geometry import Axis, Plane
from capytaine.bodies.predefined.spheres import Sphere
from capytaine.bodies.predefined.cylinders import HorizontalCylinder


def test_dof():
    nodes = np.array([[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]])
    faces = np.array([[0, 1, 2, 3]])
    body = FloatingBody(Mesh(nodes, faces), name="one_face")
    assert body.dofs == {}

    body.add_translation_dof(direction=(1.0, 0.0, 0.0), name="1")
    assert np.allclose(body.dofs["1"], np.array([1.0, 0.0, 0.0]))

    body.add_translation_dof(direction=(0.0, 1.0, 0.0), name="2")
    assert np.allclose(body.dofs["2"], np.array([0.0, 1.0, 0.0]))

    body.add_rotation_dof(Axis(vector=(0.0, 0.0, 1.0)), name="3")
    body.add_rotation_dof(Axis(point=(0.5, 0, 0), vector=(0.0, 0.0, 1.0)), name="4")


def test_dof_name_inference():
    body = HorizontalCylinder()
    body.add_translation_dof(direction=(1, 0, 0), name="Surge_1")
    for dofname in ['Surge', 'SURGE', 'surge']:
        body.add_translation_dof(name=dofname)
        assert np.allclose(body.dofs[dofname], body.dofs['Surge_1'])

    body.add_rotation_dof(name="Pitch")
    body.add_rotation_dof(name="yaw")

    body.dofs.clear()
    body.add_all_rigid_body_dofs()


def test_cropping_body_with_manual_dof():
    # https://github.com/capytaine/capytaine/issues/204
    sphere = Sphere()
    sphere.dofs["Surge"] = [(1, 0, 0) for face in sphere.mesh.faces]
    sphere.keep_immersed_part()


def test_bodies():
    body = Sphere(name="sphere", axial_symmetry=False)
    assert str(body) == "sphere"
    repr(body)
    assert np.allclose(body.geometric_center, (0, 0, 0))
    body.add_translation_dof(name="Surge")
    body.add_translation_dof(name="Heave")

    # Extract faces
    body.extract_faces(np.where(body.mesh.faces_centers[:, 2] < 0)[0])

    # Clipping
    body.keep_immersed_part(inplace=False)

    # Mirror of the dofs
    mirrored = body.mirrored(Plane(point=(1, 0, 0), normal=(1, 0, 0)))
    assert np.allclose(mirrored.geometric_center, np.array([2, 0, 0]))
    assert np.allclose(body.dofs['Surge'], -mirrored.dofs['Surge'])

    # Rotation of the dofs
    sideways = body.rotated(Axis(point=(0, 0, 0), vector=(0, 1, 0)), np.pi/2)
    assert np.allclose(sideways.dofs['Heave'][0], np.array([1, 0, 0]))

    upside_down = body.rotated(Axis(point=(0, 0, 0), vector=(0, 1, 0)), np.pi)
    assert np.allclose(body.dofs['Heave'], -upside_down.dofs['Heave'])

    # Copy of the body
    copy_of_body = body.copy(name="copy_of_sphere")
    copy_of_body.translate_x(10.0)
    copy_of_body.add_translation_dof(name="Heave")

    # Join bodies
    both = body.join_bodies(copy_of_body)
    assert set(both.dofs) == {'sphere__Surge', 'copy_of_sphere__Surge', 'sphere__Heave', 'copy_of_sphere__Heave'}


@pytest.mark.parametrize("z_center", [0, 2, -2])
@pytest.mark.parametrize("as_collection_of_meshes", [True, False])
def test_clipping_of_dofs(z_center, as_collection_of_meshes):
    """Check that clipping a body with a dof is the same as clipping the body ant then adding the dof."""
    full_sphere = Sphere(center=(0, 0, z_center), name="sphere", axial_symmetry=as_collection_of_meshes, clip_free_surface=False)
    axis = Axis(point=(1, 0, 0), vector=(1, 0, 0))

    full_sphere.add_rotation_dof(axis, name="test_dof")
    clipped_sphere = full_sphere.keep_immersed_part(free_surface=0.0, sea_bottom=-np.infty, inplace=False)

    other_clipped_sphere = FloatingBody(mesh=clipped_sphere.mesh, name="other_sphere")
    other_clipped_sphere.add_rotation_dof(axis, name="test_dof")

    if clipped_sphere.mesh.nb_faces > 0:
        assert np.allclose(clipped_sphere.dofs['test_dof'], other_clipped_sphere.dofs['test_dof'])
    else:
        assert len(clipped_sphere.dofs['test_dof']) == 0


def test_mincing():
    body = HorizontalCylinder(length=10, radius=0.5, translation_symmetry=False, reflection_symmetry=False)
    body = body.minced((4, 1, 1))
    assert len(body.mesh) == 2
    assert np.all(body.mesh[0].faces_centers[:, 0] < 0)
    assert isinstance(body.mesh[0][0], Mesh)
    body = body.minced((1, 2, 2))
    assert isinstance(body.mesh[0][0][0][0], Mesh)


def test_assemble_regular_array():
    body = Sphere()
    body.add_all_rigid_body_dofs()
    array = body.assemble_regular_array(distance=2.0, nb_bodies=(2, 3))
    assert array.mesh.nb_faces == 6*body.mesh.nb_faces

    # Check name and order of the dofs
    assert list(array.dofs.keys())[0:3] == ["0_0__Surge", "0_0__Sway", "0_0__Heave"]
    assert "2_1__Heave" not in array.dofs.keys()

    # Check that the dofs coresponds to the right panels
    faces_1_0 = np.where(array.dofs["1_0__Heave"] != 0.0)[0]
    fc_1_0 = array.mesh.merged().faces_centers[faces_1_0, :]
    assert np.all(1.0 <= fc_1_0[:, 0]) and np.all(fc_1_0[:, 0] <= 3.0)  #   1 < x < 3
    assert np.all(-1.0 <= fc_1_0[:, 1]) and np.all(fc_1_0[:, 1] <= 1.0) #  -1 < y < 1


r = 1.0
locations = np.array([[1,0],[-1,0],[0,1]])*r*2
n_bodies = locations.shape[0]

@pytest.fixture
def fb_array():
    
    sphere = Sphere(
                    radius=r,            # Dimension
                    center=(0, 0, 0),    # Position
                    nphi=4, ntheta=10,   # Fineness of the mesh
    )
    my_axis = Axis((0, 1, 0), 
                   point=(0,0,0))
    sphere.add_rotation_dof(axis=my_axis)
    sphere.keep_immersed_part()
    
    return sphere.assemble_arbitrary_array(locations)


def test_consistent_dofs_to_faces(fb_array):
    num_active_faces = []
    for fb_dof in fb_array.dofs.items():
        num_active_faces.append(np.count_nonzero(np.count_nonzero(fb_dof[1],axis=1)))

    ma = np.array(num_active_faces)

    tot_faces = fb_array.mesh.nb_faces
    exp_active_faces_per_dof = int(tot_faces/n_bodies)

    assert np.all(ma==exp_active_faces_per_dof)

def test_solve_hydrodynamics(fb_array):
    solver = BEMSolver()
    test_matrix = xr.Dataset(coords={
          'rho': 1e3,                         
          'water_depth': [np.infty],          
          'omega': np.pi * 2 / 1,
          'wave_direction': 0,
          'radiating_dof': list(fb_array.dofs.keys()),
          })
    data = solver.fill_dataset(test_matrix, [fb_array],
                                 mesh=True,
                                 wavelength=True,
                                 wavenumber=True)
    assert data.influenced_dof.size == n_bodies
    assert data.radiating_dof.size == n_bodies
    assert data.added_mass.notnull().all()
    assert data.radiation_damping.notnull().all()
    assert data.diffraction_force.notnull().all()
    assert data.Froude_Krylov_force.notnull().all()
