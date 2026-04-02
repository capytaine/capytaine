from functools import lru_cache
import logging
import json
from pathlib import Path

import pytest

import numpy as np
import capytaine as cpt

from capytaine.meshes.predefined import mesh_sphere, mesh_horizontal_cylinder, mesh_vertical_cylinder

######################################################################

################################
# Immersed part and waterplane #
################################

@lru_cache
def immersed_sphere():
    return mesh_sphere(
        radius=1.0,
        center=(0, 0, -2),
        resolution=(30, 30)
    )

def test_wet_surface_area_of_submerged_sphere():
    assert np.isclose(
        immersed_sphere().wet_surface_area,
        4*np.pi*1**2,
        rtol=5e-2
    )

def test_volume_of_submerged_sphere():
    assert np.isclose(
        immersed_sphere().volume,
        4/3*np.pi*1**3,
        rtol=5e-2
    )

def test_disp_mass_of_submerged_sphere():
    analytical_volume = 4/3*np.pi*1.0**3
    assert np.isclose(
        immersed_sphere().disp_mass(rho=1000),
        1000*analytical_volume,
        rtol=2e-2
    )

def test_waterplane_area_of_submerged_sphere():
    assert np.isclose(
        immersed_sphere().waterplane_area,
        0.0
    )

def test_center_of_buoyancy_of_submerged_sphere():
    assert np.allclose(
            immersed_sphere().center_of_buoyancy,
            [0, 0, -2]
            )

def test_waterplane_center_of_submerged_sphere():
    assert immersed_sphere().waterplane_center is None

######################################################################

@lru_cache
def floating_sphere():
    return mesh_sphere(
        radius=1.0,
        center=(0, 0, 0),
        resolution=(30, 30)
    )

def test_wet_surface_area_of_floating_sphere():
    assert np.isclose(
        floating_sphere().wet_surface_area,
        4*np.pi*1**2/2,
        rtol=5e-2
    )

def test_volume_of_floating_sphere():
    assert np.isclose(
        floating_sphere().volume,
        4/3*np.pi*1**3,
        rtol=5e-2
    )

def test_waterplane_area_of_floating_sphere():
    assert np.allclose(
        floating_sphere().waterplane_area,
        np.pi*1**2,
        atol=5e-2
    )

def test_disp_mass_of_floating_sphere():
    analytical_volume = 4/3*np.pi*1.0**3/2
    assert np.isclose(
        floating_sphere().disp_mass(rho=1000),
        1000*analytical_volume,
        rtol=2e-2
    )

def test_center_of_buoyancy_of_floating_sphere():
    cob = floating_sphere().center_of_buoyancy
    assert np.allclose(cob[:2], [0, 0])
    assert cob[2] < -1e-2  # Strictly below zero

def test_waterplane_center_of_floating_sphere():
    assert np.allclose(
        floating_sphere().waterplane_center,
        [0.0, 0.0]
    )

#############
# Stiffness #
#############

def test_stiffness_when_no_dofs():
    body = cpt.FloatingBody(mesh=floating_sphere())
    with pytest.raises(AttributeError, match=".* no dof .*"):
        body.compute_hydrostatic_stiffness()

def test_stiffness_no_center_of_mass():
    body = cpt.FloatingBody(
        mesh=floating_sphere(),
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.3))
    )
    with pytest.raises(ValueError, match=".*no center of mass.*"):
        body.compute_hydrostatic_stiffness()

@lru_cache
def rigid_body():
    rigid_body = cpt.FloatingBody(
        mesh=floating_sphere(),
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.3)),
        center_of_mass=(0, 0, -0.3)
    )
    return rigid_body

@lru_cache
def custom_dof_body():
    mesh = floating_sphere()
    dof = np.array([(0, 0, z) for (x, y, z) in mesh.faces_centers])
    custom_dof_body = cpt.FloatingBody(
        mesh=mesh,
        dofs={"elongate_in_z": dof},
        center_of_mass=mesh.center_of_buoyancy,
    )
    return custom_dof_body

def test_stiffness_dof_ordering():
    K = rigid_body().compute_hydrostatic_stiffness()
    assert np.all(K.coords["radiating_dof"].values == np.array(['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']))

@pytest.mark.parametrize("body", [rigid_body, custom_dof_body])
def test_stiffness_invariance_by_clipping_at_free_surface(body):
    K1 = body().compute_hydrostatic_stiffness()
    K2 = body().immersed_part().compute_hydrostatic_stiffness()
    assert np.allclose(K1, K2)

@pytest.mark.parametrize("body", [rigid_body, custom_dof_body])
def test_stiffness_invariance_by_translation(body):
    K1 = body().compute_hydrostatic_stiffness()
    K2 = body().translated([1.0, 0.0, 0.0]).compute_hydrostatic_stiffness()
    assert np.allclose(K1, K2)

def test_stiffness_single_rotation_dof():
    # Was actually not possible in version 2.x
    mesh = cpt.mesh_parallelepiped()
    body_1 = cpt.FloatingBody(
        mesh=mesh,
        dofs=cpt.rigid_body_dofs(only=["Pitch"], rotation_center=(0, 0, -2)),
        center_of_mass=(0, 0, 0)
    )
    K_1 = body_1.compute_hydrostatic_stiffness()
    assert K_1.shape == (1, 1)

    ref_body = cpt.FloatingBody(
        mesh=mesh,
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -2)),
        center_of_mass=(0, 0, 0)
    )
    full_K = ref_body.compute_hydrostatic_stiffness()
    assert np.isclose(K_1.values[0, 0], full_K.values[3, 3])

    body_2 = cpt.FloatingBody(mesh=mesh, center_of_mass=(0, 0, 0))
    body_2.add_rotation_dof(rotation_center=(0, 0, -2), direction=(0, 1, 0), name="pitch oh mon pitch")
    # Non standard name, but still a rigid body rotation
    K_2 = body_2.compute_hydrostatic_stiffness()
    assert K_2.shape == (1, 1)
    assert np.isclose(K_2.values, K_1.values)

    # Generalized dof approach despite standard name
    body_3 = cpt.FloatingBody(
        mesh=mesh,
        dofs={"Pitch": body_1.dofs["Pitch"].evaluate_motion(mesh)},
        center_of_mass=(0, 0, 0)
    )
    K_3 = body_3.compute_hydrostatic_stiffness()
    assert K_3.shape == (1, 1)
    assert not np.isclose(K_3.values, K_1.values)


# DIVERGENCE

def test_stiffness_elastic_dof_with_divergence():
    body = custom_dof_body().immersed_part()
    hs_1 = body.compute_hydrostatic_stiffness()
    hs_2 = body.compute_hydrostatic_stiffness(divergence={"elongate_in_z": np.ones(body.mesh.nb_faces)})
    assert hs_1.values[0, 0] != hs_2.values[0, 0]
    analytical_hs = - 1000.0 * 9.81 * (4 * body.volume * body.center_of_buoyancy[2])
    assert np.isclose(hs_2.values[0, 0], analytical_hs)

def test_stiffness_with_divergence_not_clipped():
    body = custom_dof_body()
    with pytest.raises(NotImplementedError):
        body.compute_hydrostatic_stiffness(divergence={"elongate_in_z": np.ones(body.mesh.nb_faces)})

def test_stiffness_with_malformed_divergence(caplog):
    body = custom_dof_body().immersed_part()
    hs_1 = body.compute_hydrostatic_stiffness()
    with caplog.at_level(logging.WARNING):
        hs_2 = body.compute_hydrostatic_stiffness(
            divergence={
                "foobar": np.ones(body.mesh.nb_faces)  # Wrong dof name ignored!
            }
        )
    assert hs_1.values[0, 0] == hs_2.values[0, 0]
    assert "without the divergence" in caplog.text

# MULTIBODY

# No caching here because we assign attributes
# (hydrostatic_stiffness and inertia_matrix) to the objects in the tests below
# @lru_cache
def several_bodies():
    a = cpt.FloatingBody(
            mesh=floating_sphere(),
            name='foo'
            )
    a.add_translation_dof(name="Heave")
    b = cpt.FloatingBody(
            mesh=floating_sphere().translated_x(3.0),
            dofs=cpt.rigid_body_dofs(),
            name='bar'
            )
    c = cpt.FloatingBody(
            mesh=floating_sphere().translated_y(-2.0),
            dofs=cpt.rigid_body_dofs(),
            name='baz'
            )
    return a, b, c

def test_joining_stiffness_when_joining_bodies():
    a, b, _ = several_bodies()
    a.hydrostatic_stiffness = a.add_dofs_labels_to_matrix(np.random.rand(1, 1))
    b.hydrostatic_stiffness = b.add_dofs_labels_to_matrix(np.random.rand(6, 6))
    both = a + b
    assert both.hydrostatic_stiffness.shape == (7, 7)

######################################################################

##################
# Inertia matrix #
##################

# Values

@pytest.mark.parametrize("mesh_and_immersed_volume", [
    (immersed_sphere(), 4/3*np.pi*1.0**3),
    (floating_sphere(), 2/3*np.pi*1.0**3),
], ids=['immersed', 'floating'])
@pytest.mark.parametrize("rho", [1000, 444])
def test_inertia_matrix_values_of_sphere(rho, mesh_and_immersed_volume):
    mesh, immersed_volume = mesh_and_immersed_volume
    center = np.mean(mesh.vertices, axis=0)
    sphere = cpt.FloatingBody(
        mesh=mesh,
        dofs=cpt.rigid_body_dofs(rotation_center=center),
        center_of_mass=center,
    )
    mass = rho*immersed_volume
    assert np.isclose(sphere.disp_mass(rho=rho), mass, rtol=5e-2)
    m = sphere.compute_rigid_body_inertia(rho=rho)
    inertia_moment = 2/5 * mass * 1.0**2
    np.testing.assert_allclose(
        np.diag(m.values),
        [mass, mass, mass, inertia_moment, inertia_moment, inertia_moment],
        rtol=5e-2
    )
    np.testing.assert_allclose(m.values - np.diag(np.diag(m.values)), 0.0, atol=1e-5)

# Missing input

def test_inertia_when_no_dofs():
    sphere = cpt.FloatingBody(
        mesh=floating_sphere(),
        center_of_mass=(0, 0, -0.3),
    )
    m = sphere.compute_rigid_body_inertia()
    assert m.shape == (0, 0)

def test_inertia_no_cog():
    sphere = cpt.FloatingBody(
        mesh=floating_sphere(),
        dofs=cpt.rigid_body_dofs()
    )
    with pytest.raises(ValueError, match=".*no center of mass.*"):
        sphere.compute_rigid_body_inertia()

def test_inertia_invariance_by_translation():
    sphere = rigid_body()
    M1 = sphere.compute_rigid_body_inertia()
    M2 = sphere.translated([1.0, 0.0, 0.0]).compute_rigid_body_inertia()
    assert np.allclose(M1, M2)

def test_inertia_single_rotation_dof():
    # Was actually not possible in version 2.x
    mesh = cpt.mesh_parallelepiped()
    body_1 = cpt.FloatingBody(
        mesh=mesh,
        dofs=cpt.rigid_body_dofs(only=["Pitch"], rotation_center=(0, 0, -2)),
        center_of_mass=(0, 0, 0)
    )
    M_1 = body_1.compute_rigid_body_inertia()
    assert M_1.shape == (1, 1)

    ref_body = cpt.FloatingBody(
        mesh=mesh,
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -2)),
        center_of_mass=(0, 0, 0)
    )
    full_M = ref_body.compute_rigid_body_inertia()
    assert np.isclose(M_1.values[0, 0], full_M.values[3, 3])

    body_2 = cpt.FloatingBody(mesh=mesh, center_of_mass=(0, 0, 0))
    body_2.add_rotation_dof(rotation_center=(0, 0, -2), direction=(0, 1, 0), name="pitch oh mon pitch")
    # Non standard name, but still a rigid body rotation
    M_2 = body_2.compute_rigid_body_inertia()
    assert M_2.shape == (1, 1)
    assert np.allclose(M_2.values, M_1.values)

    # Generalized dof approach despite standard name
    body_3 = cpt.FloatingBody(
        mesh=mesh,
        dofs={"Pitch": body_1.dofs["Pitch"].evaluate_motion(mesh)},
        center_of_mass=(0, 0, 0)
    )
    K_3 = body_3.compute_rigid_body_inertia()
    assert K_3.shape == (1, 1)
    assert np.isnan(K_3.values[0, 0])

# Multibody

def test_joining_inertia_when_joining_bodies():
    a, b, _ = several_bodies()
    a.inertia_matrix = a.add_dofs_labels_to_matrix(np.random.rand(1, 1))
    b.inertia_matrix = b.add_dofs_labels_to_matrix(np.random.rand(6, 6))
    both = a + b
    assert both.inertia_matrix.shape == (7, 7)

def test_inertia_joined_bodies_with_missing_inertia():
    a, b, _ = several_bodies()
    # No inertia matrix
    b.inertia_matrix = b.add_dofs_labels_to_matrix(np.random.rand(6, 6))
    both = a + b
    assert not hasattr(both, "inertia_matrix")

def test_inertia_joined_bodies_associativity():
    a, b, c = several_bodies()
    a.inertia_matrix = a.add_dofs_labels_to_matrix(np.array([[1.0]]))
    b.inertia_matrix = b.add_dofs_labels_to_matrix(np.random.rand(6, 6))
    c.inertia_matrix = c.add_dofs_labels_to_matrix(np.random.rand(6, 6))
    i1 = ((a + b) + c).inertia_matrix
    i2 = (a + (b + c)).inertia_matrix
    i3 = cpt.FloatingBody.join_bodies(a, b, c).inertia_matrix
    assert np.allclose(i1.values, i2.values, i3.values)

def test_inertia_joined_bodies_commutativity():
    a, b, _ = several_bodies()
    a.inertia_matrix = a.add_dofs_labels_to_matrix(np.array([[1.0]]))
    b.inertia_matrix = b.add_dofs_labels_to_matrix(np.random.rand(6, 6))
    i1 = (a + b).inertia_matrix
    i2 = (b + a).inertia_matrix
    assert not np.allclose(i1.values, i2.values)
    # but their are the same when the coordinates are sorted in the same way:
    assert np.allclose(i1.sel(influenced_dof=i2.influenced_dof, radiating_dof=i2.radiating_dof).values, i2.values)

######################################################################

##########################
# Full hydrostatics dict #
##########################

def test_hydrostatics_no_rigid_dof():
    mesh = floating_sphere()
    body = cpt.FloatingBody(mesh=mesh, dofs={"Bulge": mesh.faces_normals}, center_of_mass=(0, 0, 0))
    hs = body.compute_hydrostatics()
    assert 'inertia_matrix' not in hs

@pytest.mark.slow
def test_all_hydrostatics():
    density = 1000
    gravity = 9.80665

    sphere = mesh_sphere(
        radius=10.0,
        center=(0,0,-1),
        resolution=(100, 50)
    )
    horizontal_cylinder = mesh_horizontal_cylinder(
        length=10.0, radius=5.0,
        center=(0,10,-1),
        resolution=(100, 10, 100)
    )
    vertical_cylinder = mesh_vertical_cylinder(
        length=10.0, radius=5.0,
        center=(10,0,0),
        resolution=(200, 10, 200)
    )
    full_mesh = sphere + horizontal_cylinder + vertical_cylinder
    body = cpt.FloatingBody(
            mesh=full_mesh,
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            center_of_mass=(1.300, 1.300, -0.870),
            )

    capy_hsdb = body.compute_hydrostatics(rho=density, g=gravity)

    stiff_compare_dofs = ["Heave", "Roll", "Pitch"]
    capy_hsdb["stiffness_matrix"] = capy_hsdb["hydrostatic_stiffness"].sel(
        influenced_dof=stiff_compare_dofs, radiating_dof=stiff_compare_dofs
        ).values

    mass_compare_dofs = ["Roll", "Pitch", "Yaw"]
    capy_hsdb["inertia_matrix"] = capy_hsdb["inertia_matrix"].sel(
        influenced_dof=mass_compare_dofs, radiating_dof=mass_compare_dofs
        ).values


    # =============================================================================
    # Meshmagick
    # =============================================================================
    case_dir = Path(__file__).parent / "Hydrostatics_cases"
    # case_dir.mkdir(parents=True, exist_ok=True)
    # import meshmagick.mesh as mmm
    # import meshmagick.hydrostatics as mmhs
    # body_mesh = mmm.Mesh(full_mesh.vertices, full_mesh.faces, name=body.name)
    # mm_hsdb = mmhs.compute_hydrostatics(body_mesh, body.center_of_mass, density, gravity)
    # rho_body = 1000*body.mesh.disp_volume/body.mesh.volume
    # mm_hsdb["inertia_matrix"] = body_mesh.eval_plain_mesh_inertias(
    #                                          rho_medium=rho_body).inertia_matrix
    # mm_hsdb["mesh"] = ""
    # with open(f'{case_dir}/sphere__hor_cyl__ver_cyl.pkl.json', 'w') as convert_file:
    #     mm_hsdb_json = {key:(value.tolist() if type(value)==np.ndarray else value)
    #                         for key, value in mm_hsdb.items() }
    #     convert_file.write(json.dumps(mm_hsdb_json))

    with open(f'{case_dir}/sphere__hor_cyl__ver_cyl.pkl.json', 'r',
              encoding="UTF-8") as f:
        mm_hsdb = json.load(f)

    # =============================================================================
    # Testing
    # =============================================================================

    for var in capy_hsdb.keys():
        if var in mm_hsdb.keys():
            # if not np.isclose(capy_hsdb[var], mm_hsdb[var],
            #                   rtol=1e-2, atol=1e-3).all():
            #     print(f"{var}:")
            #     print(f"    Capytaine  - {capy_hsdb[var]}")
            #     print(f"    Meshmagick - {mm_hsdb[var]}")
            assert np.isclose(capy_hsdb[var], mm_hsdb[var],
                              rtol=1e-2, atol=1e-3).all()


######################################################################

###########################
#  Non-neutrally buoyant  #
###########################

def test_mass_joined_bodies():
    a = cpt.FloatingBody(mass=100, name="body_1")
    b = cpt.FloatingBody(mass=300, name="body_2")
    assert (a + b).as_FloatingBody().mass == 400

def test_mass_joined_bodies_with_missing_mass():
    a = cpt.FloatingBody(name="body_1")
    b = cpt.FloatingBody(mass=300, name="body_2")
    assert (a + b).as_FloatingBody().mass is None

def test_center_of_mass_joined_bodies():
    a = cpt.FloatingBody(mass=100, center_of_mass=(0, 0, 0), name="body_1")
    b = cpt.FloatingBody(mass=300, center_of_mass=(1, 0, 0), name="body_2")
    assert np.allclose((a + b).as_FloatingBody().center_of_mass, (0.75, 0, 0))

def test_center_of_mass_joined_bodies_with_missing_mass():
    a = cpt.FloatingBody(name="body_1")
    b = cpt.FloatingBody(mass=300, center_of_mass=(1, 0, 0), name="body_2")
    assert (a + b).as_FloatingBody().center_of_mass is None

def test_not_single_rigid_and_non_neutrally_buoyant_body():
    m = mesh_sphere()
    a = cpt.FloatingBody(
        mesh=m, mass=100, center_of_mass=(0, 0, 0),
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
        name="body_1",
    )
    b = cpt.FloatingBody(
        mesh=m.translated_x(2.0), mass=300, center_of_mass=(2, 0, 0),
        dofs=cpt.rigid_body_dofs(rotation_center=(2, 0, 0)),
        name="body_2",
    )
    (a + b).compute_hydrostatic_stiffness()

@lru_cache
def non_neutrally_buoyant_body(length=2.0):
    mesh = mesh_vertical_cylinder(
        radius=1.0,
        length=length,
        center=(0.0, 0.0, 0.0),
        resolution=(20, 40, 40)
    )
    return cpt.FloatingBody(
            mesh=mesh,
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            mass = 500 * mesh.disp_volume,
            center_of_mass = (0.0, 0.0, -0.25),
            )

def test_non_neutrally_buoyant_stiffness():
    body = non_neutrally_buoyant_body(length=2.0)
    K = body.compute_hydrostatic_stiffness()
    dofs = ["Heave", "Roll", "Pitch", "Yaw"]
    K = K.sel(influenced_dof=dofs, radiating_dof=dofs).values
    # K is now a 4x4 np.array in correct order

    rho_g = 1000*9.81
    K_ref = rho_g * np.array([
        [3.137, -0.382e-3, -0.613e-4, 0.0       ],
        [0.0,   -0.392,    -0.276e-4, -0.448e-4 ],
        [0.0,   0.0,       -0.392,    0.313e-3  ],
        [0.0,   0.0,       0.0,       0.0       ],
        ])  # Computed with WAMIT
    print(K)
    print(K_ref)
    assert np.allclose(K, K_ref, atol=rho_g*1e-2)

@pytest.mark.parametrize("data", [
        dict(body_density=1000, z_cog=-0.00, K55=-7.70e3),
        dict(body_density=1000, z_cog=-0.25, K55=5.0),
        dict(body_density=1000, z_cog=-0.50, K55=7.67e3),
        dict(body_density=500,  z_cog=-0.00, K55=-7.67e3),
        dict(body_density=500,  z_cog=-0.25, K55=-3.83e3),
        dict(body_density=500,  z_cog=-0.50, K55=5.0),
        ])
def test_non_neutrally_buoyant_K55(data):
    mesh = mesh_vertical_cylinder(
        radius=1.0,
        length=2.0,
        center=(0.0, 0.0, 0.0),
        resolution=(20, 40, 40)
    )
    body = cpt.FloatingBody(
            mesh=mesh,
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            mass = data["body_density"] * mesh.disp_volume,
            center_of_mass = (0.0, 0.0, data["z_cog"]),
            )
    rho_g = 1000*9.81
    K55 = body.compute_hydrostatic_stiffness().sel(
        influenced_dof="Pitch",
        radiating_dof="Pitch"
    ).values
    assert np.isclose(K55, data["K55"], atol=rho_g*1e-2)

def test_non_neutrally_buoyant_stiffness_invariance_by_translation():
    body = non_neutrally_buoyant_body()
    K1 = body.compute_hydrostatic_stiffness()
    K2 = body.translated([1.0, 0.0, 0.0]).compute_hydrostatic_stiffness()
    assert np.allclose(K1, K2)

def test_non_neutrally_buoyant_inertia():
    mesh = mesh_vertical_cylinder(
        radius=1.0,
        length=1.0,
        center=(0.0, 0.0, -0.5),
        resolution=(20, 40, 20)
    )
    body = cpt.FloatingBody(
            mesh=mesh,
            dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
            mass = 500 * mesh.disp_volume,
            center_of_mass = (0.0, 0.0, -0.25),
            )
    M = body.compute_rigid_body_inertia().values
    assert np.allclose(np.diag(M), np.array([1570, 1570, 1570, 936, 936, 801]), rtol=5e-2)

def test_non_neutrally_buoyant_inertia_invariance_by_translation():
    body = non_neutrally_buoyant_body()
    M1 = body.compute_rigid_body_inertia()
    M2 = body.translated([1.0, -1.0, 0.0]).compute_rigid_body_inertia()
    assert np.allclose(M1, M2)
