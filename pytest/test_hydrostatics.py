"""
Tests for the functions that computes hydrostatic from the mesh vertices
and faces
"""

import logging
import pytest
import json
from pathlib import Path

import pandas as pd
import numpy as np
import capytaine as cpt

# TODO: Can I use pytest fixtures to avoid regenerating so many spheres?
# Does pytest fixture do a copy of the object, since I modifying the sphere in-place?

def test_mesh_properties():
    sphere = cpt.Sphere(radius=1.0, center=(0, 0, -2), nphi=50, ntheta=50)
    assert np.allclose(sphere.center_of_buoyancy, sphere.mesh.center_of_buoyancy)

def test_disp_mass_of_sphere():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,-2), nphi=50, ntheta=50)
    analytical_volume = 4/3*np.pi*1.0**3
    assert np.isclose(sphere.disp_mass(rho=1000), 1000*analytical_volume, rtol=2e-2)

def test_waterplane_area_of_submerged_sphere():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,-2), nphi=20, ntheta=20)
    assert np.isclose(sphere.waterplane_area, 0.0)

def test_waterplane_center_of_submerged_sphere():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,-2), nphi=20, ntheta=20)
    assert sphere.waterplane_center is None

def test_waterplane_center_of_sphere_at_surface():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    assert np.allclose(sphere.waterplane_center, [0.0, 0.0])

def test_infer_rotation_center():
    sphere = cpt.Sphere().keep_immersed_part()
    sphere.rotation_center = (7, 8, 9)
    sphere.add_all_rigid_body_dofs()
    assert np.allclose(sphere._infer_rotation_center(), (7, 8, 9))
    del sphere.rotation_center
    assert np.allclose(sphere._infer_rotation_center(), (7, 8, 9))

def test_stiffness_when_no_dofs():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.center_of_mass = np.array([0, 0, -0.3])
    with pytest.raises(AttributeError, match=".* no dof .*"):
        sphere.compute_hydrostatic_stiffness()

def test_stiffness_no_cog():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.add_all_rigid_body_dofs()
    with pytest.raises(ValueError, match=".*no center of mass.*"):
        sphere.compute_hydrostatic_stiffness()

def test_stiffness_dof_ordering():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20)
    sphere.keep_immersed_part()
    sphere.add_all_rigid_body_dofs()
    sphere.center_of_mass = np.array([0.0, 0.0, -0.2])
    K = sphere.compute_hydrostatic_stiffness()
    assert np.all(K.coords["radiating_dof"].values == np.array(['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']))

def test_stifness_rigid_body_invariance_by_translation():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20)
    sphere.keep_immersed_part()
    sphere.center_of_mass = np.array([0.0, 0.0, -0.2])
    sphere.add_all_rigid_body_dofs()
    K1 = sphere.compute_hydrostatic_stiffness()
    K2 = sphere.translated([1.0, 0.0, 0.0]).compute_hydrostatic_stiffness()
    assert np.allclose(K1, K2)

def test_stifness_generalized_dof_invariance_by_translation():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20)
    sphere.keep_immersed_part()
    # Not really a generalized dof, but since the dofs have non-standard names
    # the code will not recognize the rigid body dof.
    sphere.add_translation_dof(direction=(1, 0, 0), name="cavalement")
    sphere.add_rotation_dof(axis=cpt.Axis(vector=(0, 1, 0), point=(0, 0, 0)), name="tangage")
    K1 = sphere.compute_hydrostatic_stiffness()
    K2 = sphere.translated([1.0, 0.0, 0.0]).compute_hydrostatic_stiffness()
    assert np.allclose(K1, K2)

def test_stiffness_with_divergence():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.center_of_mass = np.array([0, 0, -0.3])
    sphere.dofs["elongate_in_z"] = np.array([(0, 0, z) for (x, y, z) in sphere.mesh.faces_centers])
    hs_1 = sphere.compute_hydrostatic_stiffness()
    hs_2 = sphere.compute_hydrostatic_stiffness(divergence={"elongate_in_z": np.ones(sphere.mesh.nb_faces)})
    assert hs_1.values[0, 0] != hs_2.values[0, 0]
    analytical_hs = - 1000.0 * 9.81 * (4 * sphere.volume * sphere.center_of_buoyancy[2])
    assert np.isclose(hs_2.values[0, 0], analytical_hs)

def test_stiffness_with_malformed_divergence(caplog):
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.center_of_mass = np.array([0, 0, -0.3])
    sphere.dofs["elongate_in_z"] = np.array([(0, 0, z) for (x, y, z) in sphere.mesh.faces_centers])
    hs_1 = sphere.compute_hydrostatic_stiffness()
    with caplog.at_level(logging.WARNING):
        hs_2 = sphere.compute_hydrostatic_stiffness(divergence={"foobar": np.ones(sphere.mesh.nb_faces)})
    assert hs_1.values[0, 0] == hs_2.values[0, 0]
    assert "without the divergence" in caplog.text

def test_stiffness_joined_bodies():
    a = cpt.Sphere(name="foo")
    a.add_translation_dof(name="Heave")
    a.hydrostatic_stiffness = a.add_dofs_labels_to_matrix(np.array([[1.0]]))
    b = cpt.RectangularParallelepiped(name="bar")
    b.add_all_rigid_body_dofs()
    b.hydrostatic_stiffness = b.add_dofs_labels_to_matrix(np.random.rand(6, 6))
    both = a + b
    assert both.hydrostatic_stiffness.shape == (7, 7)

def test_mass_of_sphere_for_non_default_density():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,-2), nphi=50, ntheta=50)
    sphere.add_all_rigid_body_dofs()
    sphere.center_of_mass = np.array([0, 0, -2])
    m = sphere.compute_rigid_body_inertia(rho=500)
    analytical_volume = 4/3*np.pi*1.0**3
    assert np.isclose(m.values[0, 0], 500*analytical_volume, rtol=2e-2)

def test_inertia_rigid_body_dofs():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.center_of_mass = np.array([0, 0, -0.3])
    sphere.add_all_rigid_body_dofs()
    assert np.all(sphere.compute_rigid_body_inertia(output_type="rigid_dofs")
            == sphere.compute_rigid_body_inertia(output_type="all_dofs"))
    assert np.all(sphere.compute_rigid_body_inertia(output_type="rigid_dofs")
            == sphere.compute_rigid_body_inertia(output_type="body_dofs"))

def test_inertia_invariance_by_translation():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20)
    sphere.keep_immersed_part()
    sphere.center_of_mass = np.array([0.0, 0.0, -0.2])
    sphere.add_all_rigid_body_dofs()
    M1 = sphere.compute_rigid_body_inertia()
    M2 = sphere.translated([1.0, 0.0, 0.0]).compute_rigid_body_inertia()
    assert np.allclose(M1, M2)

def test_inertia_wrong_output_type():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.center_of_mass = np.array([0, 0, -0.3])
    sphere.add_all_rigid_body_dofs()
    with pytest.raises(ValueError):
        sphere.compute_rigid_body_inertia(output_type="foo")

def test_inertia_when_no_dofs():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.rotation_center = np.array([0, 0, 0])
    sphere.center_of_mass = np.array([0, 0, -0.3])
    m =sphere.compute_rigid_body_inertia()
    assert m.shape == (0, 0)

def test_inertia_no_cog():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20).keep_immersed_part()
    sphere.add_all_rigid_body_dofs()
    with pytest.raises(ValueError, match=".*no center of mass.*"):
        sphere.compute_rigid_body_inertia()

def test_inertia_joined_bodies():
    a = cpt.Sphere(name="foo")
    a.add_translation_dof(name="Heave")
    a.inertia_matrix = a.add_dofs_labels_to_matrix(np.array([[1.0]]))
    b = cpt.RectangularParallelepiped(name="bar")
    b.add_all_rigid_body_dofs()
    b.inertia_matrix = b.add_dofs_labels_to_matrix(np.random.rand(6, 6))
    both = a + b
    assert both.inertia_matrix.shape == (7, 7)

def test_inertia_joined_bodies_with_missing_inertia():
    a = cpt.Sphere(name="foo")
    a.add_translation_dof(name="Heave")
    # No inertia matrix
    b = cpt.RectangularParallelepiped(name="bar")
    b.add_all_rigid_body_dofs()
    b.inertia_matrix = b.add_dofs_labels_to_matrix(np.random.rand(6, 6))
    both = a + b
    assert not hasattr(both, "inertia_matrix")

def test_inertia_joined_bodies_associativity():
    a = cpt.Sphere(name="foo")
    a.add_translation_dof(name="Heave")
    a.inertia_matrix = a.add_dofs_labels_to_matrix(np.array([[1.0]]))
    b = cpt.RectangularParallelepiped(name="bar")
    b.add_all_rigid_body_dofs()
    b.inertia_matrix = b.add_dofs_labels_to_matrix(np.random.rand(6, 6))
    c = cpt.VerticalCylinder(name="baz")
    c.add_rotation_dof(name="Pitch")
    c.inertia_matrix = c.add_dofs_labels_to_matrix(np.array([[2.0]]))
    i1 = ((a + b) + c).inertia_matrix
    i2 = (a + (b + c)).inertia_matrix
    i3 = cpt.FloatingBody.join_bodies(a, b, c).inertia_matrix
    assert np.allclose(i1.values, i2.values, i3.values)

def test_inertia_joined_bodies_commutativity():
    a = cpt.Sphere(name="foo")
    a.add_translation_dof(name="Heave")
    a.inertia_matrix = a.add_dofs_labels_to_matrix(np.array([[1.0]]))
    b = cpt.RectangularParallelepiped(name="bar")
    b.add_all_rigid_body_dofs()
    b.inertia_matrix = b.add_dofs_labels_to_matrix(np.random.rand(6, 6))
    i1 = (a + b).inertia_matrix
    i2 = (b + a).inertia_matrix
    assert not np.allclose(i1.values, i2.values)
    # but their are the same when the coordinates are sorted in the same way:
    assert np.allclose(i1.sel(influenced_dof=i2.influenced_dof, radiating_dof=i2.radiating_dof).values, i2.values)


def test_hydrostatics_of_submerged_sphere():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,-2), nphi=20, ntheta=20)
    sphere.add_all_rigid_body_dofs()
    sphere.center_of_mass = np.array([0, 0, -2])
    sphere.compute_hydrostatics()
    assert sphere.inertia_matrix.shape == (6, 6)

def test_hydrostatics_no_rigid_dof():
    mesh = cpt.mesh_sphere()
    body = cpt.FloatingBody(mesh=mesh, dofs={"Bulge": mesh.faces_normals}, center_of_mass=(0, 0, 0))
    hs = body.compute_hydrostatics()
    assert 'inertia_matrix' not in hs

def test_all_hydrostatics():
    density = 1000
    gravity = 9.80665

    sphere = cpt.Sphere(
        radius=10.0,
        center=(0,0,-1),
        nphi=100, ntheta=50,
    )
    horizontal_cylinder = cpt.HorizontalCylinder(
        length=10.0, radius=5.0,
        center=(0,10,-1),
        nr=100, nx=100, ntheta=10,
    )
    vertical_cylinder = cpt.VerticalCylinder(
        length=10.0, radius=5.0,
        center=(10,0,0),
        nr=200, nx=200, ntheta=10,
    )
    body = sphere + horizontal_cylinder + vertical_cylinder
    body.rotation_center = np.array([0, 0, 0])
    body.center_of_mass = body.center_of_buoyancy
    body.add_all_rigid_body_dofs()

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
    # body_mesh = mmm.Mesh(body.mesh.vertices, body.mesh.faces, name=body.mesh.name)
    # mm_hsdb = mmhs.compute_hydrostatics(body_mesh, body.center_of_mass, density, gravity)
    # mm_hsdb["inertia_matrix"] = body_mesh.eval_plain_mesh_inertias(
    # rho_medium=density).inertia_matrix
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


def test_vertical_elastic_dof():

    bodies = [
        cpt.VerticalCylinder(
            length=10.0, radius=5.0,
            center=(0.0,0.0,0.0),
            nr=100, nx=100, ntheta=100,
        ),

        cpt.Sphere(
            radius=100,
            center=(0,0,0),
            nphi=50, ntheta=50,
        ),

        cpt.HorizontalCylinder(
            length=5.0, radius=1.0,
            center=(0,10,0),
            reflection_symmetry=False,
            nr=20, nx=20, ntheta=10,
        )
    ]

    for body in bodies:
        body.center_of_mass = body.center_of_buoyancy
        body.keep_immersed_part()
        faces_centers = body.mesh.faces_centers
        body.dofs["elongate"] = np.zeros_like(faces_centers)
        body.dofs["elongate"][:,2] = faces_centers[:,2]

        divergence = np.ones(body.mesh.faces_centers.shape[0])

        density = 1000
        gravity = 9.80665
        capy_hs = body.each_hydrostatic_stiffness("elongate", "elongate",
                    influenced_dof_div=divergence, rho=density, g=gravity).values[0][0]

        analytical_hs = - (density * gravity * 4 * body.volume
                        * body.center_of_buoyancy[2])

        assert np.isclose(capy_hs, analytical_hs)


###########################
#  Non-neutrally buoyant  #
###########################

def test_mass_joined_bodies():
    a = cpt.FloatingBody(mass=100)
    b = cpt.FloatingBody(mass=300)
    assert (a + b).mass == 400

def test_mass_joined_bodies_with_missing_mass():
    a = cpt.FloatingBody()
    b = cpt.FloatingBody(mass=300)
    assert (a + b).mass is None

def test_center_of_mass_joined_bodies():
    a = cpt.FloatingBody(mass=100, center_of_mass=(0, 0, 0))
    b = cpt.FloatingBody(mass=300, center_of_mass=(1, 0, 0))
    assert np.allclose((a + b).center_of_mass, (0.75, 0, 0))

def test_center_of_mass_joined_bodies_with_missing_mass():
    a = cpt.FloatingBody()
    b = cpt.FloatingBody(mass=300, center_of_mass=(1, 0, 0))
    assert (a + b).center_of_mass is None

def test_not_single_rigid_and_non_neutrally_buoyant_body():
    m = cpt.mesh_sphere()
    a = cpt.FloatingBody(
        mesh=m, mass=100, center_of_mass=(0, 0, 0),
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
    )
    b = cpt.FloatingBody(
        mesh=m.translated_x(2.0), mass=300, center_of_mass=(2, 0, 0),
        dofs=cpt.rigid_body_dofs(rotation_center=(2, 0, 0)),
    )
    with pytest.raises(NotImplementedError):
        (a + b).compute_hydrostatic_stiffness()

def test_non_neutrally_buoyant_stiffness():
    body = cpt.VerticalCylinder(radius=1.0, length=2.0, center=(0.0, 0.0, 0.0), nx=40, ntheta=40, nr=20)
    body.rotation_center = (0, 0, 0)
    body.add_all_rigid_body_dofs()
    body.keep_immersed_part()
    body.mass = 500 * body.volume
    body.center_of_mass = (0.0, 0.0, -0.25)

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


def test_non_neutrally_buoyant_K55():
    body = cpt.VerticalCylinder(radius=1.0, length=2.0, center=(0.0, 0.0, 0.0), nx=40, ntheta=40, nr=20)
    body.rotation_center = (0, 0, 0)
    body.add_all_rigid_body_dofs()
    body.keep_immersed_part()

    ref_data = pd.DataFrame([
        dict(body_density=1000, z_cog=-0.00, K55=-7.70e3),
        dict(body_density=1000, z_cog=-0.25, K55=5.0),
        dict(body_density=1000, z_cog=-0.50, K55=7.67e3),
        dict(body_density=500,  z_cog=-0.00, K55=-7.67e3),
        dict(body_density=500,  z_cog=-0.25, K55=-3.83e3),
        dict(body_density=500,  z_cog=-0.50, K55=5.0),
        ])
    ref_data['mass'] = body.volume * ref_data['body_density']

    rho_g = 1000*9.81
    for (i, case) in ref_data.iterrows():
        body.mass = case.mass
        body.center_of_mass = (0.0, 0.0, case.z_cog)
        K55 = body.compute_hydrostatic_stiffness().sel(influenced_dof="Pitch", radiating_dof="Pitch").values
        assert np.isclose(K55, case.K55, atol=rho_g*1e-2)


def test_non_neutrally_buoyant_stiffness_invariance_by_translation():
    sphere = cpt.Sphere(radius=1.0, center=(0,0,0), nphi=20, ntheta=20)
    sphere.keep_immersed_part()
    sphere.mass = 300 * sphere.volume
    sphere.center_of_mass = np.array([0.0, 0.0, -0.2])
    sphere.add_all_rigid_body_dofs()
    K1 = sphere.compute_hydrostatic_stiffness()
    K2 = sphere.translated([1.0, 0.0, 0.0]).compute_hydrostatic_stiffness()
    assert np.allclose(K1, K2)

def test_non_neutrally_buoyant_inertia():
    body = cpt.VerticalCylinder(radius=1.0, length=1.0, center=(0.0, 0.0, -0.5), nx=20, ntheta=40, nr=20)
    body.rotation_center = (0, 0, 0)
    body.add_all_rigid_body_dofs()
    body.keep_immersed_part()
    body.center_of_mass = (0.0, 0.0, -0.25)
    body.mass = body.volume * 500
    M = body.compute_rigid_body_inertia().values
    assert np.allclose(np.diag(M), np.array([1570, 1570, 1570, 936, 936, 801]), rtol=5e-2)

def test_non_neutrally_buoyant_inertia_invariance_by_translation():
    sphere = cpt.Sphere(radius=1.0, center=(0, 0, 0), nphi=20, ntheta=20)
    sphere.rotation_center = np.array([0.0, 0.0, 0.0])
    sphere.mass = 300 * sphere.volume
    sphere.center_of_mass = np.array([0.0, 0.0, -0.2])
    sphere.add_all_rigid_body_dofs()
    M1 = sphere.compute_rigid_body_inertia()
    M2 = sphere.translated([1.0, -1.0, 10.0]).compute_rigid_body_inertia()
    print(M1.values)
    print(M2.values)
    assert np.allclose(M1, M2)
