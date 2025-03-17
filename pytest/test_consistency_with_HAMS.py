import numpy as np
import capytaine as cpt


vertices = np.array([
    [-0.6, -0.2, -0.5],
    [-0.6, -0.4, -0.5],
    [-0.8, -0.4, -0.5],
    [-0.8, -0.2, -0.5],
    [-0.6,  0.0, -0.5],
    [-0.8,  0.0, -0.5],
    ])
faces = np.array([[0, 1, 2, 3], [4, 0, 3, 5]])
mesh = cpt.Mesh(vertices, faces)


def test_S_matrix_pure_Rankine():
    HAMS_S_Rankine = np.array([
        [0.70509886961563439, 0.20760994718095135],
        [0.20760994718095122, 0.70509886961563462],
        ])
    S_rankine, _ = cpt.Delhommeau().evaluate(mesh, mesh, free_surface=np.inf, water_depth=np.inf, wavenumber=0.0, adjoint_double_layer=False)
    np.testing.assert_allclose(-4*np.pi*S_rankine, HAMS_S_Rankine, rtol=1e-3)


def test_dispersion_roots():
    k = 1.0
    h = 10.0
    HAMS_roots = np.array([0.99996435, 0.17434080, 0.51912325, 0.85620867, 1.1866127, 1.5123763, 1.8350620, 2.1556817, 2.4748751, 2.7930533, 3.1104880, 3.4273639, 3.7438106, 4.0599209, 4.3757631, 4.6913883, 5.0068356, 5.3221352, 5.6373110, 5.9523820,])
    omega2_over_g = k * np.tanh(k*h)
    roots = cpt.FinGreen3D().fortran_core.green_wave.compute_dispersion_roots(20, omega2_over_g, 10.0)
    np.testing.assert_allclose(HAMS_roots, roots, rtol=1e-3)

def test_S_matrix_finite_depth_k1():
    k = 1.0
    h = 10.0
    HAMS_S_matrix = np.array([
        [(0.73247863590752549 + 1j * 2.3114545333287739) , (-0.43898677720635493 + 1j * 2.2883993524884838) ,],
        [(-0.43898677720635493 + 1j * 2.2883993524884838),   (0.73247863590752549 + 1j * 2.3114545333287739), ],
        ]) * mesh.faces_areas[0]
    S, _ = cpt.FinGreen3D().evaluate(mesh, mesh, free_surface=0.0, water_depth=h, wavenumber=k, adjoint_double_layer=False)
    S_rankine, _ = cpt.Delhommeau().evaluate(mesh, mesh, free_surface=np.inf, water_depth=np.inf, wavenumber=0.0, adjoint_double_layer=False)
    np.testing.assert_allclose(-4*np.pi*(S - S_rankine), HAMS_S_matrix, rtol=1e-3)


def test_S_matrix_finite_depth_k2():
    k = 2.0
    h = 10.0
    HAMS_S_matrix = np.array([
                [(-0.36177831640996905 + 1j * 1.7007340213300215),   (-1.7129415732152193 + 1j * 1.6333866990767409)],
                [(-1.7129415732152193 + 1j * 1.6333866990767409) , (-0.36177831640996905 + 1j * 1.7007340213300215)],
                  ]) * mesh.faces_areas[0]
    S, _ = cpt.FinGreen3D().evaluate(mesh, mesh, free_surface=0.0, water_depth=h, wavenumber=k, adjoint_double_layer=False)
    S_rankine, _ = cpt.Delhommeau().evaluate(mesh, mesh, free_surface=np.inf, water_depth=np.inf, wavenumber=0.0, adjoint_double_layer=False)
    np.testing.assert_allclose(-4*np.pi*(S - S_rankine), HAMS_S_matrix, rtol=1e-3)
