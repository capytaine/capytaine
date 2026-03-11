import pytest

import numpy as np
import capytaine as cpt 

@pytest.mark.parametrize("radius", [0.8, 1., 2.8])
def test_integration_on_circle(radius):
    mesh_sphere = cpt.mesh_sphere(radius=radius, resolution=(10, 50), axial_symmetry=True).immersed_part()
    integral_sphere = mesh_sphere.water_line_integral(np.ones(mesh_sphere.nb_edges_waterline))
    assert np.isclose(integral_sphere, 2*np.pi*radius, rtol=1e-3)


@pytest.mark.parametrize("l, L, h", [(2, 3, 4), (1, 1, 1), (4.8, 0.1, 2.1)])
def test_integration_on_rectangle(l, L, h): 
    mesh_rectangle = cpt.mesh_parallelepiped(size=(l,L,h), reflection_symmetry=True).immersed_part()
    integral_rectangle = mesh_rectangle.water_line_integral(np.ones(mesh_rectangle.nb_edges_waterline))
    assert np.isclose(integral_rectangle, 2*(l+L))