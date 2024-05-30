import pytest

import numpy as np
import capytaine as cpt


def test_lid_below_free_surface():
    mesh = cpt.AxialSymmetricMesh.from_profile(lambda z: (z + 1.0)**2, np.linspace(-1.0, 0.0, 10)).merged()
    lid_mesh = mesh.generate_lid(z=-0.5)
    x, y, z = lid_mesh.faces_centers.T
    assert np.all(np.hypot(x, y) <= (z + 1.0)**2)
