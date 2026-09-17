# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import pytest

import numpy as np
import capytaine as cpt

@pytest.mark.parametrize("radius", [0.8, 1., 2.8])
def test_integration_on_circle(radius):
    mesh_sphere = cpt.mesh_sphere(radius=radius, resolution=(10, 50), axial_symmetry=True).immersed_part()
    integral_sphere = mesh_sphere.waterline_integral(np.ones(mesh_sphere.nb_edges_waterline))
    assert np.isclose(integral_sphere, 2*np.pi*radius, rtol=1e-3)


@pytest.mark.parametrize("l, L, h", [(2, 3, 4), (1, 1, 1), (4.8, 0.1, 2.1)])
def test_integration_on_rectangle(l, L, h):
    mesh_rectangle = cpt.mesh_parallelepiped(size=(l,L,h), reflection_symmetry=True).immersed_part()
    integral_rectangle = mesh_rectangle.waterline_integral(np.ones(mesh_rectangle.nb_edges_waterline))
    assert np.isclose(integral_rectangle, 2*(l+L))
