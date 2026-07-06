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
from capytaine.bem.airy_waves import froude_krylov_force


def test_Froude_Krylov():
    sphere = cpt.FloatingBody(cpt.mesh_sphere(radius=1.0, resolution=(6, 12)).immersed_part())
    sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")

    problem = cpt.DiffractionProblem(body=sphere, omega=1.0, water_depth=np.inf)
    assert np.isclose(froude_krylov_force(problem)['Heave'], 27596, rtol=1e-3)

    problem = cpt.DiffractionProblem(body=sphere, omega=2.0, water_depth=np.inf)
    assert np.isclose(froude_krylov_force(problem)['Heave'], 22491, rtol=1e-3)

    problem = cpt.DiffractionProblem(body=sphere, omega=1.0, water_depth=10.0)
    assert np.isclose(froude_krylov_force(problem)['Heave'], 27610, rtol=1e-3)
