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

import numpy as np

import capytaine as cpt
from capytaine.bem.problems_checks import _check_wavelength_and_water_depth
from capytaine.meshes.predefined import mesh_parallelepiped


def test_warning_mesh_resolution(caplog):
    mesh = cpt.mesh_sphere(radius=1.0, resolution=(4, 4)).immersed_part()
    sphere = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    solver = cpt.BEMSolver()
    pb = cpt.RadiationProblem(body=sphere, wavelength=0.1*sphere.minimal_computable_wavelength)
    with caplog.at_level("WARNING"):
        solver.solve(pb)
    assert "resolution " in caplog.text


def test_warning_for_deep_water_single_problem(caplog):
    """Test warning for a single problem with very deep water."""
    mesh = mesh_parallelepiped(size=(1, 1, 1), center=(0, 0, -0.5))
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    pb = cpt.DiffractionProblem(body=body, omega=1.0, water_depth=500.0)
    solver = cpt.BEMSolver()
    with caplog.at_level("WARNING"):
        solver.solve(pb)
    assert "Water depth" in caplog.text
    assert "water_depth=500.0" in caplog.text


def test_warning_for_deep_water_multiple_problems(caplog):
    """Test warning for multiple problems with deep water."""
    mesh = mesh_parallelepiped(size=(1, 1, 1), center=(0, 0, -0.5))
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    problems = [
        cpt.DiffractionProblem(body=body, omega=1.0, water_depth=500.0),
        cpt.DiffractionProblem(body=body, omega=1.5, water_depth=500.0),
    ]
    solver = cpt.BEMSolver()
    with caplog.at_level("WARNING"):
        solver.solve_all(problems)
    assert "Water depth for 2 problems" in caplog.text
    assert "wavelength ranging from" in caplog.text
