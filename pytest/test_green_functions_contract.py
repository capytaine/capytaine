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
"""Tests of the contract shared by all the Green functions: shapes of the matrices returned by `evaluate`."""

from functools import lru_cache

import pytest

import numpy as np
import capytaine as cpt


# Each Green function with the configurations (free_surface, water_depth) in which it is expected to work.
INF_DEPTH = {"free_surface": 0.0, "water_depth": np.inf}
FINITE_DEPTH = {"free_surface": 0.0, "water_depth": 5.0}
NO_FREE_SURFACE = {"free_surface": np.inf, "water_depth": np.inf}

GREEN_FUNCTIONS = [
    (cpt.Delhommeau, INF_DEPTH),
    (cpt.Delhommeau, FINITE_DEPTH),
    (cpt.Delhommeau, NO_FREE_SURFACE),
    (cpt.LiangWuNoblesseGF, INF_DEPTH),
    (cpt.FinGreen3D, FINITE_DEPTH),
    (cpt.HAMS_GF, INF_DEPTH),
    (cpt.HAMS_GF, FINITE_DEPTH),
]


def _id(param):
    gf_class, config = param
    return f"{gf_class.__name__}-depth={config['water_depth']}-fs={config['free_surface']}"


@lru_cache
def sphere_mesh():
    return cpt.mesh_sphere(radius=1, center=(0, 0, -2), resolution=(6, 6)).immersed_part()


@pytest.mark.parametrize("gf_class, config", GREEN_FUNCTIONS, ids=[_id(p) for p in GREEN_FUNCTIONS])
@pytest.mark.parametrize("adjoint_double_layer", [True, False])
@pytest.mark.parametrize("early_dot_product", [True, False])
def test_shapes_of_S_and_K(gf_class, config, adjoint_double_layer, early_dot_product):
    mesh = sphere_mesh()
    n = mesh.nb_faces
    S, K = gf_class().evaluate(
        mesh, mesh, wavenumber=1.0,
        adjoint_double_layer=adjoint_double_layer, early_dot_product=early_dot_product,
        **config,
    )
    assert S.shape == (n, n)
    assert K.shape == ((n, n) if early_dot_product else (3, n, n))


@pytest.mark.parametrize("gf_class, config", GREEN_FUNCTIONS, ids=[_id(p) for p in GREEN_FUNCTIONS])
@pytest.mark.parametrize("early_dot_product", [True, False])
def test_evaluate_at_array_of_points(gf_class, config, early_dot_product):
    mesh = sphere_mesh()
    points = np.array([[2.0, 0.0, -1.0], [0.0, 3.0, -2.0], [1.0, 1.0, -4.0]])
    # The diagonal term requires the normals of the receiving mesh, so it cannot be used with points.
    S, K = gf_class().evaluate(
        points, mesh, wavenumber=1.0, early_dot_product=early_dot_product,
        adjoint_double_layer=True, diagonal_term_in_double_layer=False,
        **config,
    )
    assert S.shape == (3, mesh.nb_faces)
    assert K.shape == ((3, mesh.nb_faces) if early_dot_product else (3, 3, mesh.nb_faces))


@pytest.mark.parametrize("gf_class", [cpt.Delhommeau, cpt.LiangWuNoblesseGF, cpt.FinGreen3D, cpt.HAMS_GF])
def test_identical_green_functions_have_same_hash(gf_class):
    assert hash(gf_class()) == hash(gf_class())


def test_different_green_functions_have_different_hash():
    assert hash(cpt.FinGreen3D(nb_dispersion_roots=100)) != hash(cpt.FinGreen3D(nb_dispersion_roots=200))
    assert hash(cpt.Delhommeau(tabulation_nr=100)) != hash(cpt.Delhommeau())


@pytest.mark.parametrize("gf, expected_str, expected_repr", [
    (cpt.LiangWuNoblesseGF(), "LiangWuNoblesseGF()", "LiangWuNoblesseGF()"),
    (cpt.HAMS_GF(), "HAMS_GF()", "HAMS_GF()"),
    (cpt.FinGreen3D(), "FinGreen3D()", "FinGreen3D(nb_dispersion_roots=200)"),
    (cpt.FinGreen3D(nb_dispersion_roots=100), "FinGreen3D(nb_dispersion_roots=100)", "FinGreen3D(nb_dispersion_roots=100)"),
    (cpt.Delhommeau(), "Delhommeau()", None),
    (cpt.Delhommeau(gf_singularities="high_freq"), "Delhommeau(gf_singularities='high_freq')", None),
])
def test_str_and_repr(gf, expected_str, expected_repr):
    assert str(gf) == expected_str
    if expected_repr is not None:
        assert repr(gf) == expected_repr
    else:  # repr shows all the settings
        assert repr(gf).startswith(f"{gf.__class__.__name__}(tabulation_nr=676")
