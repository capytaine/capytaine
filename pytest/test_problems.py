#!/usr/bin/env python
# coding: utf-8

import logging

import pytest

import numpy as np

from capytaine.problems import *
from capytaine.results import *
from capytaine.import_export import import_cal_file
from capytaine.reference_bodies import generate_sphere


def test_LinearPotentialFlowProblem():
    pb = LinearPotentialFlowProblem(omega=1.0)
    assert pb.omega == 1.0
    assert pb.wavenumber == 1.0/9.81
    assert pb.wavelength == 9.81*2*np.pi

    res = pb.make_results_container()
    assert isinstance(res, LinearPotentialFlowResult)
    assert res.problem is pb
    assert res.omega == pb.omega
    assert res.dimensionless_omega == pb.dimensionless_omega
    assert res.body is pb.body
    assert res.forces == {}

    assert LinearPotentialFlowProblem(free_surface=np.infty, sea_bottom=-np.infty).depth == np.infty
    assert LinearPotentialFlowProblem(free_surface=0.0, sea_bottom=-np.infty).depth == np.infty
    assert LinearPotentialFlowProblem(free_surface=0.0, sea_bottom=-1.0).depth == 1.0

    with pytest.raises(ValueError):
        LinearPotentialFlowProblem(free_surface=0.0, sea_bottom=1.0)

    with pytest.raises(TypeError):
        LinearPotentialFlowProblem(angle=1.0)

    with pytest.raises(TypeError):
        LinearPotentialFlowProblem(radiating_dof="Heave")

    # LinearPotentialFlowProblem(boundary_condition=[1, 1, 1])


def test_diffraction_problem():
    assert DiffractionProblem().body is None

    sphere = generate_sphere(radius=1.0, ntheta=20, nphi=40)
    sphere.dofs["Heave"] = sphere.faces_normals @ (0, 0, 1)

    pb = DiffractionProblem(body=sphere, angle=1.0)
    assert len(pb.boundary_condition) == sphere.nb_faces

    with pytest.raises(TypeError):
        DiffractionProblem(boundary_conditions=[0, 0, 0])

    res = pb.make_results_container()
    assert isinstance(res, DiffractionResult)


def test_radiation_problem(caplog):
    sphere = generate_sphere(radius=1.0, ntheta=20, nphi=40)

    # with pytest.raises(ValueError):
    #     RadiationProblem(body=sphere)

    sphere.dofs["Heave"] = sphere.faces_normals @ (0, 0, 1)
    pb = RadiationProblem(body=sphere)
    assert len(pb.boundary_condition) == sphere.nb_faces

    sphere.dofs["Surge"] = sphere.faces_normals @ (1, 0, 0)
    pb2 = RadiationProblem(body=sphere, radiating_dof="Heave")
    assert np.all(pb.boundary_condition == pb2.boundary_condition)

    res = pb.make_results_container()
    assert isinstance(res, RadiationResult)
    assert res.forces is None
    assert res.added_masses == {}
    assert res.radiation_dampings == {}


# def test_import_cal_file():
#     """Test the importation of legacy Nemoh.cal files."""
#     problems = import_cal_file("examples/data/Nemoh.cal")
#     assert len(problems) == 4
#     for problem in problems:
#         assert problem.rho == 1000.0
#         assert problem.g == 9.81
#         assert problem.depth == np.infty
#         assert problem.body.nb_subbodies == 1
#         assert problem.body.nb_dofs == 6
#         assert problem.body.nb_vertices == 540
#         assert problem.body.nb_faces == 300
#         assert problem.omega == 0.1 or problem.omega == 2.0
#         if isinstance(problem, DiffractionProblem):
#             assert problem.angle == 0.0



