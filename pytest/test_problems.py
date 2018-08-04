#!/usr/bin/env python
# coding: utf-8

import os

import pytest

import numpy as np

from meshmagick.mesh import Mesh

from capytaine.tools.import_export import import_cal_file

from capytaine.bodies import FloatingBody
from capytaine.symmetries import ReflectionSymmetry
from capytaine.geometric_bodies.sphere import Sphere

from capytaine.problems import LinearPotentialFlowProblem, DiffractionProblem, RadiationProblem
from capytaine.results import LinearPotentialFlowResult, DiffractionResult, RadiationResult


def test_LinearPotentialFlowProblem():
    # Without a body
    pb = LinearPotentialFlowProblem(omega=1.0)
    assert pb.omega == 1.0
    assert pb.wavenumber == 1.0/9.81
    assert pb.wavelength == 9.81*2*np.pi

    assert LinearPotentialFlowProblem(free_surface=np.infty, sea_bottom=-np.infty).depth == np.infty
    assert LinearPotentialFlowProblem(free_surface=0.0, sea_bottom=-np.infty).depth == np.infty
    assert LinearPotentialFlowProblem(free_surface=0.0, sea_bottom=-1.0).depth == 1.0

    with pytest.raises(NotImplementedError):
        LinearPotentialFlowProblem(free_surface=2.0)

    with pytest.raises(ValueError):
        LinearPotentialFlowProblem(free_surface=0.0, sea_bottom=1.0)

    with pytest.raises(TypeError):
        LinearPotentialFlowProblem(angle=1.0)

    with pytest.raises(TypeError):
        LinearPotentialFlowProblem(radiating_dof="Heave")

    # With a body
    sphere = Sphere(center=(0, 0, -2.0))
    sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    pb = LinearPotentialFlowProblem(body=sphere,
                                    boundary_condition=sphere.mesh.faces_normals @ (1, 1, 1))
    assert list(pb.influenced_dofs.keys()) == ['Heave']

    # Test transformation to result class
    res = pb.make_results_container()
    assert isinstance(res, LinearPotentialFlowResult)
    assert res.problem is pb
    assert res.omega == pb.omega
    assert res.dimensionless_omega == pb.dimensionless_omega
    assert res.body is pb.body


def test_diffraction_problem():
    assert DiffractionProblem().body is None

    sphere = Sphere(radius=1.0, ntheta=20, nphi=40)
    sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")

    pb = DiffractionProblem(body=sphere, angle=1.0)
    assert len(pb.boundary_condition) == sphere.mesh.nb_faces

    with pytest.raises(TypeError):
        DiffractionProblem(boundary_conditions=[0, 0, 0])

    res = pb.make_results_container()
    assert isinstance(res, DiffractionResult)


def test_wamit_convention():
    sphere = Sphere()
    pb1 = DiffractionProblem(body=sphere, convention="Nemoh")
    pb2 = DiffractionProblem(body=sphere, convention="WAMIT")
    assert np.allclose(pb1.boundary_condition, np.conjugate(pb2.boundary_condition))


def test_radiation_problem(caplog):
    sphere = Sphere(radius=1.0, ntheta=20, nphi=40, clip_free_surface=True)

    # with pytest.raises(ValueError):
    #     RadiationProblem(body=sphere)

    sphere.add_translation_dof(direction=(0, 0, 1), name="Heave")
    pb = RadiationProblem(body=sphere)
    assert len(pb.boundary_condition) == sphere.mesh.nb_faces

    sphere.add_translation_dof(direction=(1, 0, 0), name="Surge")
    pb2 = RadiationProblem(body=sphere, radiating_dof="Heave")
    assert np.all(pb.boundary_condition == pb2.boundary_condition)

    res = pb.make_results_container()
    assert isinstance(res, RadiationResult)
    assert 'forces' not in res.__slots__
    assert res.added_masses == {}
    assert res.radiation_dampings == {}


def test_import_cal_file():
    """Test the importation of legacy Nemoh.cal files."""
    current_file_path = os.path.dirname(os.path.abspath(__file__))

    # Non symmetrical body
    cal_file_path = os.path.join(current_file_path, "Nemoh_verification_cases", "NonSymmetrical", "Nemoh.cal")
    problems = import_cal_file(cal_file_path)

    assert len(problems) == 6*41+41
    for problem in problems:
        assert problem.rho == 1000.0
        assert problem.g == 9.81
        assert problem.depth == np.infty
        assert isinstance(problem.body, FloatingBody)
        assert problem.body.nb_dofs == 6
        assert problem.body.mesh.nb_vertices == 351
        assert problem.body.mesh.nb_faces == 280
        assert problem.omega in np.linspace(0.1, 2.0, 41)
        if isinstance(problem, DiffractionProblem):
            assert problem.angle == 0.0

    # Symmetrical cylinder
    cal_file_path = os.path.join(current_file_path, "Nemoh_verification_cases", "Cylinder", "Nemoh.cal")
    problems = import_cal_file(cal_file_path)

    assert len(problems) == 6*2+2
    for problem in problems:
        assert problem.rho == 1000.0
        assert problem.g == 9.81
        assert problem.depth == np.infty
        assert isinstance(problem.body.mesh, ReflectionSymmetry)
        assert isinstance(problem.body.mesh.submeshes[0], Mesh)
        assert problem.body.nb_dofs == 6
        # assert problem.body.mesh.nb_vertices == 2*540
        assert problem.body.mesh.nb_faces == 2*300
        assert problem.omega == 0.1 or problem.omega == 2.0
        if isinstance(problem, DiffractionProblem):
            assert problem.angle == 0.0


