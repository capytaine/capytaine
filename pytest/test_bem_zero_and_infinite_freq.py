import pytest

import numpy as np
import xarray as xr

import capytaine as cpt
from capytaine.green_functions.abstract_green_function import GreenFunctionEvaluationError


###########
#  Setup  #
###########

def make_simple_body():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh)
    body.add_translation_dof(direction=(1, 0, 0), name="Surge")
    return body

def make_symmetric_body():
    mesh = cpt.mesh_parallelepiped(resolution=(2, 2, 2), reflection_symmetry=True).immersed_part()
    body = cpt.FloatingBody(mesh=mesh)
    body.add_translation_dof(direction=(1, 0, 0), name="Surge")
    return body

def make_lid_body():
    mesh = cpt.mesh_sphere(resolution=(4, 4)).immersed_part()
    body = cpt.FloatingBody(mesh=mesh, lid_mesh=mesh.generate_lid())
    body.add_translation_dof(direction=(1, 0, 0), name="Surge")
    return body

BODIES = [make_simple_body, make_symmetric_body, make_lid_body]

SOLVERS = [cpt.BEMSolver(method='direct'), cpt.BEMSolver(method='indirect')]


####################################################
#  Define problem with zero or infinite frequency  #
####################################################
# Check that when a problem is defined with a zero or infinite frequency, the
# other representations of the frequency are consistent.

@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_define_problem_for_zero_freq(water_depth):
    pb = cpt.RadiationProblem(
        body=make_simple_body(),
        omega=0.0,
        water_depth=water_depth,
        radiating_dof="Surge"
    )
    assert float(pb.wavelength) == np.inf


@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_define_problem_for_inf_freq(water_depth):
    pb = cpt.RadiationProblem(
        body=make_simple_body(),
        omega=np.inf,
        water_depth=water_depth,
        radiating_dof="Surge"
    )
    assert float(pb.wavelength) == 0.0


@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_define_problem_for_zero_wavelength(water_depth):
    pb = cpt.RadiationProblem(
        body=make_simple_body(),
        wavelength=0.0,
        water_depth=water_depth,
        radiating_dof="Surge"
    )
    assert float(pb.period) == 0.0
    assert float(pb.wavenumber) == np.inf
    assert float(pb.freq) == np.inf


@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_define_problem_for_inf_wavelength(water_depth):
    pb = cpt.RadiationProblem(
        body=make_simple_body(),
        wavelength=np.inf,
        water_depth=water_depth,
        radiating_dof="Surge"
    )
    assert float(pb.period) == np.inf
    assert float(pb.wavenumber) == 0.0


##############################################################
#  Solve radiation problems with zero or infinite frequency  #
##############################################################

@pytest.mark.parametrize("omega", [0.0, np.inf])
@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_radiation_bc_is_defined_at_limit_freq(omega, water_depth):
    pb = cpt.RadiationProblem(
        body=make_simple_body(),
        omega=omega,
        water_depth=water_depth,
    )
    from capytaine.tools.symbolic_multiplication import SymbolicMultiplication
    assert isinstance(pb.boundary_condition, SymbolicMultiplication)
    assert isinstance(pb.boundary_condition.value, np.ndarray)
    assert not np.any(np.isnan(pb.boundary_condition.value))

@pytest.mark.parametrize("body", BODIES)
@pytest.mark.parametrize("omega", [0.0, np.inf])
@pytest.mark.parametrize("solver", SOLVERS)
def test_solve_radiation_inf_depth(body, omega, solver):
    pb = cpt.RadiationProblem(
        body=body(),
        omega=omega,
        water_depth=np.inf
    )
    res = solver.solve(pb)
    assert isinstance(res.added_masses['Surge'], float)
    assert res.radiation_damping['Surge'] == 0.0


@pytest.mark.parametrize("body", BODIES)
@pytest.mark.parametrize("solver", SOLVERS)
def test_solve_radiation_fin_depth_infinite_frequency(body, solver):
    pb = cpt.RadiationProblem(
        body=body(),
        omega=np.inf,
        water_depth=5.0
    )
    res = solver.solve(pb)
    assert isinstance(res.added_masses['Surge'], float)
    assert res.radiation_damping['Surge'] == 0.0


@pytest.mark.parametrize("body", BODIES)
@pytest.mark.parametrize("solver", SOLVERS)
def test_solve_radiation_fin_depth_zero_frequencies(body, solver):
    pb = cpt.RadiationProblem(
        body=body(),
        omega=0.0,
        water_depth=5.0
    )
    with pytest.raises((NotImplementedError, GreenFunctionEvaluationError)):
        solver.solve(pb)

def test_no_warning_mesh_resolution_at_zero_wavelength(caplog):
    solver = SOLVERS[0]
    pb = cpt.RadiationProblem(
        body=make_simple_body(),
        wavelength=0.0,
    )
    with caplog.at_level("WARNING"):
        solver.solve(pb)
    assert "resolution " not in caplog.text


##############################################
#  Fails with correct error for diffraction  #
##############################################

@pytest.mark.parametrize("omega", [0.0, np.inf])
@pytest.mark.parametrize("water_depth", [10.0, np.inf])
@pytest.mark.parametrize("solver", SOLVERS)
def test_fail_to_solve_diffraction_problem_with_limit_freq(omega, water_depth, solver):
    pb = cpt.DiffractionProblem(
        body=make_simple_body(),
        omega=omega,
        water_depth=water_depth
    )
    with pytest.raises(ValueError, match="zero or infinite frequency"):
        solver.solve(pb)

@pytest.mark.parametrize("omega", [0.0, np.inf])
@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_diffraction_bc_is_nan_at_limit_freq(omega, water_depth):
    pb = cpt.DiffractionProblem(
        body=make_simple_body(),
        omega=omega,
        water_depth=water_depth
    )
    assert np.all(np.isnan(pb.boundary_condition))

@pytest.mark.parametrize("omega", [0.0, np.inf])
@pytest.mark.parametrize("water_depth", [10.0, np.inf])
def test_froude_krylov_is_nan_at_limit_freq(omega, water_depth):
    from capytaine.bem.airy_waves import froude_krylov_force
    pb = cpt.DiffractionProblem(
        body=make_simple_body(),
        omega=omega,
        water_depth=water_depth
    )
    assert np.all(np.isnan(froude_krylov_force(pb)['Surge']))

##########################################
#  Fill datasets with limit frequencies  #
##########################################

def test_dataset_with_limit_frequency_including_radiation_and_diffraction(caplog):
    sphere = make_simple_body()
    test_matrix = xr.Dataset(coords={
        'omega': [0.0, 1.0, np.inf],
        'wave_direction': [0.0],
        'radiating_dof': list(sphere.dofs),
    })
    solver = cpt.BEMSolver()
    with caplog.at_level("INFO"):
        ds = solver.fill_dataset(test_matrix, sphere)
    assert np.all(np.isnan(ds.diffraction_force.sel(omega=0.0)))
    assert not np.any(np.isnan(ds.diffraction_force.sel(omega=1.0)))
    assert np.all(np.isnan(ds.diffraction_force.sel(omega=np.inf)))
    assert not np.any(np.isnan(ds.added_mass))
    assert not np.any(np.isnan(ds.radiation_damping))
    assert "Diffraction problems at zero or infinite frequency are not defined" in caplog.text
