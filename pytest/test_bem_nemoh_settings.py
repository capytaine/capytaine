import pytest

import numpy as np
import xarray as xr
from numpy import pi

from capytaine import __version__
from capytaine.bem.nemoh import Nemoh
from capytaine.bem.problems_and_results import RadiationProblem
from capytaine.bodies.predefined.spheres import Sphere


def test_limit_freq():
    solver = Nemoh(linear_solver='gmres', use_symmetries=False, matrix_cache_size=0)

    sphere = Sphere(radius=1.0, ntheta=3, nphi=12, clip_free_surface=True)
    sphere.add_translation_dof(direction=(0, 0, 1), name="Surge")
    solver.solve(RadiationProblem(body=sphere, omega=0.0, sea_bottom=-np.infty))

    with pytest.raises(NotImplementedError):
        solver.solve(RadiationProblem(body=sphere, omega=0.0, sea_bottom=-1.0))

    solver.solve(RadiationProblem(body=sphere, omega=np.infty, sea_bottom=-np.infty))

    with pytest.raises(NotImplementedError):
        solver.solve(RadiationProblem(body=sphere, omega=np.infty, sea_bottom=-10))


def test_fill_dataset_with_kochin():
    body = Sphere(clip_free_surface=True)
    body.add_translation_dof(name="Heave")
    test_matrix = xr.Dataset(coords={
        'omega': [1.0],
        'theta': np.linspace(0, 2*pi, 5),
        'radiating_dof': list(body.dofs.keys()),
    })
    ds = Nemoh().fill_dataset(test_matrix, [body])
    print(ds)
    assert 'start_of_computation' in ds.attrs
    assert 'cache_rankine_matrices' in ds.attrs
    assert ds.attrs['capytaine_version'] == __version__


    # # Test various things on dataset.
    # assert 'Froude_Krylov_force' in data
    # wavenumbers = wavenumber_data_array(results)
    # assert isinstance(wavenumbers, xr.DataArray)
    #
    # naked_data = data.drop(["added_mass", "radiation_damping", "diffraction_force", "Froude_Krylov_force"])
    # recomputed_data = solver.fill_dataset(naked_data, [both])
    # assert "added_mass" in recomputed_data
    # assert np.allclose(recomputed_data["added_mass"].data, data["added_mass"].data)

