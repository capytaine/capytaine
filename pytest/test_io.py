
import xarray as xr
import numpy as np

from capytaine.io.xarray import separate_complex_values, merge_complex_values


def test_remove_complex_values():
    original_dataset = xr.Dataset(
        data_vars={
            'variable1': (['x', 'y'], np.random.rand(4, 5)),
            'variable2': (['x', 'y', 'z'], np.random.rand(4, 5, 3) + 1j * np.random.rand(4, 5, 3)),
            'variable3': (['y', 'z'], np.random.rand(5, 3) + 1j * np.random.rand(5, 3)),
        },
        coords={
            'x': (['x'], np.linspace(0, 10, 4))
        })

    real_dataset = separate_complex_values(original_dataset)
    assert np.allclose(real_dataset['variable3'].sel(complex='re').data,
                       np.real(original_dataset['variable3'].data),
                       atol=0.1)
    assert set(original_dataset.dims) == set(real_dataset.dims) - {'complex'}

    complex_dataset = merge_complex_values(real_dataset)
    assert set(original_dataset.dims) == set(complex_dataset.dims)

