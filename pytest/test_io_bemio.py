import os
import shutil
import pytest
import numpy as np
import capytaine as cpt


def test_dataset_from_wamit_via_bemio():
    bemio = pytest.importorskip("bemio.io.wamit", reason="Bemio not installed, test skipped.")
    current_file_path = os.path.dirname(os.path.abspath(__file__))
    out_file = os.path.join(current_file_path, "Bemio_verification_cases", "sphere.out")
    bemio_data = bemio.read(out_file)

    new_dataset = cpt.assemble_dataset(bemio_data)
    assert (np.moveaxis(bemio_data.body[0].am.all, 2, 0) * bemio_data.body[0].rho == \
        new_dataset['added_mass'].values).all()


def test_dataset_from_nemoh_via_bemio(tmp_path):
    bemio = pytest.importorskip("bemio.io.nemoh", reason="Bemio not installed, test skipped.")
    current_file_path = os.path.dirname(os.path.abspath(__file__))
    nemoh_data = os.path.join(current_file_path, "Nemoh_verification_cases", "Cylinder")
    shutil.copytree(nemoh_data, tmp_path, dirs_exist_ok=True)
    shutil.move(tmp_path / "reference_results", tmp_path / "Results")
    bemio_data = bemio.read(tmp_path)
    cpt.assemble_dataset(bemio_data)
