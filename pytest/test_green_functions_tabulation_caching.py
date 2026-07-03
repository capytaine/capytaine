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
import os
import logging
from glob import glob

import pytest

import capytaine as cpt
from capytaine.tools.cache_on_disk import cache_directory


def test_tabulation_caching_write_files():
    dir = cache_directory()
    cpt.Delhommeau(tabulation_nr=3, tabulation_nz=3)
    assert len(glob(os.path.join(dir, "tabulation*.npz"))) > 0


def test_tabulation_caching_write_files(caplog):
    dir = cache_directory()
    filename = cpt.Delhommeau()._create_or_load_tabulation(
        tabulation_nr=3,
        tabulation_rmax=1.0,
        tabulation_nz=3,
        tabulation_zmin=-1.0,
        tabulation_nb_integration_points=10,
        tabulation_cache_dir=dir
    )
    os.remove(os.path.join(dir, filename))
    with open(os.path.join(dir, filename), 'wb') as f:
        f.write(b"foooooo")
        ## Corrupting file
    with caplog.at_level(logging.WARNING):
        filename = cpt.Delhommeau()._create_or_load_tabulation(
            tabulation_nr=3,
            tabulation_rmax=1.0,
            tabulation_nz=3,
            tabulation_zmin=-1.0,
            tabulation_nb_integration_points=10,
            tabulation_cache_dir=dir
        )
    assert 'Error loading tabulation' in caplog.text
