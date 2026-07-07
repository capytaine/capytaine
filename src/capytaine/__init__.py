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

from .__about__ import (
    __title__, __description__, __version__, __author__, __uri__, __license__, __build_info__
)

from capytaine.meshes.meshes import Mesh
from capytaine.meshes.io import load_mesh
from capytaine.meshes.symmetric_meshes import ReflectionSymmetricMesh, RotationSymmetricMesh

from capytaine.meshes.predefined.cylinders import mesh_disk, mesh_horizontal_cylinder, mesh_vertical_cylinder
from capytaine.meshes.predefined.spheres import mesh_sphere
from capytaine.meshes.predefined.rectangles import mesh_rectangle, mesh_parallelepiped

from capytaine.bodies.bodies import FloatingBody
from capytaine.bodies.multibodies import Multibody
from capytaine.bodies.dofs import rigid_body_dofs

from capytaine.bem.problems_and_results import RadiationProblem, DiffractionProblem
from capytaine.bem.solver import BEMSolver
from capytaine.bem.engines import DefaultMatrixEngine
from capytaine.green_functions.delhommeau import Delhommeau, XieDelhommeau
from capytaine.green_functions.hams import LiangWuNoblesseGF, FinGreen3D, HAMS_GF

from capytaine.io.xarray import assemble_dataframe, assemble_dataset, assemble_matrices, compute_hydrostatics_dataset, export_dataset

from capytaine.ui.rich import set_logging

set_logging(level="WARNING")
