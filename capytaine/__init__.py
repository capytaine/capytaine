# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

from .__about__ import (
    __title__, __description__, __version__, __author__, __uri__, __license__
)

from capytaine.meshes.geometry import Axis, Plane, xOz_Plane, yOz_Plane, xOy_Plane
from capytaine.meshes.meshes import Mesh
from capytaine.meshes.collections import CollectionOfMeshes
from capytaine.meshes.symmetric import ReflectionSymmetricMesh, TranslationalSymmetricMesh, AxialSymmetricMesh
from capytaine.meshes.predefined.cylinders import mesh_disk, mesh_horizontal_cylinder, mesh_vertical_cylinder
from capytaine.meshes.predefined.spheres import mesh_sphere
from capytaine.meshes.predefined.rectangles import mesh_rectangle, mesh_parallelepiped

from capytaine.bodies.bodies import FloatingBody
from capytaine.bodies.dofs import rigid_body_dofs

from capytaine.bodies.predefined.spheres import Sphere
from capytaine.bodies.predefined.cylinders import VerticalCylinder, HorizontalCylinder, Disk
from capytaine.bodies.predefined.rectangles import Rectangle, RectangularParallelepiped, OpenRectangularParallelepiped

from capytaine.bem.problems_and_results import RadiationProblem, DiffractionProblem
from capytaine.bem.solver import BEMSolver
from capytaine.bem.engines import BasicMatrixEngine, HierarchicalToeplitzMatrixEngine, HierarchicalPrecondMatrixEngine
from capytaine.green_functions.delhommeau import Delhommeau, XieDelhommeau

from capytaine.post_pro.free_surfaces import FreeSurface

from capytaine.io.mesh_loaders import load_mesh
from capytaine.io.xarray import assemble_dataset

from capytaine.ui.rich import set_logging

set_logging(level="WARNING")
