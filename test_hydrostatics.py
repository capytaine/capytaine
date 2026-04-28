import logging
from typing import Union, Sequence
import numpy as np
import capytaine as cpt
import xarray as xr
from capytaine.io.xarray import compute_hydrostatics_dataset
from capytaine.bodies.abstract_bodies import AbstractBody

LOG = logging.getLogger(__name__)


mesh = cpt.mesh_sphere().immersed_part()
body_1 = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.2)), center_of_mass=(0, 0, -0.2), name="body_1")
body_2 = cpt.FloatingBody(mesh=mesh.translated_x(5.0), dofs=cpt.rigid_body_dofs(rotation_center=(5.0, 0, -0.2)), center_of_mass=(5.0, 0, -0.2), name="body_2")
hs = compute_hydrostatics_dataset(body_1)
assert not np.any(hs.influenced_dof.str.startswith("body_"))
hs = compute_hydrostatics_dataset(body_1, rho=[1000.0, 1025.0], g=9.81)
hs = compute_hydrostatics_dataset(cpt.Multibody([body_1]))
assert np.all(hs.influenced_dof.str.startswith("body_"))
hs = compute_hydrostatics_dataset(body_1 + body_2)
assert np.all(hs.influenced_dof.str.startswith("body_"))
hs = compute_hydrostatics_dataset(body_1 + body_2, rho=[1000.0, 1025.0], g=9.81)

# mesh = cpt.mesh_sphere().immersed_part()
# body_1 = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(only=["Heave"], rotation_center=(0, 0, -0.2)), center_of_mass=(0, 0, -0.2), name="body_1")
# hs = compute_hydrostatics_dataset(body_1)
# hs = compute_hydrostatics_dataset(body_1, rho=[1000.0, 1025.0], g=9.81)

# mesh = cpt.mesh_sphere().immersed_part()
# body_1 = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.2)), name="body_1")
# hs = compute_hydrostatics_dataset(body_1)
# hs = compute_hydrostatics_dataset(body_1, rho=[1000.0, 1025.0], g=9.81)

# mesh = cpt.mesh_sphere().immersed_part()
# body_1 = cpt.FloatingBody(mesh=mesh, center_of_mass=(0, 0, -0.2), name="body_1")
# body_2 = cpt.FloatingBody(mesh=mesh.translated_x(5.0), dofs=cpt.rigid_body_dofs(rotation_center=(5.0, 0, -0.2)), center_of_mass=(5.0, 0, -0.2), name="body_2")
# hs = compute_hydrostatics_dataset(body_1)
# hs = compute_hydrostatics_dataset(body_1, rho=[1000.0, 1025.0], g=9.81)
# hs = compute_hydrostatics_dataset(body_1 + body_2)
