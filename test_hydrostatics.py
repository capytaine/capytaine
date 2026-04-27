import logging
from typing import Union, Sequence
import numpy as np
import capytaine as cpt
import xarray as xr
from capytaine.io.xarray import _squeeze_dimensions
from capytaine.bodies.dofs import DofOnSubmesh
from capytaine.bodies.abstract_bodies import AbstractBody

LOG = logging.getLogger(__name__)

def rotation_center(body: cpt.FloatingBody) -> tuple[float, float, float]:
    centers = set()
    for dof_name, dof in body.dofs.items():
        if isinstance(dof, DofOnSubmesh) and hasattr(dof.dof, "rotation_center"):
            centers.add(tuple(dof.dof.rotation_center))
        if hasattr(dof, "rotation_center"):
            centers.add(tuple(dof.rotation_center))
    if len(centers) == 0:
        return (np.nan, np.nan, np.nan)
    if len(centers) > 1:
        LOG.warning(f"Body {body.name} has no uniquely defined rotation center. Returning an arbitrary one Found {centers}. Returning an arbitrary one.")
    return next(iter(centers))


def _compute_hydrostatics_dataset(mb: cpt.Multibody, *, g: float = 9.81, rho: float = 1000.0):
    hs = xr.Dataset()
    hs.coords["space_coordinate"] = xr.DataArray(["x", "y", "z"], dims=["space_coordinate"])
    hs.coords["g"] = xr.DataArray([g], dims=["g"])
    hs.coords["rho"] = xr.DataArray([rho], dims=["rho"])
    hs.coords["body"] = xr.DataArray([b.name for b in mb.bodies], dims=["body"])

    # Putting them first such that we have an helpful error message if the dofs or the center of mass is not defined
    hs.coords["influenced_dof"] = xr.DataArray(list(mb.dofs.keys()), dims=["influenced_dof"])
    hs.coords["radiating_dof"] = xr.DataArray(list(mb.dofs.keys()), dims=["radiating_dof"])
    hs["hydrostatic_stiffness"] = xr.DataArray([[mb.compute_hydrostatic_stiffness(rho=rho, g=g)]], dims=["g", "rho", "influenced_dof", "radiating_dof"])
    hs["inertia_matrix"] = xr.DataArray([mb.compute_rigid_body_inertia(rho=rho)], dims=["rho", "influenced_dof", "radiating_dof"])

    hs.coords["center_of_mass"] = xr.DataArray(list(mb.center_of_mass.values()), dims=("body", "space_coordinate"))
    rotation_centers = np.array([rotation_center(b) for b in mb.bodies])
    if not np.all(np.isnan(rotation_centers)):
        hs.coords["rotation_center"] = xr.DataArray(rotation_centers, dims=("body", "space_coordinate"))
    else:
        LOG.debug("No rotation center could be found for any body. Skipping the `rotation_center` coordinate.")
    hs["center_of_buoyancy"] = xr.DataArray(list(mb.center_of_buoyancy.values()), dims=("body", "space_coordinate"))
    hs["disp_mass"] = xr.DataArray([[b.disp_mass(rho=rho) for b in mb.bodies]], dims=("rho", "body"), attrs=dict(long_name="Diplaced mass", units="kg"))
    hs["draught"] = xr.DataArray([np.abs(b.mesh.vertices[:, 2].min()) for b in mb.bodies], dims=("body"))
    return hs


def compute_hydrostatics_dataset(
        body: AbstractBody,
        *,
        g: Union[float, Sequence[float]] = 9.81,
        rho: Union[float, Sequence[float]] = 1000.0
        ):

    if isinstance(g, (int, float)):
        g = [g]
    if isinstance(rho, (int, float)):
        rho = [rho]
    if isinstance(body, cpt.FloatingBody):
        body = cpt.Multibody([body])

    datasets = []
    for g_ in g:
        for rho_ in rho:
            datasets.append(_compute_hydrostatics_dataset(body, g=g_, rho=rho_))
    hs = xr.merge(datasets, compat="override", join="outer")

    optional_dims = ["g", "rho", "body"]
    hs = _squeeze_dimensions(hs, dimensions=optional_dims)
    return hs


# mesh = cpt.mesh_sphere().immersed_part()
# body_1 = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.2)), center_of_mass=(0, 0, -0.2), name="body_1")
# body_2 = cpt.FloatingBody(mesh=mesh.translated_x(5.0), dofs=cpt.rigid_body_dofs(rotation_center=(5.0, 0, -0.2)), center_of_mass=(5.0, 0, -0.2), name="body_2")
# hs = compute_hydrostatics_dataset(body_1)
# hs = compute_hydrostatics_dataset(body_1, rho=[1000.0, 1025.0], g=9.81)
# hs = compute_hydrostatics_dataset(body_1 + body_2)
# hs = compute_hydrostatics_dataset(body_1 + body_2, rho=[1000.0, 1025.0], g=9.81)

# mesh = cpt.mesh_sphere().immersed_part()
# body_1 = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(only=["Heave"], rotation_center=(0, 0, -0.2)), center_of_mass=(0, 0, -0.2), name="body_1")
# hs = compute_hydrostatics_dataset(body_1)
# hs = compute_hydrostatics_dataset(body_1, rho=[1000.0, 1025.0], g=9.81)

mesh = cpt.mesh_sphere().immersed_part()
body_1 = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, -0.2)), name="body_1")
hs = compute_hydrostatics_dataset(body_1)
hs = compute_hydrostatics_dataset(body_1, rho=[1000.0, 1025.0], g=9.81)

mesh = cpt.mesh_sphere().immersed_part()
body_1 = cpt.FloatingBody(mesh=mesh, center_of_mass=(0, 0, -0.2), name="body_1")
body_2 = cpt.FloatingBody(mesh=mesh.translated_x(5.0), dofs=cpt.rigid_body_dofs(rotation_center=(5.0, 0, -0.2)), center_of_mass=(5.0, 0, -0.2), name="body_2")
hs = compute_hydrostatics_dataset(body_1)
hs = compute_hydrostatics_dataset(body_1, rho=[1000.0, 1025.0], g=9.81)
hs = compute_hydrostatics_dataset(body_1 + body_2)
