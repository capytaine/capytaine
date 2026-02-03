from typing import Optional
from itertools import cycle

import numpy as np

from capytaine.tools.optional_imports import import_optional_dependency
from capytaine.new_meshes.visualization import show_pyvista as show_pyvista_mesh

def show_pyvista(
    body,
    *,
    plotter=None,
    **kwargs
) -> Optional["pv.Plotter"]:  # noqa: F821
    pv = import_optional_dependency("pyvista")

    if plotter is None:
        default_plotter = True
        plotter = pv.Plotter()
    else:
        default_plotter = False

    plotter = show_pyvista_mesh(body.mesh, plotter=plotter, **kwargs)

    def show_dof(dof_name):
        mini = body.mesh.vertices.min()
        maxi = body.mesh.vertices.max()
        if 'dof' in plotter.actors:
            plotter.remove_actor('dof')
        if 'dof_label' in plotter.actors:
            plotter.remove_actor('dof_label')
        if dof_name is not None:
            plotter.add_arrows(
                    body.mesh.faces_centers,
                    body.dofs[dof_name],
                    name='dof',
                    mag=0.04*(maxi-mini),
                    show_scalar_bar=False
                    )
            plotter.add_text(
                    dof_name,
                    position="upper_right",
                    name='dof_label',
                    )

    dofs_iterators = cycle([*body.dofs.keys(), None])
    def show_next_dof():
        show_dof(next(dofs_iterators))

    plotter.add_key_event("d", lambda : show_next_dof())

    if default_plotter:
        plotter.show()
    else:
        return plotter


if __name__ == "__main__":
    import capytaine as cpt
    from capytaine.new_meshes.meshes import to_new_mesh
    mesh = to_new_mesh(cpt.mesh_horizontal_cylinder())
    body = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs())
    import pyvista as pv
    plotter = pv.Plotter()
    plotter = show_pyvista(body, plotter=plotter)
    plotter.show()
