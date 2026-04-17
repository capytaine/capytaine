from typing import Optional
from itertools import cycle

from capytaine.meshes.visualization import PyVistaViewer, show_matplotlib as show_mesh_matplotlib
from capytaine.bodies.dofs import AbstractDof


def show_3d(body, *, backend=None, **kwargs):
    """Dispatch the 3D viewing to one of the available backends below."""
    backends_functions = {
            "pyvista": show_pyvista,
            "matplotlib": show_matplotlib,
            }
    if backend is not None:
        if backend in backends_functions:
            return backends_functions[backend](body, **kwargs)
        else:
            raise NotImplementedError(f"Backend '{backend}' is not implemented.")
    else:
        for backend in backends_functions:
            try:
                return backends_functions[backend](body, **kwargs)
            except (NotImplementedError, ImportError):
                pass
        raise NotImplementedError(f"No compatible backend found to show the body {body}.\n"
                                  "Consider installing `matplotlib` or `pyvista`.")


def show_pyvista(body, **kwargs) -> Optional["pv.Plotter"]:  # noqa: F821
    """Visualize a body using PyVista. See :class:`BodyPyVistaViewer` for parameters."""
    viewer = BodyPyVistaViewer(body, **kwargs)
    return viewer.show()


class BodyPyVistaViewer(PyVistaViewer):
    """Extends :class:`~capytaine.meshes.visualization.PyVistaViewer` with DOF visualization."""

    def __init__(self, body, **kwargs):
        self.body = body
        self._dofs_iterator = cycle([*body.dofs.keys(), None])
        self.help_text = [
                *super().help_text[:-1],
                "  d: toggle next dof display,",
                super().help_text[-1]
        ]

        super().__init__(
                body.mesh.merged(),  # TODO: better handle symmetric meshes
                ghost_meshes=[body.lid_mesh.merged()] if body.lid_mesh is not None else None,
                **kwargs
        )
        self.plotter.add_key_event("d", self.show_next_dof)

    def _show_dof(self, dof_name):
        mini = self.body.mesh.vertices.min()
        maxi = self.body.mesh.vertices.max()
        if 'dof' in self.plotter.actors:
            self.plotter.remove_actor('dof')
        if 'dof_label' in self.plotter.actors:
            self.plotter.remove_actor('dof_label')
        if dof_name is not None:
            if isinstance(self.body.dofs[dof_name], AbstractDof):
                dof_motion = self.body.dofs[dof_name].evaluate_motion(self.body.mesh)
            else:
                dof_motion = self.body.dofs[dof_name]
            self.plotter.add_arrows(
                self.body.mesh.faces_centers,
                dof_motion,
                name='dof',
                mag=0.08 * (maxi - mini),
                color="black",
                show_scalar_bar=False,
            )
            self.plotter.add_text(
                dof_name,
                position="upper_right",
                name='dof_label',
            )

    def show_next_dof(self):
        self._show_dof(next(self._dofs_iterator))


###################################################################################################

def show_matplotlib(body, **kwargs):
    return show_mesh_matplotlib(body.mesh_including_lid, **kwargs)
