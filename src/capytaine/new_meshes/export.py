import numpy as np
from typing import List

from capytaine.tools.optional_imports import import_optional_dependency

def mesh_to_pyvista(
    vertices: np.ndarray, faces: List[List[int]]
) -> "pv.UnstructuredGrid":
    """
    Build a PyVista UnstructuredGrid from a list of faces (triangles or quads).
    """
    pv = import_optional_dependency("pyvista")

    # flatten into the VTK cell‐array format: [n0, i0, i1, ..., in-1, n1, j0, j1, ...]
    flat_cells = []
    cell_types = []
    for face in faces:
        n = len(face)
        flat_cells.append(n)
        flat_cells.extend(face)
        if n == 3:
            cell_types.append(pv.CellType.TRIANGLE)
        elif n == 4:
            cell_types.append(pv.CellType.QUAD)
        else:
            # if you ever have ngons, you can map them as POLYGON:
            cell_types.append(pv.CellType.POLYGON)

    cells_array = np.array(flat_cells, dtype=np.int64)
    cell_types = np.array(cell_types, dtype=np.uint8)

    return pv.UnstructuredGrid(cells_array, cell_types, vertices.astype(np.float32))
