import numpy as np
import xarray as xr

from capytaine.tools.optional_imports import import_optional_dependency


def export_mesh(mesh, format: str):
    format = format.lower()
    if format == "pyvista":
        return export_to_pyvista(mesh)
    elif format == "xarray":
        return export_to_xarray(mesh)
    elif format == "meshio":
        return export_to_meshio(mesh)
    elif format == "trimesh":
        return export_to_trimesh(mesh)
    else:
        raise ValueError(f"Unrecognized output format: {format}")


def export_to_pyvista(mesh):
    """
    Build a PyVista UnstructuredGrid from a list of faces (triangles or quads).
    """
    pv = import_optional_dependency("pyvista")

    # flatten into the VTK cell‐array format: [n0, i0, i1, ..., in-1, n1, j0, j1, ...]
    flat_cells = []
    cell_types = []
    for face in mesh._faces:
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

    return pv.UnstructuredGrid(cells_array, cell_types, mesh.vertices.astype(np.float32))


def export_to_xarray(mesh):
    return xr.Dataset(
            {
                "mesh_vertices": (
                    ["face", "vertices_of_face", "space_coordinate"],
                    mesh.as_array_of_faces()
                    )
                },
            coords={
                "space_coordinate": ["x", "y", "z"],
                })


def export_to_meshio(mesh):
    meshio = import_optional_dependency("meshio")

    quads = [f for f in mesh._faces if len(f) == 4]
    tris = [f for f in mesh._faces if len(f) == 3]

    cells = []
    if quads:
        cells.append(meshio.CellBlock("quad", np.array(quads, dtype=np.int32)))
    if tris:
        cells.append(meshio.CellBlock("triangle", np.array(tris, dtype=np.int32)))

    return meshio.Mesh(points=mesh.vertices, cells=cells)


def export_to_trimesh(mesh):
    trimesh = import_optional_dependency("trimesh")
    triangle_faces = []
    for face in mesh._faces:
        if len(face) == 4 and face[3] != face[2]:
            triangle_faces.append([face[0], face[1], face[2]])
            triangle_faces.append([face[0], face[2], face[3]])
        else:
            triangle_faces.append(face[:3])
    return trimesh.Trimesh(
        vertices=mesh.vertices,
        faces=np.array(triangle_faces),
        process=False
    )
