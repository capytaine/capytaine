# Copyright 2025 Mews Labs
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

from typing import Union, Any, Optional
from io import IOBase
from pathlib import Path
import tempfile
import logging

import numpy as np
import xarray as xr

from capytaine.tools.optional_imports import silently_import_optional_dependency
from capytaine.new_meshes import Mesh, ReflectionSymmetricMesh

meshio = silently_import_optional_dependency("meshio")
trimesh = silently_import_optional_dependency("trimesh")

LOG = logging.getLogger(__name__)

_MESHIO_EXTS = {"inp", "msh", "avs", "cgns", "xml", "e", "exo", "f3grid", "h5m", "mdpa", "mesh", "meshb", "med", "bdf", "fem", "nas", "vol", "vol.gz", "obj", "off", "post", "post.gz", "dato", "dato.gz", "ply", "stl", "dat", "node", "ele", "svg", "su2", "ugrid", "vtk", "vtu", "wkt", "xdmf", "xmf", "gmsh"}

_TRIMESH_EXTS = {"obj", "stl", "ply", "glb", "gltf", "off"}

_XARRAY_EXTS = {"nc", "netcdf"}

_BUILTIN_EXTS = {"pnl", "hst", "mar", "gdf", "nemoh", "wamit", "hydrostar", "hams"}

_ALL_EXTS = _MESHIO_EXTS | _TRIMESH_EXTS | _XARRAY_EXTS | _BUILTIN_EXTS

# Structure of the function calls:
#
#                                  load_mesh
#                                  /     | \
#                         if path /      |  \
#                                /       |   \
#                  _load_from_path       |    --------------------\
#                    |         \         |                        |
#                    |          \        | if opened file         | if object from external lib
#   if implemented   |           \       |                        |
#   by external lib  |         _read_mesh_from_file_like_object   |
#                    |                   |  \                     |
#                    |   if implemented  |   \ if built-in        |
#                    |   by external lib |    \                   |
#                    \                   |    _read_***           |
#                     \                  |                        |
#                      \                 |        ----------------/
#                       \                |       /
#                      _import_from_***_mesh_class
#


def load_mesh(mesh_to_be_loaded, file_format=None, *, backend=None) -> Union[Mesh, ReflectionSymmetricMesh]:
    """Load a mesh from a file or file-like object.

    This function can load mesh data from many file formats. It supports file
    paths and file-like objects.

    Parameters
    ----------
    mesh_to_be_loaded : str, pathlib.Path, or file-like object
        Can be either:
            - a path to a mesh file
            - a file-like object with mesh file data
            - a mesh object from one of the compatible external libraries (meshio, trimesh)
    file_format : str, optional
        Format hint used only when loading from file-like objects, since the
        filename extension is unavailable.
        Valid values include ``"stl"``, ``"obj"``, ``"hst"``, ``"pnl"``, etc.
        Can be automatically inferred from file extension when a path is provided as first argument.

    Returns
    -------
    Mesh or ReflectionSymmetricMesh
        A Mesh object, or a ReflectionSymmetricMesh if symmetry information
        is present in the file format.

    Raises
    ------
    ValueError
        If ``file_format`` is not provided when loading from a file-like
        object, or if an unsupported format is encountered.

    Examples
    --------
    Load a mesh from a file path.

    >>> load_mesh("model.stl")

    Load from a gzip-compressed file.

    >>> import gzip
    >>> with gzip.open("model.stl.gz") as handler:
    ...     load_mesh(handler, file_format="stl")

    Load from any file-like object.

    >>> with open("model.obj", "rb") as handler:
    ...     load_mesh(handler, file_format="obj")

    """
    if isinstance(mesh_to_be_loaded, IOBase):  # A file already opened
        return _read_mesh_from_file_like_object(mesh_to_be_loaded, file_format, backend=backend)

    elif backend in {None, "xarray"} and isinstance(mesh_to_be_loaded, xr.Dataset):
        return _import_from_xarray_dataset(mesh_to_be_loaded)

    elif trimesh is not None and backend in {None, "trimesh"} and isinstance(mesh_to_be_loaded, trimesh.base.Trimesh):
        return _import_from_trimesh_mesh_class(mesh_to_be_loaded)

    elif meshio is not None and backend in {None, "meshio"} and isinstance(mesh_to_be_loaded, meshio.Mesh):
        return _import_from_meshio_mesh_class(mesh_to_be_loaded)

    elif isinstance(mesh_to_be_loaded, (str, Path)):
        return _load_from_path(Path(mesh_to_be_loaded), file_format, backend=backend)

    else:
        raise TypeError(f"load_mesh() can't interpret the type of input of {repr(mesh_to_be_loaded)}"
                        + "" if backend is None else f" when using backend={backend}")



def _normalise_format_hint(value: Any) -> str:
    text = str(value).strip().lower()
    if text.startswith("."):
        text = text[1:]
    return text


def _load_from_path(path: Path, file_format: Optional[str] = None, *, backend: Optional[str] = None):
    if not path.exists():
        raise FileNotFoundError(f"Mesh file not found: {path}")

    if file_format is None:
        if len(path.suffixes) == 1:
            file_format = path.suffixes[0]
        else:
            raise ValueError(f"No file format has been provided, nor can it be inferred from path: {path}")
    fmt = _normalise_format_hint(file_format)

    if trimesh is not None and backend in {None, "trimesh"} and fmt in _TRIMESH_EXTS:
        trimesh_mesh = trimesh.load(path, force="mesh", file_type=fmt)
        return _import_from_trimesh_mesh_class(trimesh_mesh)
        # We could skip this, as _read_mesh_from_file_like_object would also try to pass the opened file to trimesh.
        # However, some file format might supported by trimesh require the file to be opened in byte mode, and we don't want to handle that.

    if meshio is not None and backend in {None, "meshio"} and fmt in _MESHIO_EXTS:
        meshio_mesh = meshio.read(path, file_format=fmt)
        return _import_from_meshio_mesh_class(meshio_mesh)
        # We could skip this, as _read_mesh_from_file_like_object would also try to pass the opened file to meshio.
        # But meshio does not support reading from an opened file, so we avoid the need for a temporary file by loading from the path here.

    if backend in {None, "xarray"} and fmt in _XARRAY_EXTS:
        dataset = xr.load_dataset(path)
        return _import_from_xarray_dataset(dataset)

    with open(path, 'r') as f:
        return _read_mesh_from_file_like_object(f, file_format)


def _read_mesh_from_file_like_object(
    file_obj: IOBase,
    file_format: str,
    *,
    backend: Optional[str] = None
) -> Union[Mesh, ReflectionSymmetricMesh]:
    fmt = _normalise_format_hint(file_format)

    if fmt in _BUILTIN_EXTS and backend is None:
        if fmt == "mar" or fmt == "nemoh":
            return _read_mar(file_obj)
        elif fmt == "gdf" or fmt == "wamit":
            return _read_gdf(file_obj)
        elif fmt == "hst" or fmt == "hydrostar":
            return _read_hst(file_obj)
        elif fmt == "pnl" or fmt == "hams":
            return _read_pnl(file_obj)
    elif trimesh is not None and backend in {None, "trimesh"} and fmt in _TRIMESH_EXTS:
        trimesh_mesh = trimesh.load(file_obj, force="mesh", file_type=fmt)
        return _import_from_trimesh_mesh_class(trimesh_mesh)
    elif meshio is not None and backend in {None, "meshio"} and fmt in _MESHIO_EXTS:
        # Meshio does not support reading from opened file, so it is written in a temporary file here
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            content = file_obj.read()
            if isinstance(content, str):
                temp_file.write(content.encode("utf-8"))
            else:
                temp_file.write(content)
        meshio_mesh = meshio.read(temp_file.name, file_format=fmt)
        return _import_from_meshio_mesh_class(meshio_mesh)
    elif backend in {None, "xarray"} and fmt in _XARRAY_EXTS:
        # Seems to be supported according to documentation of xarray,
        # but could not make it work with any backend in practice...
        dataset = xr.load_dataset(file_obj)
        return _import_from_xarray_dataset(dataset)
    else:
        raise ValueError(
                f"Unrecognized or unsupported mesh format: {file_format}. "
                f"Supported mesh formats (some may require external libraries to be installed): {sorted(_ALL_EXTS)}"
        )


def _import_from_meshio_mesh_class(meshio_mesh):
    faces = []
    if "quad" in meshio_mesh.cells_dict:
        faces.extend(list(meshio_mesh.cells_dict["quad"]))
    if "triangle" in meshio_mesh.cells_dict:
        faces.extend(list(meshio_mesh.cells_dict["triangle"]))

    if not faces:
        raise ValueError("No triangle or quad cells found in meshio mesh.")

    return Mesh(meshio_mesh.points, faces)


def _import_from_trimesh_mesh_class(trimesh_mesh):
    if not isinstance(trimesh_mesh, trimesh.base.Trimesh):
        raise TypeError(f"Expected trimesh.base.Trimesh, received {type(trimesh_mesh)}")
    return Mesh(trimesh_mesh.vertices, trimesh_mesh.faces)


def _import_from_xarray_dataset(dataset):
    return Mesh.from_list_of_faces(dataset["mesh_vertices"].values)


def _read_hst(file_obj: IOBase) -> Union[Mesh, ReflectionSymmetricMesh]:
    """HST files have a 1-based indexing"""

    lines = file_obj.readlines()

    optional_keywords = ['PROJECT', 'SYMMETRY']
    not_implemented_optional_keywords = ['USER', 'REFLENGTH', 'GRAVITY', 'RHO', 'NBBODY']

    vertices = []
    faces = []
    optional_data = {kw: None for kw in optional_keywords}
    current_context = None
    ignored_lines = []

    for i_line, line in enumerate(lines):
        line = line.lstrip()

        if line == '':
            continue

        elif line.startswith("COORDINATES"):
            current_context = 'vertices'

        elif current_context == 'vertices' and line.startswith("ENDCOORDINATES"):
            current_context = None

        elif line.startswith("PANEL"):
            panels_type = int(line[10:])
            current_context = ('panels', panels_type)

        elif (current_context == ('panels', 0) or current_context == ('panels', 1)) and line.startswith("ENDPANEL"):
            current_context = None

        elif current_context == 'vertices':  # parse vertex coordinates
            numbers = line.split()
            if len(numbers) == 4:
                i_vertex, x, y, z = numbers
                if int(i_vertex) != len(vertices) + 1:
                    raise ValueError(
                        f"HST mesh reader expected the next vertex to be indexed as {len(vertices)+1}, "
                        f"but it was actually indexed as {i_vertex} (line {i_line+1}.")
            elif len(numbers) == 3:
                x, y, z = numbers
            vertices.append([x, y, z])

        elif current_context == ('panels', 0):  # parse face definition (no index given)
            numbers = line.split()
            if len(numbers) == 3:
                v1, v2, v3 = numbers
                v4 = v3
            elif len(numbers) == 4:
                v1, v2, v3, v4 = numbers
            faces.append([v1, v2, v3, v4])

        elif current_context == ('panels', 1):  # parse face definition
            numbers = line.split()
            if len(numbers) == 4:
                i_face, v1, v2, v3 = numbers
                v4 = v3
            elif len(numbers) == 5:
                i_face, v1, v2, v3, v4 = numbers

            if int(i_face) != len(faces) + 1:
                ii = len(faces) + 1
                raise ValueError(f"HST mesh reader expected the next face to be indexed {ii},\n"
                                 f"but it was actually indexed with {i_face} (line {i_line+1}.")
            faces.append([v1, v2, v3, v4])

        elif line.startswith("ENDFILE"):
            break

        else:
            for keyword in optional_data:
                if line.startswith(keyword):
                    optional_data[keyword] = line[len(keyword)+1:].lstrip(':').strip()
                    break
            else:
                ignored_lines.append((i_line+1, line))

    if len(ignored_lines) > 0:
        formatted_ignored_lines = ["{: 4} | {}".format(i, line.strip('\n')) for (i, line) in ignored_lines]
        LOG.warning("HST mesh reader ignored the following lines:\n" + "\n".join(formatted_ignored_lines))

    vertices = np.array(vertices, dtype=float)
    faces = np.array(faces, dtype=int) - 1

    if optional_data['SYMMETRY'] == '1':
        return ReflectionSymmetricMesh(Mesh(vertices, faces), plane="xOz")
    elif optional_data['SYMMETRY'] == '2':
        return ReflectionSymmetricMesh(
            ReflectionSymmetricMesh(Mesh(vertices, faces), plane="yOz"),
            plane="xOz",
        )
    else:
        return Mesh(vertices, faces)


def _read_gdf(file_obj: IOBase) -> Union[Mesh, ReflectionSymmetricMesh]:
    title = file_obj.readline()
    ulen, grav = map(float, file_obj.readline().split()[:2])
    isx, isy = map(int, file_obj.readline().split()[:2])
    npan = int(file_obj.readline().split()[0])
    faces_vertices = np.genfromtxt(file_obj)

    faces_vertices = faces_vertices.reshape(-1, 3)
    vertices, indices = np.unique(faces_vertices, axis=0, return_inverse=True)
    faces = indices.reshape(-1, 4)

    if faces.shape[0] != npan:
        raise ValueError(
            f"In GDF file, npan value: {npan} is not equal to face count: \
                {faces.shape[0]}."
        )

    if isx == 1 and isy == 1:
        return ReflectionSymmetricMesh(
            ReflectionSymmetricMesh(Mesh(vertices, faces), plane="xOz"),
            plane="yOz",
        )
    elif isx == 1:
        return ReflectionSymmetricMesh(Mesh(vertices, faces), plane="yOz")
    elif isy == 1:
        return ReflectionSymmetricMesh(Mesh(vertices, faces), plane="xOz")
    else:
        return Mesh(vertices, faces)


def _read_mar(file_obj) -> Union[Mesh, ReflectionSymmetricMesh]:
    vertices = []
    faces = []

    # Read header: "n_faces symmetry_flag"
    header = file_obj.readline().split()
    symmetry_flag = int(header[1]) if len(header) > 1 else 0

    # Read vertices until "0" marker
    for line in file_obj:
        tokens = line.split()
        if tokens[0] == "0":
            break
        vertices.append(list(map(float, tokens[1:])))

    # Read faces until "0" marker
    for line in file_obj:
        tokens = line.split()
        if tokens[0] == "0":
            break
        faces.append(list(map(int, tokens)))

    # Convert to numpy arrays and adjust indices (Fortran 1-based to Python 0-based)
    vertices = np.array(vertices, dtype=float)
    faces = np.array(faces, dtype=int) - 1

    if symmetry_flag == 1:
        return ReflectionSymmetricMesh(Mesh(vertices, faces), plane="xOz")
    else:
        return Mesh(vertices, faces)


def _read_pnl(file_obj) -> Union[Mesh, ReflectionSymmetricMesh]:
    # Skip 3 title lines
    file_obj.readline()
    file_obj.readline()
    file_obj.readline()
    # Read header data
    nb_faces, nb_vertices, x_sym, y_sym = map(int, file_obj.readline().split())
    # Skip 2 more lines
    file_obj.readline()
    file_obj.readline()
    vertices = np.genfromtxt((file_obj.readline() for _ in range(nb_vertices)), usecols=(1, 2, 3))
    # Skip 3 more lines
    file_obj.readline()
    file_obj.readline()
    file_obj.readline()

    faces = np.zeros((nb_faces, 4), dtype=int)
    for i in range(nb_faces):
        index, nb_corners, *data = map(int, file_obj.readline().split())
        assert i+1 == index
        if nb_corners == 3:  # Triangle
            assert len(data) == 3
            faces[i, 0:3] = data
            faces[i, 3] = faces[i, 2]  # Convention for triangles in Capytaine: repeat last vertex
        elif int(nb_corners) == 4:  # Quadrangle
            assert len(data) == 4
            faces[i, :] = data
    faces = faces - 1  # Going from Fortran 1-based indices to Numpy 0-based indices

    if x_sym == 1 and y_sym == 0:
        return ReflectionSymmetricMesh(Mesh(vertices, faces), plane="yOz")
    elif x_sym == 0 and y_sym == 1:
        return ReflectionSymmetricMesh(Mesh(vertices, faces), plane="xOz")
    elif x_sym == 1 and y_sym == 1:
        half_mesh = ReflectionSymmetricMesh(Mesh(vertices, faces), plane="xOz")
        return ReflectionSymmetricMesh(half_mesh, plane="yOz")
    else:
        return Mesh(vertices, faces)
