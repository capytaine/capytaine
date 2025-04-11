"""Functions to load meshes from different file formats.
Based on meshmagick <https://github.com/LHEEA/meshmagick> by François Rongère.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin, based on the work of François Rongère
# See LICENSE file at <https://github.com/mancellin/capytaine>

import os
import logging
import numpy as np

from capytaine.meshes.meshes import Mesh
from capytaine.meshes.symmetric import ReflectionSymmetricMesh
from capytaine.meshes.geometry import xOz_Plane, yOz_Plane
from capytaine.tools.optional_imports import import_optional_dependency, silently_import_optional_dependency

LOG = logging.getLogger(__name__)

real_str = r'[+-]?(?:\d+\.\d*|\d*\.\d+)(?:[Ee][+-]?\d+)?'  # Regex for floats


def _check_file(filename, name=None):
    if not os.path.isfile(filename):
        raise IOError("file %s not found" % filename)
    return


def load_mesh(mesh, file_format=None, name=None):
    """Driver function that loads every mesh file format known by meshmagick.
    Dispatch to one of the other function depending on file_format.

    Parameters
    ----------
    mesh: str or meshio object
        Either the path to the mesh on disk
        or a meshio object to be loaded with the dedicated method
    file_format: str, optional
        format of the mesh defined in the extension_dict dictionary
    name: str, optional
        name for the created mesh object

    Returns
    -------
    Mesh or SymmetricMesh
        the loaded mesh
    """
    meshio = silently_import_optional_dependency("meshio")
    if meshio is not None and isinstance(mesh, meshio._mesh.Mesh):
        from capytaine.io.meshio import load_from_meshio
        return load_from_meshio(mesh, name=name)

    filename = mesh

    _check_file(filename)

    if file_format is None:
        _, file_format = os.path.splitext(filename)
        file_format = file_format.strip('.').lower()

    if file_format not in extension_dict:
        raise IOError('Extension ".%s" is not known' % file_format)

    loader = extension_dict[file_format]

    if name is None:
        name = filename

    return loader(filename, name)


def load_RAD(filename, name=None):
    """Loads RADIOSS mesh files. This export file format may be chosen in ICEM meshing program.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh

    Note
    ----
    RAD files have a 1-indexing
    """

    import re
    _check_file(filename)
    ifile = open(filename, 'r')
    data = ifile.read()
    ifile.close()

    # node_line = r'\s*\d+(?:\s*' + real_str + '){3}'
    node_line = r'\s*\d+\s*(' + real_str + r')\s*(' + real_str + r')\s*(' + real_str + ')'
    node_section = r'((?:' + node_line + ')+)'

    elem_line = r'^\s*(?:\d+\s+){6}\d+\s*[\r\n]+'
    elem_section = r'((?:' + elem_line + '){3,})'

    pattern_node_line = re.compile(node_line, re.MULTILINE)
    # pattern_node_line_group = re.compile(node_line, re.MULTILINE)
    pattern_elem_line = re.compile(elem_line, re.MULTILINE)
    pattern_node_section = re.compile(node_section, re.MULTILINE)
    pattern_elem_section = re.compile(elem_section, re.MULTILINE)

    vertices = []
    node_section = pattern_node_section.search(data).group(1)
    for node in pattern_node_line.finditer(node_section):
        vertices.append(list(map(float, list(node.groups()))))
    vertices = np.asarray(vertices, dtype=float)

    faces = []
    elem_section = pattern_elem_section.search(data).group(1)
    for elem in pattern_elem_line.findall(elem_section):
        faces.append(list(map(int, elem.strip().split()[3:])))
    faces = np.asarray(faces, dtype=int) - 1

    return Mesh(vertices, faces, name)


def load_HST(filename, name=None):
    """Loads HYDROSTAR (Bureau Veritas (c)) mesh files.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh

    Note
    ----
    HST files have a 1-indexing
    """

    _check_file(filename)

    with open(filename, 'r') as f:
        lines = f.readlines()

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
                        f"but it was actually indexed as {i_vertex} (line {i_line+1} of {filename}).")
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
                                 f"but it was actually indexed with {i_face} (line {i_line+1} of file {filename}).")
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
        LOG.warning(f"HST mesh reader ignored the following lines from file {filename}:\n" + "\n".join(formatted_ignored_lines))

    vertices = np.array(vertices, dtype=float)
    faces = np.array(faces, dtype=int) - 1

    if name is None: name = optional_data['PROJECT']

    if optional_data['SYMMETRY'] == '1':
        return ReflectionSymmetricMesh(Mesh(vertices, faces, f"half_of_{name}"), xOz_Plane, name)
    elif optional_data['SYMMETRY'] == '2':
        return ReflectionSymmetricMesh(ReflectionSymmetricMesh(Mesh(vertices, faces, f"quarter_of_{name}"), yOz_Plane, f"half_of_{name}"), xOz_Plane, name)
    else:
        return Mesh(vertices, faces, name)


def load_DAT(filename, name=None):
    """Not implemented.
    Intended to load .DAT files used in DIODORE (PRINCIPIA (c))
    """
    _check_file(filename)
    raise NotImplementedError


def load_INP(filename, name=None):
    """Loads DIODORE (PRINCIPIA (c)) configuration file format.

    It parses the .INP file and extracts meshes defined in subsequent .DAT files using the different information
    contained in the .INP file.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh

    Note
    ----
    INP/DAT files use a 1-indexing
    """
    _check_file(filename)
    import re

    with open(filename, 'r') as f:
        text = f.read()

    # Retrieving frames into a dictionary frames
    pattern_frame_str = r'^\s*\*FRAME,NAME=(.+)[\r\n]+(.*)'
    pattern_frame = re.compile(pattern_frame_str, re.MULTILINE)

    frames = {}
    for match in pattern_frame.finditer(text):
        frame_name = match.group(1).strip()
        frame_vector = re.split(r'[, ]', match.group(2).strip())
        frames[frame_name] = np.asarray(list(map(float, frame_vector)))

    # Storing the inp layout into a list of dictionary
    pattern_node_elements = re.compile(r'^\s*\*(NODE|ELEMENT),(.*)', re.MULTILINE)
    layout = []
    mesh_files = {}
    for match in pattern_node_elements.finditer(text):
        field_dict = dict()
        field_dict['type'] = match.group(1)
        if field_dict['type'] == 'NODE':
            field_dict['INCREMENT'] = 'NO'
        opts = match.group(2).split(',')
        for opt in opts:
            key, pair = opt.split('=')
            field_dict[key] = pair.strip()

        # Retrieving information on mesh files and their usage
        file = field_dict['INPUT']
        if file in mesh_files:
            mesh_files[file][field_dict['type'] + '_CALL_INP'] += 1
        else:
            mesh_files[file] = {}
            mesh_files[file]['NODE_CALL_INP'] = 0
            mesh_files[file]['ELEMENT_CALL_INP'] = 0
            mesh_files[file][field_dict['type'] + '_CALL_INP'] += 1

        layout.append(field_dict)

        # RETRIEVING DATA SECTIONS FROM MESHFILES
        # patterns for recognition of sections
    node_line = r'\s*\d+(?:\s+' + real_str + '){3}'
    node_section = r'((?:' + node_line + ')+)'
    elem_line = r'^ +\d+(?: +\d+){3,4}[\r\n]+'  # 3 -> triangle, 4 -> quadrangle
    elem_section = r'((?:' + elem_line + ')+)'
    pattern_node_line = re.compile(node_line, re.MULTILINE)
    pattern_elem_line = re.compile(elem_line, re.MULTILINE)
    pattern_node_section = re.compile(node_section, re.MULTILINE)
    pattern_elem_section = re.compile(elem_section, re.MULTILINE)

    for file in mesh_files:
        try:
            meshfile = open(os.path.join(os.path.dirname(filename), file + '.DAT'), 'r')
        except:
            raise IOError('File {0:s} not found'.format(file + '.DAT'))
        data = meshfile.read()
        meshfile.close()

        node_section = pattern_node_section.findall(data)
        if len(node_section) > 1:
            raise IOError("""Several NODE sections into a .DAT file is not supported by meshmagick
                              as it is considered as bad practice""")
        node_array = []
        idx_array = []
        for node in pattern_node_line.findall(node_section[0]):
            node = node.split()

            node[0] = int(node[0])
            idx_array.append(node[0])
            node[1:] = list(map(float, node[1:]))
            node_array.append(node[1:])

        mesh_files[file]['NODE_SECTION'] = node_array

        # Detecting renumberings to do
        real_idx = 0
        # renumberings = []
        id_new = - np.ones(max(idx_array) + 1, dtype=int)
        # FIXME: cette partie est tres buggee !!!
        for i, idx in enumerate(idx_array):
            id_new[idx] = i+1

        mesh_files[file]['ELEM_SECTIONS'] = []
        for elem_section in pattern_elem_section.findall(data):

            elem_array = []
            for elem in pattern_elem_line.findall(elem_section):
                elem = list(map(int, elem.split()))
                # for node in elem[1:]:
                elem = id_new[elem[1:]].tolist()
                if len(elem) == 3:  # Case of a triangle, we repeat the first node at the last position
                    elem.append(elem[0])

                elem_array.append(list(map(int, elem)))
            mesh_files[file]['ELEM_SECTIONS'].append(elem_array)
        mesh_files[file]['nb_elem_sections'] = len(mesh_files[file]['ELEM_SECTIONS'])

        mesh_files[file]['nb_elem_sections_used'] = 0

    nb_nodes = 0
    nb_elems = 0
    for field in layout:
        file = field['INPUT']
        if field['type'] == 'NODE':
            nodes = np.asarray(mesh_files[file]['NODE_SECTION'], dtype=float)
            # Translation of nodes according to frame option id any
            nodes += frames[field['FRAME']]  # TODO: s'assurer que frame est une options obligatoire...

            if nb_nodes == 0:
                vertices = nodes.copy()
                nb_nodes = vertices.shape[0]
                increment = False
                continue

            if field['INCREMENT'] == 'NO':
                vertices[idx, :] = nodes.copy()
                increment = False
            else:
                vertices = np.concatenate((vertices, nodes))
                nb_nodes = vertices.shape[0]
                increment = True
        else:  # this is an ELEMENT section
            elem_section = np.asarray(mesh_files[file]['ELEM_SECTIONS'][mesh_files[file]['nb_elem_sections_used']],
                                      dtype=int)

            mesh_files[file]['nb_elem_sections_used'] += 1
            if mesh_files[file]['nb_elem_sections_used'] == mesh_files[file]['nb_elem_sections']:
                mesh_files[file]['nb_elem_sections_used'] = 0

            # Updating to new id of nodes
            elems = elem_section
            if increment:
                elems += nb_nodes

            if nb_elems == 0:
                faces = elems.copy()
                nb_elems = faces.shape[0]
                continue
            else:
                faces = np.concatenate((faces, elems))
                nb_elems = faces.shape[0]

    return Mesh(vertices, faces-1, name)


def load_TEC(filename, name=None):
    """Loads TECPLOT (Tecplot (c)) mesh files.

    It relies on the tecplot file reader from the VTK library.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh

    Note
    ----
    TEC files have a 1-indexing
    """

    import re

    _check_file(filename)

    data_pattern = re.compile(
                    r'ZONE.*\s*N\s*=\s*(\d+)\s*,\s*E=\s*(\d+)\s*,\s*F\s*=\s*FEPOINT\s*,\s*ET\s*=\s*QUADRILATERAL\s+'
                    + r'(^(?:\s*' + real_str + r'){3,})\s+'
                    + r'(^(?:\s*\d+)*)', re.MULTILINE)

    with open(filename, 'r') as f:
        data = f.read()

    nv, nf, vertices, faces = data_pattern.search(data).groups()
    nv = int(nv)
    nf = int(nf)

    vertices = np.asarray(list(map(float, vertices.split())), dtype=float).reshape((nv, -1))[:, :3]
    faces = np.asarray(list(map(int, faces.split())), dtype=int).reshape((nf, 4))-1

    return Mesh(vertices, faces, name)


def load_VTU(filename, name=None):
    """Loads VTK file format in the new XML format (vtu file extension for unstructured meshes).

    It relies on the reader from the VTK library.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh

    Note
    ----
    VTU files have a 0-indexing
    """

    _check_file(filename)

    vtk = import_optional_dependency("vtk")

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(filename))
    reader.Update()
    vtk_mesh = reader.GetOutput()

    vertices, faces = _dump_vtk(vtk_mesh)
    return Mesh(vertices, faces, name)


def load_VTP(filename, name=None):
    """Loads VTK file format in the new XML format (vtp file extension for polydata meshes).

    It relies on the reader from the VTK library.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh

    Note
    ----
    VTP files have a 0-indexing
    """
    _check_file(filename)

    vtk = import_optional_dependency("vtk")

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(filename))
    reader.Update()
    vtk_mesh = reader.GetOutput()

    vertices, faces = _dump_vtk(vtk_mesh)
    return Mesh(vertices, faces, name)


def load_VTK(filename, name=None):
    """Loads VTK file format in the legacy format (vtk file extension).

    It relies on the reader from the VTK library.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh

    Note
    ----
    VTU files have a 0-indexing
    """
    _check_file(filename)

    vtk = import_optional_dependency("vtk")

    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(str(filename))
    reader.Update()
    vtk_mesh = reader.GetOutput()

    vertices, faces = _dump_vtk(vtk_mesh)
    return Mesh(vertices, faces, name)


def _dump_vtk(vtk_mesh):
    """Internal driver function that uses the VTK library to read VTK polydata or vtk unstructured grid data structures

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    """

    nv = vtk_mesh.GetNumberOfPoints()
    vertices = np.zeros((nv, 3), dtype=float)
    for k in range(nv):
        vertices[k] = np.array(vtk_mesh.GetPoint(k))

    nf = vtk_mesh.GetNumberOfCells()
    faces = np.zeros((nf, 4), dtype=int)
    for k in range(nf):
        cell = vtk_mesh.GetCell(k)
        nv_facet = cell.GetNumberOfPoints()
        for l in range(nv_facet):
            faces[k][l] = cell.GetPointId(l)
        if nv_facet == 3:
            faces[k][3] = faces[k][0]

    return vertices, faces


def load_STL(filename, name=None):
    """Loads STL file format.

    It relies on the reader from the VTK library. As STL file format maintains a redundant set of vertices for each
    faces of the mesh, it returns a merged list of nodes and connectivity array by using the merge_duplicates function.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh

    Note
    ----
    STL files have a 0-indexing
    """
    vtk = import_optional_dependency("vtk")

    from capytaine.meshes.quality import merge_duplicate_rows

    _check_file(filename)

    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(filename))
    reader.Update()

    data = reader.GetOutputDataObject(0)

    nv = data.GetNumberOfPoints()
    vertices = np.zeros((nv, 3), dtype=float)
    for k in range(nv):
        vertices[k] = np.array(data.GetPoint(k))
    nf = data.GetNumberOfCells()
    faces = np.zeros((nf, 4), dtype=int)
    for k in range(nf):
        cell = data.GetCell(k)
        if cell is not None:
            for l in range(3):
                faces[k][l] = cell.GetPointId(l)
                faces[k][3] = faces[k][0]  # always repeating the first node as stl is triangle only

    # Merging duplicates nodes
    vertices, new_id = merge_duplicate_rows(vertices)
    faces = new_id[faces]

    return Mesh(vertices, faces, name)


def load_NAT(filename, name=None):
    """This function loads natural file format for meshes.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh

    Notes
    -----
    The file format is as follow::

        xsym    ysym
        n    m
        x1    y1    z1
        .
        .
        .
        xn    yn    zn
        i1    j1    k1    l1
        .
        .
        .
        im    jm    km    lm

    where :
    n : number of nodes
    m : number of cells
    x1 y1 z1 : cartesian coordinates of node 1
    i1 j1 k1 l1 : counterclock wise Ids of nodes for cell 1
    if cell 1 is a triangle, i1==l1

    Note
    ----
    NAT files have a 1-indexing
    """

    _check_file(filename)

    ifile = open(filename, 'r')
    ifile.readline()
    nv, nf = list(map(int, ifile.readline().split()))

    vertices = []
    for i in range(nv):
        vertices.append(list(map(float, ifile.readline().split())))
    vertices = np.array(vertices, dtype=float)

    faces = []
    for i in range(nf):
        faces.append(list(map(int, ifile.readline().split())))
    faces = np.array(faces, dtype=int)

    ifile.close()
    return Mesh(vertices, faces-1, name)


def load_GDF(filename, name=None):
    """Loads WAMIT (Wamit INC. (c)) GDF mesh files.

    As GDF file format maintains a redundant set of vertices for each faces of the mesh, it returns a merged list of
    nodes and connectivity array by using the merge_duplicates function.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh or ReflectionSymmetricMesh
        the loaded mesh

    Note
    ----
    GDF files have a 1-indexing
    """

    _check_file(filename)

    with open(str(filename)) as gdf_file:
        title = gdf_file.readline()
        ulen, grav = map(float, gdf_file.readline().split()[:2])
        isx, isy = map(int, gdf_file.readline().split()[:2])
        npan = int(gdf_file.readline().split()[0])
        faces_vertices = np.genfromtxt(gdf_file)

    faces_vertices = faces_vertices.reshape(-1, 3)
    vertices, indices = np.unique(faces_vertices, axis=0, return_inverse=True)
    faces = indices.reshape(-1, 4)

    if faces.shape[0] != npan:
        raise ValueError(
            f"In {filename} npan value: {npan} is not equal to face count: \
                {faces.shape[0]}."
        )

    if isx == 1 and isy == 1:
        return ReflectionSymmetricMesh(ReflectionSymmetricMesh(Mesh(vertices, faces, f"quarter_of_{name}"), yOz_Plane, f"half_of_{name}"), xOz_Plane, name)
    elif isx == 1:
        return ReflectionSymmetricMesh(Mesh(vertices, faces, f"half_of_{name}"), yOz_Plane, name)
    elif isy == 1:
        return ReflectionSymmetricMesh(Mesh(vertices, faces, f"half_of_{name}"), xOz_Plane, name)
    else:
        return Mesh(vertices, faces, name)


def load_MAR(filename, name=None):
    """Loads Nemoh (Ecole Centrale de Nantes) mesh files.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh or ReflectionSymmetry
        the loaded mesh

    Note
    ----
    MAR files have a 1-indexing
    """

    _check_file(filename)

    ifile = open(filename, 'r')

    header = ifile.readline()
    _, symmetric_mesh = header.split()

    vertices = []
    while 1:
        line = ifile.readline()
        line = line.split()
        if line[0] == '0':
            break
        vertices.append(list(map(float, line[1:])))

    vertices = np.array(vertices, dtype=float)
    faces = []
    while 1:
        line = ifile.readline()
        line = line.split()
        if line[0] == '0':
            break
        faces.append(list(map(int, line)))

    faces = np.array(faces, dtype=int)

    ifile.close()

    if int(symmetric_mesh) == 1:
        if name is None:
            half_mesh = Mesh(vertices, faces-1)
            return ReflectionSymmetricMesh(half_mesh, plane=xOz_Plane)
        else:
            half_mesh = Mesh(vertices, faces-1, name=f"half_of_{name}")
            return ReflectionSymmetricMesh(half_mesh, plane=xOz_Plane, name=name)
    else:
        return Mesh(vertices, faces-1, name)


def load_MSH(filename, name=None):
    """Loads .MSH mesh files generated by GMSH by C. Geuzaine and J.F. Remacle.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh

    Note
    ----
    MSH files have a 1-indexing
    """

    import re

    _check_file(filename)

    try:
        meshio = import_optional_dependency("meshio")
    except:
        with open(filename, 'r') as file:
            data = file.read()
        version = float(re.search(r'\$MeshFormat\n(\d.\d).*\n\$EndMeshFormat', data, re.DOTALL).groups()[0])

        if 4 <= version < 5:
            message = (
                f"Meshio is required to read MSH file format version 4. "
                f"Use pip or conda to install Meshio."
            )
            raise ImportError(message) from None
        else:
            nb_nodes, nodes_data = re.search(r'\$Nodes\n(\d+)\n(.+)\$EndNodes', data, re.DOTALL).groups()
            nb_elts, elts_data = re.search(r'\$Elements\n(\d+)\n(.+)\$EndElements', data, re.DOTALL).groups()

            vertices = np.asarray(list(map(float, nodes_data.split())), dtype=float).reshape((-1, 4))[:, 1:]
            vertices = np.ascontiguousarray(vertices)
            faces = []

            # Triangles
            for tri_elt in re.findall(r'(^\d+\s2(?:\s\d+)+?$)', elts_data, re.MULTILINE):
                tri_elt = list(map(int, tri_elt.split()))
                triangle = tri_elt[-3:]
                triangle.append(triangle[0])
                faces.append(triangle)

            for quad_elt in re.findall(r'(^\d+\s3(?:\s\d+)+?$)', elts_data, re.MULTILINE):
                quad_elt = list(map(int, quad_elt.split()))
                quadrangle = quad_elt[-4:]
                faces.append(quadrangle)

            faces = np.asarray(faces, dtype=int) - 1

            return Mesh(vertices, faces, name)

    msh_mesh = meshio.read(filename)
    from capytaine.io.meshio import load_from_meshio
    return load_from_meshio(msh_mesh, name)


def load_MED(filename, name=None):
    """Loads MED mesh files generated by SALOME MECA.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh

    Note
    ----
    MED files have a 1-indexing
    """

    try:
        import h5py
    except ImportError:
        raise ImportError('MED file format reader needs h5py module to be installed')

    _check_file(filename)

    file = h5py.File(filename)

    list_of_names = []
    file.visit(list_of_names.append)

    nb_quadrangles = nb_triangles = 0

    for item in list_of_names:
        if '/NOE/COO' in item:
            vertices = file[item][:].reshape((3, -1)).T
            nv = vertices.shape[0]
        if '/MAI/TR3/NOD' in item:
            triangles = file[item][:].reshape((3, -1)).T - 1
            nb_triangles = triangles.shape[0]
        if '/MAI/QU4/NOD' in item:
            quadrangles = file[item][:].reshape((4, -1)).T - 1
            nb_quadrangles = quadrangles.shape[0]

    file.close()

    if nb_triangles == 0:
        triangles = np.zeros((0, 4), dtype=int)
    else:
        triangles = np.column_stack((triangles, triangles[:, 0]))
    if nb_quadrangles == 0:
        quadrangles = np.zeros((0, 4), dtype=int)

    faces = np.row_stack([triangles, quadrangles])

    return Mesh(vertices, faces, name=name)


def load_WRL(filename, name=None):
    """Loads VRML 2.0 mesh files.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh
    """
    import re

    vtk = import_optional_dependency("vtk")

    _check_file(filename)

    # Checking version
    with open(filename, 'r') as f:
        line = f.readline()
        ver = re.search(r'#VRML\s+V(\d.\d)', line).group(1)
        if not ver == '2.0':
            raise NotImplementedError('VRML loader only supports VRML 2.0 format (version %s given)' % ver)

    importer = vtk.vtkVRMLImporter()
    importer.SetFileName(str(filename))
    importer.Update()

    actors = importer.GetRenderer().GetActors()
    actors.InitTraversal()
    dataset = actors.GetNextActor().GetMapper().GetInput()

    return _dump_vtk(dataset)


def load_NEM(filename, name=None):
    """Loads mesh files that are used by the ``Mesh`` tool included in Nemoh.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh
        the loaded mesh

    Note
    ----
    This format is different from that is used directly by Nemoh software. It is only dedicated to the Mesh tool.
    """

    _check_file(filename)

    ifile = open(filename, 'r')

    nv = int(ifile.readline())
    nf = int(ifile.readline())

    vertices = []
    for ivertex in range(nv):
        vertices.append(list(map(float, ifile.readline().split())))
    vertices = np.asarray(vertices, dtype=float)

    faces = []
    for iface in range(nf):
        faces.append(list(map(int, ifile.readline().split())))
    faces = np.asarray(faces, dtype=int)
    faces -= 1

    return Mesh(vertices, faces, name)


def load_PNL(filename, name=None):
    """Load mesh using HAMS file format.

    Parameters
    ----------
    filename: str
        name of the mesh file on disk

    Returns
    -------
    Mesh or ReflectionSymmetricMesh
        the loaded mesh
    """

    with open(filename, 'r') as f:
        # Skip 3 title lines
        f.readline()
        f.readline()
        f.readline()
        # Read header data
        nb_faces, nb_vertices, x_sym, y_sym = map(int, f.readline().split())
        # Skip 2 more lines
        f.readline()
        f.readline()
        vertices = np.genfromtxt((f.readline() for _ in range(nb_vertices)), usecols=(1, 2, 3))
        # Skip 3 more lines
        f.readline()
        f.readline()
        f.readline()
        faces = np.zeros((nb_faces, 4), dtype=int)
        for i in range(nb_faces):
            index, nb_corners, *data = map(int, f.readline().split())
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
        half_mesh = Mesh(vertices, faces, name=(f"half_of_{name}" if name is not None else None))
        return ReflectionSymmetricMesh(half_mesh, plane=yOz_Plane, name=name)
    elif x_sym == 0 and y_sym == 1:
        half_mesh = Mesh(vertices, faces, name=(f"half_of_{name}" if name is not None else None))
        return ReflectionSymmetricMesh(half_mesh, plane=xOz_Plane, name=name)
    elif x_sym == 1 and y_sym == 1:
        quarter_mesh = Mesh(vertices, faces, name=(f"quarter_of_{name}" if name is not None else None))
        half_mesh = ReflectionSymmetricMesh(quarter_mesh, plane=xOz_Plane, name=(f"half_of_{name}" if name is not None else None))
        return ReflectionSymmetricMesh(half_mesh, plane=yOz_Plane, name=name)
    else:
        return Mesh(vertices, faces, name)


extension_dict = {  # keyword, reader
    'dat': load_MAR,
    'mar': load_MAR,
    'nemoh': load_MAR,
    'wamit': load_GDF,
    'gdf': load_GDF,
    'diodore-inp': load_INP,
    'inp': load_INP,
    'diodore-dat': load_DAT,
    'hydrostar': load_HST,
    'hst': load_HST,
    'natural': load_NAT,
    'nat': load_NAT,
    'gmsh': load_MSH,
    'msh': load_MSH,
    'rad': load_RAD,
    'radioss': load_RAD,
    'stl': load_STL,
    'vtu': load_VTU,
    'vtp': load_VTP,
    'paraview-legacy': load_VTK,
    'vtk': load_VTK,
    'tecplot': load_TEC,
    'tec': load_TEC,
    'med': load_MED,
    'salome': load_MED,
    'vrml': load_WRL,
    'wrl': load_WRL,
    'nem': load_NEM,
    'nemoh_mesh': load_NEM,
    'pnl': load_PNL,
    'hams': load_PNL,
}
