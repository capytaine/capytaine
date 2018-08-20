#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import os
import time
import numpy as np

real_str = r'[+-]?(?:\d+\.\d*|\d*\.\d+)(?:[Ee][+-]?\d+)?'  # Regex for floats
# =======================================================================
# MESH LOADERS
# ======================================================================
# Contains here all functions to load meshes from different file formats


def _check_file(filename):
    if not os.path.isfile(filename):
        raise IOError("file %s not found" % filename)
    return


def load_mesh(filename, file_format=None):
    """Driver function that loads every mesh file format known by meshmagick and returns the node list and the
    connectivity array

    Parameters
    ----------
    filename: str
        name of the meh file on disk
    file_format: str, optional
        format of the mesh defined in the extension_dict dictionary

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    """
    _check_file(filename)

    if file_format is None:
        _, file_format = os.path.splitext(filename)

    if file_format not in extension_dict:
        raise IOError('Extension ".%s" is not known' % file_format)

    loader = extension_dict[file_format][0]

    vertices, faces = loader(filename)

    return vertices, faces


def load_RAD(filename):
    """Loads RADIOSS mesh files. This export file format may be chosen in ICEM meshing program.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

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
    node_line = r'\s*\d+\s*(' + real_str + ')\s*(' + real_str + ')\s*(' + real_str + ')'
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
    faces = np.asarray(faces, dtype=np.int) - 1

    return vertices, faces


def load_HST(filename):
    """Loads HYDROSTAR (Bureau Veritas (c)) mesh files.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

    Note
    ----
    HST files have a 1-indexing
    """
    _check_file(filename)

    ifile = open(filename, 'r')
    data = ifile.read()
    ifile.close()

    import re

    node_line = r'\s*\d+(?:\s+' + real_str + '){3}'
    node_section = r'((?:' + node_line + ')+)'

    elem_line = r'^\s*(?:\d+\s+){3}\d+\s*[\r\n]+'
    elem_section = r'((?:' + elem_line + ')+)'

    pattern_node_line = re.compile(node_line, re.MULTILINE)
    pattern_elem_line = re.compile(elem_line, re.MULTILINE)
    pattern_node_section = re.compile(node_section, re.MULTILINE)
    pattern_elem_section = re.compile(elem_section, re.MULTILINE)

    vertices_tmp = []
    vertices = []
    nv = 0
    for node_section in pattern_node_section.findall(data):
        for node in pattern_node_line.findall(node_section):
            vertices_tmp.append(list(map(float, node.split()[1:])))
        nv_tmp = len(vertices_tmp)
        vertices_tmp = np.asarray(vertices_tmp, dtype=np.float)
        if nv == 0:
            vertices = vertices_tmp.copy()
            nv = nv_tmp
        else:
            vertices = np.concatenate((vertices, vertices_tmp))
            nv += nv_tmp

    faces_tmp = []
    faces = []
    nf = 0
    for elem_section in pattern_elem_section.findall(data):
        for elem in pattern_elem_line.findall(elem_section):
            faces_tmp.append(list(map(int, elem.split())))
        nf_tmp = len(faces_tmp)
        faces_tmp = np.asarray(faces_tmp, dtype=np.int)
        if nf == 0:
            faces = faces_tmp.copy()
            nf = nf_tmp
        else:
            faces = np.concatenate((faces, faces_tmp))
            nf += nf_tmp

    return vertices, faces-1


def load_DAT(filename):
    """Not implemented.
    Intended to load .DAT files used in DIODORE (PRINCIPIA (c))
    """
    _check_file(filename)
    raise NotImplementedError


def load_INP(filename):
    """Loads DIODORE (PRINCIPIA (c)) configuration file format.

    It parses the .INP file and extract meshes defined in subsequent .DAT files using the different informations
    contained in the .INP file.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

    Note
    ----
    INP/DAT files use a 1-indexing
    """
    _check_file(filename)
    import re

    with open(filename, 'r') as f:
        text = f.read()

    # Retrieving frames into a dictionnary frames
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
        id_new = - np.ones(max(idx_array) + 1, dtype=np.int)
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
            nodes = np.asarray(mesh_files[file]['NODE_SECTION'], dtype=np.float)
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
                                      dtype=np.int)

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

    return vertices, faces-1


def load_TEC(filename):
    """Loads TECPLOT (Tecplot (c)) mesh files.

    It relies on the tecplot file reader from the VTK library.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

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

    vertices = np.asarray(list(map(float, vertices.split())), dtype=np.float).reshape((nv, -1))[:, :3]
    faces = np.asarray(list(map(int, faces.split())), dtype=np.int).reshape((nf, 4))-1

    return vertices, faces


def load_VTU(filename):
    """Loads VTK file format in the new XML format (vtu file extension for unstructured meshes).

    It relies on the reader from the VTK library.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

    Note
    ----
    VTU files have a 0-indexing
    """

    _check_file(filename)

    from vtk import vtkXMLUnstructuredGridReader
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    vtk_mesh = reader.GetOutput()

    vertices, faces = _dump_vtk(vtk_mesh)
    return vertices, faces


def load_VTP(filename):
    """Loads VTK file format in the new XML format (vtp file extension for polydata meshes).

    It relies on the reader from the VTK library.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

    Note
    ----
    VTP files have a 0-indexing
    """
    _check_file(filename)

    from vtk import vtkXMLPolyDataReader
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    vtk_mesh = reader.GetOutput()

    vertices, faces = _dump_vtk(vtk_mesh)
    return vertices, faces


def load_VTK(filename):
    """Loads VTK file format in the legacy format (vtk file extension).

    It relies on the reader from the VTK library.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

    Note
    ----
    VTU files have a 0-indexing
    """
    _check_file(filename)

    from vtk import vtkPolyDataReader
    reader = vtkPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    vtk_mesh = reader.GetOutput()

    vertices, faces = _dump_vtk(vtk_mesh)
    return vertices, faces


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
    vertices = np.zeros((nv, 3), dtype=np.float)
    for k in range(nv):
        vertices[k] = np.array(vtk_mesh.GetPoint(k))

    nf = vtk_mesh.GetNumberOfCells()
    faces = np.zeros((nf, 4), dtype=np.int)
    for k in range(nf):
        cell = vtk_mesh.GetCell(k)
        nv_facet = cell.GetNumberOfPoints()
        for l in range(nv_facet):
            faces[k][l] = cell.GetPointId(l)
        if nv_facet == 3:
            faces[k][3] = faces[k][0]

    return vertices, faces


def load_STL(filename):
    """Loads STL file format.

    It relies on the reader from the VTK library. As STL file format maintains a redundant set of vertices for each
    faces of the mesh, it returns a merged list of nodes and connectivity array by using the merge_duplicates function.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

    Note
    ----
    STL files have a 0-indexing
    """
    from vtk import vtkSTLReader
    from .tools import merge_duplicate_rows

    _check_file(filename)

    reader = vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()

    data = reader.GetOutputDataObject(0)

    nv = data.GetNumberOfPoints()
    vertices = np.zeros((nv, 3), dtype=np.float)
    for k in range(nv):
        vertices[k] = np.array(data.GetPoint(k))
    nf = data.GetNumberOfCells()
    faces = np.zeros((nf, 4), dtype=np.int)
    for k in range(nf):
        cell = data.GetCell(k)
        if cell is not None:
            for l in range(3):
                faces[k][l] = cell.GetPointId(l)
                faces[k][3] = faces[k][0]  # always repeating the first node as stl is triangle only

    # Merging duplicates nodes
    vertices, new_id = merge_duplicate_rows(vertices, return_index=True)
    faces = new_id[faces]

    return vertices, faces


def load_NAT(filename):
    """This function loads natural file format for meshes.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        array of the coordinates of the mesh's nodes
    faces: ndarray
        array of the faces nodes connectivity

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
    vertices = np.array(vertices, dtype=np.float)

    faces = []
    for i in range(nf):
        faces.append(list(map(int, ifile.readline().split())))
    faces = np.array(faces, dtype=np.int)

    ifile.close()
    return vertices, faces-1


def load_GDF(filename):
    """Loads WAMIT (Wamit INC. (c)) GDF mesh files.

    As GDF file format maintains a redundant set of vertices for each faces of the mesh, it returns a merged list of
    nodes and connectivity array by using the merge_duplicates function.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

    Note
    ----
    GDF files have a 1-indexing
    """

    _check_file(filename)

    ifile = open(filename, 'r')

    ifile.readline()  # skip one header line
    line = ifile.readline().split()
    ulen = line[0]
    grav = line[1]

    line = ifile.readline().split()
    isx = line[0]
    isy = line[1]

    line = ifile.readline().split()
    nf = int(line[0])

    vertices = np.zeros((4 * nf, 3), dtype=np.float)
    faces = np.zeros((nf, 4), dtype=np.int)

    iv = -1
    for icell in range(nf):

        for k in range(4):
            iv += 1
            vertices[iv, :] = np.array(ifile.readline().split())
            faces[icell, k] = iv

    ifile.close()

    return vertices, faces


def load_MAR(filename):
    """Loads Nemoh (Ecole Centrale de Nantes) mesh files.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

    Note
    ----
    MAR files have a 1-indexing
    """

    _check_file(filename)

    ifile = open(filename, 'r')

    ifile.readline()  # Skipping the first line of the file
    vertices = []
    while 1:
        line = ifile.readline()
        line = line.split()
        if line[0] == '0':
            break
        vertices.append(list(map(float, line[1:])))

    vertices = np.array(vertices, dtype=np.float)
    faces = []
    while 1:
        line = ifile.readline()
        line = line.split()
        if line[0] == '0':
            break
        faces.append(list(map(int, line)))

    faces = np.array(faces, dtype=np.int)

    ifile.close()

    return vertices, faces-1


def load_MSH(filename):
    """Loads .MSH mesh files generated by GMSH by C. Geuzaine and J.F. Remacle.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

    Note
    ----
    MSH files have a 1-indexing
    """

    import re

    _check_file(filename)

    with open(filename, 'r') as file:
        data = file.read()

    nb_nodes, nodes_data = re.search(r'\$Nodes\n(\d+)\n(.+)\$EndNodes', data, re.DOTALL).groups()
    nb_elts, elts_data = re.search(r'\$Elements\n(\d+)\n(.+)\$EndElements', data, re.DOTALL).groups()

    vertices = np.asarray(list(map(float, nodes_data.split())), dtype=np.float).reshape((-1, 4))[:, 1:]
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

    faces = np.asarray(faces, dtype=np.int) - 1

    return vertices, faces


def load_MED(filename):
    """Loads MED mesh files generated by SALOME MECA.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

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

    # TODO: gerer les cas ou on a que des tris ou que des quads...
    nb_quadrangles = nb_triangles = 0

    for item in list_of_names:
        if '/NOE/COO' in item:
            vertices = file.get(item).value.reshape((3, -1)).T
            nv = vertices.shape[0]
        if '/MAI/TR3/NOD' in item:
            triangles = file.get(item).value.reshape((3, -1)).T - 1
            nb_triangles = triangles.shape[0]
        if '/MAI/QU4/NOD' in item:
            quadrangles = file.get(item).value.reshape((4, -1)).T - 1
            nb_quadrangles = quadrangles.shape[0]

    file.close()

    if nb_triangles == 0:
        triangles = np.zeros((0, 4), dtype=np.int)
    else:
        triangles = np.column_stack((triangles, triangles[:, 0]))
    if nb_quadrangles == 0:
        quadrangles = np.zeros((0, 4), dtype=np.int)

    faces = np.zeros((nb_triangles+nb_quadrangles, 4), dtype=np.int)
    faces[:nb_triangles] = triangles
    # faces[:nb_triangles, -1] = triangles[:, 0]
    faces[nb_triangles:] = quadrangles

    vertices = np.ascontiguousarray(vertices)
    return vertices, faces


def load_WRL(filename):
    """Loads VRML 2.0 mesh files.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    """

    from vtk import vtkVRMLImporter
    import re

    _check_file(filename)

    # Checking version
    with open(filename, 'r') as f:
        line = f.readline()
        ver = re.search(r'#VRML\s+V(\d.\d)', line).group(1)
        if not ver == '2.0':
            raise NotImplementedError('VRML loader only supports VRML 2.0 format (version %s given)' % ver)

    importer = vtkVRMLImporter()
    importer.SetFileName(filename)
    importer.Update()

    actors = importer.GetRenderer().GetActors()
    actors.InitTraversal()
    dataset = actors.GetNextActor().GetMapper().GetInput()

    return _dump_vtk(dataset)


def load_NEM(filename):
    """Loads mesh files that are used by the ``Mesh`` tool included in Nemoh.

    Parameters
    ----------
    filename: str
        name of the meh file on disk

    Returns
    -------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

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
    vertices = np.asarray(vertices, dtype=np.float)

    faces = []
    for iface in range(nf):
        faces.append(list(map(int, ifile.readline().split())))
    faces = np.asarray(faces, dtype=np.int)
    faces -= 1

    return vertices, faces


#=======================================================================
#                             MESH WRITERS
#=======================================================================
# Contains here all functions to write meshes in different file formats

def write_mesh(filename, vertices, faces, file_format):
    """Driver function that writes every mesh file file_format known by meshmagick

    Parameters
    ----------
    filename: str
        name of the mesh file to be written on disk
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    file_format: str
        file_format of the mesh defined in the extension_dict dictionary

    """

    if file_format not in extension_dict:
        raise IOError('Extension "%s" is not known' % file_format)

    writer = extension_dict[file_format][1]

    writer(filename, vertices, faces)


def write_DAT(filename, vertices, faces):
    """Writes .DAT file format for the DIODORE (PRINCIPA (c)) software.

    It also displays suggestions for inclusion into the .INP configuration
    file.

    Parameters
    ----------
    filename: str
        name of the mesh file to be written on disk
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    """

    import os

    root_filename, ext = os.path.splitext(filename)
    filename = root_filename + ext.upper()
    ofile = open(filename, 'w')

    ofile.write('$\n$ Data for DIODORE input file : {0}\n'.format(root_filename.upper()))
    ofile.write('$ GENERATED BY MESHMAGICK ON {0}\n$\n'.format(time.strftime('%c')))

    ofile.write('$ NODE\n')
    vertex_block = \
        ''.join(
            (
                '\n'.join(
                    ''.join(
                        (
                            '{:8d}'.format(idx+1),
                            ''.join('{:13.5E}'.format(elt) for elt in node)
                        )
                    ) for (idx, node) in enumerate(vertices)
                ),

                '\n*RETURN\n'
            )
        )
    ofile.write(vertex_block)

    quad_block = '$\n$ ELEMENT,TYPE=Q4C000,ELSTRUCTURE={0}'.format(root_filename.upper())
    tri_block = '$\n$ ELEMENT,TYPE=T3C000,ELSTRUCTURE={0}'.format(root_filename.upper())
    nq = 0
    nt = 0
    for (idx, cell) in enumerate(faces+1):
        if cell[0] != cell[-1]:
            # quadrangle
            nq += 1
            quad_block = ''.join(
                (
                    quad_block,
                    '\n',
                    '{:8d}'.format(idx+1),
                    ''.join('{:8d}'.format(node_id) for node_id in cell)
                )
            )

        else:
            # Triangle
            nt += 1
            tri_block = ''.join(
                (
                    tri_block,
                    '\n',
                    '{:8d}'.format(idx+1),
                    ''.join('{:8d}'.format(node_id) for node_id in cell[:3])
                )
            )

    print('-------------------------------------------------')
    print('Suggestion for .inp DIODORE input file :')
    print('')
    print('*NODE,INPUT={0},FRAME=???'.format(root_filename))

    if nq > 0:
        quad_block = ''.join((quad_block, '\n*RETURN\n'))
        ofile.write(quad_block)
        print('*ELEMENT,TYPE=Q4C000,ELSTRUCTURE={0},INPUT={0}'.format(root_filename))
    if nt > 0:
        tri_block = ''.join((tri_block, '\n*RETURN\n'))
        ofile.write(tri_block)
        print('*ELEMENT,TYPE=T3C000,ELSTRUCTURE={0},INPUT={0}'.format(root_filename))

    print('')
    print('-------------------------------------------------')
    ofile.close()


def write_HST(filename, vertices, faces):
    """Writes .HST file format for the HYDROSTAR (Bureau Veritas (c)) software.

    Parameters
    ----------
    filename: str
        name of the mesh file to be written on disk
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    """
    # TODO: allow many bodies

    ofile = open(filename, 'w')

    ofile.write(''.join((
        'PROJECT:\n',
        'USERS:   meshmagick\n\n'
        'NBODY   1\n'
        'RHO   1025.0\n'
        'GRAVITY   9.81\n\n'
    )))

    coordinates_block = ''.join((  # block
            'COORDINATES\n',
            '\n'.join(  # line
                ''.join(
                    (
                        '{:10d}'.format(idx+1),  # index
                        ''.join('{:16.6E}'.format(elt) for elt in node)  # node coordinates
                    )
                ) for (idx, node) in enumerate(vertices)
            ),
            '\nENDCOORDINATES\n\n'
    ))

    ofile.write(coordinates_block)

    cells_coordinates = ''.join((  # block
        'PANEL TYPE 0\n',
        '\n'.join(  # line
            ''.join(
                '{:10d}'.format(node_idx) for node_idx in cell
            ) for cell in faces + 1
        ),
        '\nENDPANEL\n\n'
    ))

    ofile.write(cells_coordinates)

    ofile.write('ENDFILE\n')

    ofile.close()


def write_TEC(filename, vertices, faces):
    """Writes .TEC file format for the TECPLOT (Tecplot (c)) visualisation software.

    It relies on the VTK library for its writer.

    Parameters
    ----------
    filename: str
        name of the mesh file to be written on disk
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    """

    ofile = open(filename, 'w')

    nv = vertices.shape[0]
    nf = faces.shape[0]

    ofile.write('TITLE = \" THIS FILE WAS GENERATED BY MESHMAGICK - FICHIER : {} \" \n'.format(filename))

    ofile.write('VARIABLES = \"X\",\"Y\",\"Z\" \n')
    ofile.write('ZONE T=\"MESH\" \n')
    ofile.write('N={nv:10d} ,E={nf:10d} , F=FEPOINT, ET=QUADRILATERAL\n'.format(nv=nv, nf=nf))

    node_block = '\n'.join( # block
        ''.join(
            ''.join('{:16.6E}'.format(elt) for elt in node)
        ) for node in vertices
    ) + '\n'
    ofile.write(node_block)

    cells_block = '\n'.join(  # block
        ''.join(
            ''.join('{:10d}'.format(node_id) for node_id in cell)
        ) for cell in faces + 1
    ) + '\n'
    ofile.write(cells_block)

    ofile.close()

    return 1


def write_VTU(filename, vertices, faces):
    """Writes .vtu file format for the paraview (Kitware (c)) visualisation software.

    It relies on the VTK library for its writer. VTU files use the last XML file format of the VTK library.

    Parameters
    ----------
    filename: str
        name of the mesh file to be written on disk
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    """

    from vtk import vtkXMLUnstructuredGridWriter, VTK_MAJOR_VERSION
    writer = vtkXMLUnstructuredGridWriter()
    writer.SetDataModeToAscii()
    writer.SetFileName(filename)

    unstructured_grid = _build_vtkUnstructuredGrid(vertices, faces)
    if VTK_MAJOR_VERSION <= 5:
        writer.SetInput(unstructured_grid)
    else:
        writer.SetInputData(unstructured_grid)
    writer.Write()


def write_VTP(filename, vertices, faces):
    """Writes .vtp file format for the Paraview (Kitware (c)) visualisation software.

    It relies on the VTK library for its writer. VTP files use the last XML file format of the VTK library and
    correspond to polydata.

    Parameters
    ----------
    filename: str
        name of the mesh file to be written on disk
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    """

    from vtk import vtkXMLPolyDataWriter, VTK_MAJOR_VERSION
    writer = vtkXMLPolyDataWriter()
    writer.SetDataModeToAscii()
    writer.SetFileName(filename)

    polydata = _build_vtkPolyData(vertices, faces)
    if VTK_MAJOR_VERSION <= 5:
        writer.SetInput(polydata)
    else:
        writer.SetInputData(polydata)
    writer.Write()


def write_VTK(filename, vertices, faces):
    """Writes .vtk file format for the Paraview (Kitware (c)) visualisation software.

    It relies on the VTK library for its writer. VTK files use the legagy ASCII file format of the VTK library.

    Parameters
    ----------
    filename: str
        name of the mesh file to be written on disk
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    """

    nv = vertices.shape[0]
    nf = faces.shape[0]

    triangle_mask = (faces[:, 0] == faces[:, -1])
    quadrangles_mask = np.invert(triangle_mask)
    nb_triangles = len(np.where(triangle_mask)[0])
    nb_quandrangles = len(np.where(quadrangles_mask)[0])

    with open(filename, 'w') as f:

        f.write('# vtk DataFile Version 4.0\n')
        f.write('vtk file generated by meshmagick on %s\n' % time.strftime('%c'))
        f.write('ASCII\n')
        f.write('DATASET POLYDATA\n')
        f.write('POINTS %u float\n' % nv)

        for vertex in vertices:
            f.write('%f %f %f\n' % (vertex[0], vertex[1], vertex[2]))

        f.write('POLYGONS %u %u\n' % (nf, 4*nb_triangles+5*nb_quandrangles))

        for face in faces:
            if face[0] == face[-1]:  # Triangle
                f.write('3 %u %u %u\n' % (face[0], face[1], face[2]))
            else:  # Quadrangle
                f.write('4 %u %u %u %u\n' % (face[0], face[1], face[2], face[3]))


def _build_vtkUnstructuredGrid(vertices, faces):
    """Internal function that builds a VTK object for manipulation by the VTK library.

    Parameters
    ----------
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

    Returns
    -------
    vtkObject
    """

    import vtk

    nv = max(np.shape(vertices))
    nf = max(np.shape(faces))

    vtk_mesh = vtk.vtkUnstructuredGrid()
    vtk_mesh.Allocate(nf, nf)

    # Building the vtkPoints data structure
    vtk_points = vtk.vtkPoints()
    vtk_points.SetNumberOfPoints(nv)
    for idx, vertex in enumerate(vertices):
        vtk_points.SetPoint(idx, vertex)

    vtk_mesh.SetPoints(vtk_points)  # Storing the points into vtk_mesh

    # Building the vtkCell data structure
    for cell in faces:
        if cell[-1] in cell[:-1]:
            vtk_cell = vtk.vtkTriangle()
            nc = 3
        else:
            # #print 'quadrangle'
            vtk_cell = vtk.vtkQuad()
            nc = 4

        for k in range(nc):
            vtk_cell.GetPointIds().SetId(k, cell[k])

        vtk_mesh.InsertNextCell(vtk_cell.GetCellType(), vtk_cell.GetPointIds())
    return vtk_mesh


def _build_vtkPolyData(vertices, faces):
    """Builds a vtkPolyData object from vertices and faces"""

    import vtk

    # Create a vtkPoints object and store the points in it
    points = vtk.vtkPoints()
    for point in vertices:
        points.InsertNextPoint(point)

    # Create a vtkCellArray to store faces
    cell_array = vtk.vtkCellArray()
    for face_ids in faces:
        if face_ids[0] == face_ids[-1]:
            # Triangle
            curface = face_ids[:3]
            vtk_face = vtk.vtkTriangle()
        else:
            # Quadrangle
            curface = face_ids[:4]
            vtk_face = vtk.vtkQuad()

        for idx, id in enumerate(curface):
            vtk_face.GetPointIds().SetId(idx, id)

        cell_array.InsertNextCell(vtk_face)

    polydata_mesh = vtk.vtkPolyData()
    polydata_mesh.SetPoints(points)
    polydata_mesh.SetPolys(cell_array)

    return polydata_mesh


def write_NAT(filename, vertices, faces):
    """Writes .nat file format as defined into the load_NAT function.

    Parameters
    ----------
    filename: str
        name of the mesh file to be written on disk
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

    See Also
    --------
    load_NAT
    """

    ofile = open(filename, 'w')

    nv = max(np.shape(vertices))
    nf = max(np.shape(faces))

    ofile.write('%6u%6u\n' % (0, 0))  # lire les symmetries dans args...
    ofile.write('%6u%6u\n' % (nv, nf))
    for vertex in vertices:
        ofile.write('%15.6E%15.6E%15.6E\n' % (vertex[0], vertex[1], vertex[2]))
    for cell in faces+1:
        ofile.write('%10u%10u%10u%10u\n' % (cell[0], cell[1], cell[2], cell[3]))

    ofile.close()


def write_NEM(filename, vertices, faces):
    """Writes mesh files used by the Mesh tool included in Nemoh

    Parameters
    ----------
    filename : str
        name of the mesh file to be written on disk
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities

    Note
    ----
    This file format is different from that used by Nemoh itself. It is only used by the Mesh tool.
    """
    ofile = open(filename, 'w')

    ofile.write('%u\n' % vertices.shape[0])
    ofile.write('%u\n' % faces.shape[0])

    for vertex in vertices:
        ofile.write('%15.6f\t%15.6f\t%15.6f\n' % (vertex[0], vertex[1], vertex[2]))

    for face in faces+1:
        ofile.write('%10u\t%10u\t%10u\t%10u\n' % (face[0], face[1], face[2], face[3]))

    ofile.close()


def write_GDF(filename, vertices, faces):
    """Writes .gdf file format for the WAMIT (Wamit INC. (c)) BEM software.

    Parameters
    ----------
    filename: str
        name of the mesh file to be written on disk
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    """

    nf = max(np.shape(faces))

    ofile = open(filename, 'w')

    ofile.write('GDF file generated by meshmagick on %s\n' % time.strftime('%c'))

    ofile.write('%16.6f%16.6f\n' % (100.0, 9.81))
    ofile.write('%12u%12u\n' % (0, 1))  # TODO : mettre les symetries en argument
    ofile.write('%12u\n' % nf)

    for cell in faces:
        for k in range(4):
            cur_vertices = vertices[cell[k], :]
            ofile.write('%16.6E%16.6E%16.6E\n' % (cur_vertices[0], cur_vertices[1], cur_vertices[2]))

    ofile.close()


def write_MAR(filename, vertices, faces):
    """Writes mesh files to be used with Nemoh BEM software (Ecole Centrale de Nantes)

    Parameters
    ----------
    filename: str
        name of the mesh file to be written on disk
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    """

    # TODO: detect symmetry in Oxz plane

    ofile = open(filename, 'w')

    ofile.write('{0:6d}{1:6d}\n'.format(2, 0))  # TODO : mettre les symetries en argument

    for (idx, vertex) in enumerate(vertices):
        ofile.write('{0:6d}{1:16.6f}{2:16.6f}{3:16.6f}\n'.format(idx+1, vertex[0], vertex[1], vertex[2]))

    ofile.write('{0:6d}{1:6d}{2:6d}{3:6d}{4:6d}\n'.format(0, 0, 0, 0, 0))

    cell_block = '\n'.join(
        ''.join('{0:10d}'.format(elt) for elt in cell)
        for cell in faces + 1
    ) + '\n'
    ofile.write(cell_block)
    ofile.write('%6u%6u%6u%6u\n' % (0, 0, 0, 0))

    ofile.close()

    print('WARNING: if you described only one part of the mesh using symmetry for Nemoh, you may manually modify the ' \
          'file header accordingly')


def write_RAD(filename, vertices, faces):
    raise NotImplementedError


def write_STL(filename, vertices, faces):
    """Writes .stl file format. It relies on the VTK library for its writer.

    Parameters
    ----------
    filename: str
        name of the mesh file to be written on disk
    vertices: ndarray
        numpy array of the coordinates of the mesh's nodes
    faces: ndarray
        numpy array of the faces' nodes connectivities
    """

    # TODO : replace this implementation by using the vtk functionalities

    # Triangulating quads
    t1 = (0, 1, 2)
    t2 = (0, 2, 3)

    quads_ids = np.where(faces[:, 0] != faces[:, -1])[0]

    new_faces = faces[quads_ids].copy()
    new_faces[:, :3] = new_faces[:, t1]
    new_faces[:, -1] = new_faces[:, 0]

    faces[quads_ids, :3] = faces[:, t2][quads_ids]
    faces[quads_ids, -1] = faces[quads_ids, 0]

    faces = np.concatenate((faces, new_faces))

    # Writing file
    ofile = open(filename, 'w')

    ofile.write('solid meshmagick\n')

    for face in faces:
        if face[0] != face[3]:
            raise RuntimeError("""Only full triangle meshes are accepted in STL files.
              Please consider using the --triangulate-quadrangles option (-tq) to
              perform a prior triangulation of the mesh""")

        # Computing normal
        v0 = vertices[face[0], :]
        v1 = vertices[face[1], :]
        v2 = vertices[face[2], :]

        n = np.cross(v1 - v0, v2 - v0)
        n /= np.linalg.norm(n)

        block_facet = ''.join(['  facet normal ', ''.join('%15.6e' % ni for ni in n) + '\n',
                               '    outer loop\n',
                               '      vertex', ''.join('%15.6e' % Vi for Vi in v0) + '\n',
                               '      vertex', ''.join('%15.6e' % Vi for Vi in v1) + '\n',
                               '      vertex', ''.join('%15.6e' % Vi for Vi in v2) + '\n',
                               '    endloop\n',
                               '  endfacet\n'])
        ofile.write(block_facet)
    ofile.write('endsolid meshmagick\n')
    ofile.close()


def write_INP(filename, vertices, faces):
    raise NotImplementedError('INP writer is not implementer yet')


def write_MSH(filename, vertices, faces):
    raise NotImplementedError('MSH writer is not implemented yet')


def write_MED(filename, vertices, faces):
    raise NotImplementedError('MED writer is not implemented yet')


def write_WRL(filename, vertices, faces):
    raise NotImplementedError('VRML writer is not implemented yet')


def know_extension(ext):
    return ext in extension_dict

extension_dict = {  # keyword,  reader,   writer
    'mar': (load_MAR, write_MAR),
    'nemoh': (load_MAR, write_MAR),
    'wamit': (load_GDF, write_GDF),
    'gdf': (load_GDF, write_GDF),
    'diodore-inp': (load_INP, write_INP),
    'inp': (load_INP, write_INP),
    'diodore-dat': (load_DAT, write_DAT),
    'hydrostar': (load_HST, write_HST),
    'hst': (load_HST, write_HST),
    'natural': (load_NAT, write_NAT),
    'nat': (load_NAT, write_NAT),
    'gmsh': (load_MSH, write_MSH),
    'msh': (load_MSH, write_MSH),
    'rad': (load_RAD, write_RAD),
    'radioss': (load_RAD, write_RAD),
    'stl': (load_STL, write_STL),
    'vtu': (load_VTU, write_VTU),
    'vtp': (load_VTP, write_VTP),
    'paraview-legacy': (load_VTK, write_VTK),
    'vtk': (load_VTK, write_VTK),
    'tecplot': (load_TEC, write_TEC),
    'tec': (load_TEC, write_TEC),
    'med': (load_MED, write_MED),
    'salome': (load_MED, write_MED),
    'vrml': (load_WRL, write_WRL),
    'wrl': (load_WRL, write_WRL),
    'nem': (load_NEM, write_NEM),
    'nemoh_mesh': (load_NEM, write_NEM)
}
