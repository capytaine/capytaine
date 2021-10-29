#!/usr/bin/env python
# coding: utf-8
"""Tools for mesh quality and mesh healing.
Based on meshmagick <https://github.com/LHEEA/meshmagick> by François Rongère.
"""
# Copyright (C) 2017-2019 Matthieu Ancellin, based on the work of François Rongère
# See LICENSE file at <https://github.com/mancellin/capytaine>

import logging

import numpy as np

from capytaine.meshes.geometry import inplace_transformation

LOG = logging.getLogger(__name__)


def merge_duplicates(mesh, atol=1e-8):
    """Merges the duplicate vertices of the mesh in place.

    Parameters
    ----------
    atol : float, optional
        Absolute tolerance. default is 1e-8

    Returns
    -------
    new_id : ndarray
        Array of indices that merges the vertices.
    """
    uniq, new_id = merge_duplicate_rows(mesh.vertices, atol=atol)

    nv_init = mesh.nb_vertices

    # Updating mesh data
    mesh.vertices = uniq
    mesh.faces = new_id[mesh.faces] # Faces vertices ids are updated here

    nv_final = mesh.nb_vertices

    LOG.debug("* Merging duplicate vertices that lie in an absolute proximity of %.1E...", atol)
    delta_n = nv_init - nv_final
    if delta_n == 0:
        LOG.debug("\t--> No duplicate vertices have been found")
    else:
        LOG.debug("\t--> Initial number of vertices : %u", nv_init)
        LOG.debug("\t--> Final number of vertices   : %u", nv_final)
        LOG.debug("\t--> %u vertices have been merged", delta_n)

    # if mesh._has_connectivity():
    #     mesh._remove_connectivity()

    return new_id


def merge_duplicate_rows(arr, atol=1e-8):
    """Returns a new node array where close nodes have been merged into one node (following atol).

    Parameters
    ----------
    arr : array_like
        array of the coordinates of the mesh's nodes
    atol : float, optional
        the tolerance used to define nodes that are coincident and
        that have to be merged

    Returns
    -------
    arr : ndarray
        array of the coordinates of the mesh's nodes where
        every node is different
    newID : ndarray
        array of the new new vertices IDs
    """
    # This function is a bottleneck in the clipping routines
    # TODO: use np.unique to cluster groups --> acceleration !!

    # atol = pow(10, -decimals)

    arr = np.asarray(arr)

    nv, nbdim = arr.shape

    levels = [0, nv]
    iperm = np.arange(nv)

    for dim in range(nbdim):
        # Sorting the first dimension
        values = arr[:, dim].copy()
        if dim > 0:
            values = values[iperm]
        levels_tmp = []
        for (ilevel, istart) in enumerate(levels[:-1]):
            istop = levels[ilevel+1]

            if istop-istart > 1:
                level_values = values[istart:istop]
                iperm_view = iperm[istart:istop]

                iperm_tmp = level_values.argsort()

                level_values[:] = level_values[iperm_tmp]
                iperm_view[:] = iperm_view[iperm_tmp]

                levels_tmp.append(istart)
                vref = values[istart]

                for idx in range(istart, istop):
                    cur_val = values[idx]
                    if np.abs(cur_val - vref) > atol:
                        levels_tmp.append(idx)
                        vref = cur_val

            else:
                levels_tmp.append(levels[ilevel])
        if len(levels_tmp) == nv:
            # No duplicate rows
            # if verbose:
            # LOG.debug "\t -> No duplicate _vertices detected :)"
            newID = np.arange(nv)

        levels_tmp.append(nv)
        levels = levels_tmp

    else:
        # Building the new merged node list
        arr_tmp = []
        newID = np.arange(nv)
        for (ilevel, istart) in enumerate(levels[:-1]):
            istop = levels[ilevel+1]

            arr_tmp.append(arr[iperm[istart]])
            newID[iperm[list(range(istart, istop))]] = ilevel
        arr = np.array(arr_tmp, dtype=float)
        # Applying renumbering to cells
        # if F is not None:
        #     for cell in F:
        #         cell[:] = newID[cell]

        # if verbose:
        # nv_new = arr.shape[0]
        # LOG.debug "\t -> Initial number of nodes : {:d}".format(nv)
        # LOG.debug "\t -> New number of nodes     : {:d}".format(nv_new)
        # LOG.debug "\t -> {:d} nodes have been merged".format(nv-nv_new)

    # if F is not None:
    #     if return_index:
    #         return arr, F, newID
    #     else:
    #         return arr, F
    # else:
    return arr, newID


@inplace_transformation
def heal_normals(mesh):
    """Heals the mesh's normals orientations so that they have a consistent orientation and try to make them outward.
    """
    # TODO: return the different groups of a mesh in case it is made of several unrelated groups

    nv = mesh.nb_vertices
    nf = mesh.nb_faces
    faces = mesh._faces

    # Building connectivities
    v_v = mesh.vv
    v_f = mesh.vf
    f_f = mesh.ff
    boundaries = mesh.boundaries

    if len(boundaries) > 0:
        mesh_closed = False
    else:
        mesh_closed = True

    # Flooding the mesh to find inconsistent normals
    type_cell = np.zeros(nf, dtype=np.int32)
    type_cell[:] = 4
    type_cell[mesh.triangles_ids] = 3

    f_vis = np.zeros(nf, dtype=bool)
    f_vis[0] = True
    stack = [0]
    nb_reversed = 0
    while 1:
        if len(stack) == 0:
            if np.any(np.logical_not(f_vis)):
                iface = np.where(np.logical_not(f_vis))[0][0]
                stack.append(iface)
                f_vis[iface] = True
            else:
                break

        iface = stack.pop()
        face = faces[iface]
        s1 = set(face)

        for iadj_f in f_f[iface]:
            if f_vis[iadj_f]:
                continue
            f_vis[iadj_f] = True
            # Removing the other pointer
            f_f[iadj_f].remove(iface)  # So as it won't go from iadj_f to iface in the future

            # Shared vertices
            adjface = faces[iadj_f]
            s2 = set(adjface)
            # try:
            common_vertices = list(s1 & s2)
            if len(common_vertices) == 2:
                i_v1, i_v2 = common_vertices
            else:
                LOG.warning('faces %u and %u have more than 2 vertices in common !', iface, iadj_f)
                continue

            # Checking normal consistency
            face_ref = np.roll(face[:type_cell[iface]], -np.where(face == i_v1)[0][0])
            adj_face_ref = np.roll(adjface[:type_cell[iadj_f]], -np.where(adjface == i_v1)[0][0])

            if face_ref[1] == i_v2:
                i = 1
            else:
                i = -1

            if adj_face_ref[i] == i_v2:
                # Reversing normal
                nb_reversed += 1
                faces[iadj_f] = np.flipud(faces[iadj_f])

            # Appending to the stack
            stack.append(iadj_f)

    LOG.debug("* Healing normals to make them consistent and if possible outward")
    if nb_reversed > 0:
        LOG.debug('\t--> %u faces have been reversed to make normals consistent across the mesh' % (nb_reversed))
    else:
        LOG.debug("\t--> Normals orientations are consistent")

    mesh._faces = faces

    # Checking if the normals are outward
    if mesh_closed:
        zmax = np.max(mesh._vertices[:, 2])

        areas = mesh.faces_areas
        normals = mesh.faces_normals
        centers = mesh.faces_centers
        # areas, normals, centers = get_all_faces_properties(vertices, faces)

        hs = (np.array([(centers[:, 2] - zmax) * areas, ] * 3).T * normals).sum(axis=0)

        tol = 1e-9
        if np.fabs(hs[0]) > tol or np.fabs(hs[1]) > tol:
            LOG.warning("\t--> the mesh does not seem watertight although marked as closed...")

        if hs[2] < 0:
            flipped = True
            mesh.flip_normals()
        else:
            flipped = False

        if flipped:
            LOG.debug('\t--> Every normals have been reversed to be outward')

    else:
        LOG.info("\t--> Mesh is not closed, meshmagick cannot test if the normals are outward")

    return mesh


@inplace_transformation
def remove_unused_vertices(mesh):
    """Removes unused vertices in the mesh in place.

    Those are vertices that are not used by any face connectivity.
    """
    # TODO: implementer return_index !!
    nv = mesh.nb_vertices
    vertices, faces = mesh._vertices, mesh._faces

    used_v = np.zeros(nv, dtype=bool)
    used_v[sum(list(map(list, faces)), [])] = True
    nb_used_v = sum(used_v)

    if nb_used_v < nv:
        new_id__v = np.arange(nv)
        new_id__v[used_v] = np.arange(nb_used_v)
        faces = new_id__v[faces]
        vertices = vertices[used_v]

    mesh._vertices, mesh._faces = vertices, faces

    LOG.debug("* Removing unused vertices in the mesh:")
    if nb_used_v < nv:
        unused_v = np.where(np.logical_not(used_v))[0]
        vlist_str = '[' + ', '.join(str(iV) for iV in unused_v) + ']'
        LOG.debug("\t--> %u unused vertices have been removed" % (nv - nb_used_v))
    else:
        LOG.debug("\t--> No unused vertices")

    return mesh


@inplace_transformation
def heal_triangles(mesh):
    """Makes the triangle connectivity consistent (in place).

    A general face is stored internally as a 4 integer array. It allows to describe indices of a quadrangle's vertices. For triangles, the first index should be equal to the last. This method ensures that this rule is applied everywhere and correct bad triangles description.
    """
    faces = mesh._faces

    quads = faces[:, 0] != faces[:, -1]
    nquads_init = sum(quads)

    faces[quads] = np.roll(faces[quads], 1, axis=1)
    quads = faces[:, 0] != faces[:, -1]

    faces[quads] = np.roll(faces[quads], 1, axis=1)
    quads = faces[:, 0] != faces[:, -1]

    faces[quads] = np.roll(faces[quads], 1, axis=1)
    quads = faces[:, 0] != faces[:, -1]
    nquads_final = sum(quads)

    mesh._faces = faces

    LOG.debug("* Ensuring consistent definition of triangles:")
    if nquads_final < nquads_init:
        LOG.debug("\t--> %u triangles were described the wrong way and have been corrected" % (
        nquads_init - nquads_final))
    else:
        LOG.debug("\t--> Triangle description is consistent")

    return mesh


@inplace_transformation
def remove_degenerated_faces(mesh, rtol=1e-5):
    """Removes tiny triangles from the mesh (in place).

    Tiny triangles are those whose area is lower than the mean triangle area in the mesh times the relative
    tolerance given.

    Parameters
    ----------
    rtol : float, optional
        Positive relative tolerance
    """

    assert 0 < rtol

    # TODO: implementer un retour d'index des faces extraites
    areas = mesh.faces_areas
    area_threshold = areas.mean() * float(rtol)

    # Detecting faces that have null area
    faces = mesh._faces[np.logical_not(areas < area_threshold)]
    nb_removed = mesh.nb_faces - faces.shape[0]
    LOG.debug('* Removing degenerated faces')
    if nb_removed > 0:
        LOG.debug('\t-->%u degenerated faces have been removed' % nb_removed)
    else:
        LOG.debug('\t--> No degenerated faces')

    mesh._faces = faces

    return mesh


def print_quality(mesh):
    """Returns data on the mesh quality.
    Needs to be tested...

    It uses VTK and is reproduced from
    http://vtk.org/gitweb?p=VTK.git;a=blob;f=Filters/Verdict/Testing/Python/MeshQuality.py
    """
    # This function is reproduced from
    # http://vtk.org/gitweb?p=VTK.git;a=blob;f=Filters/Verdict/Testing/Python/MeshQuality.py
    polydata = mesh._vtk_polydata()
    quality = vtk.vtkMeshQuality()
    quality.SetInputData(polydata)

    def DumpQualityStats(iq, arrayname):
        an = iq.GetOutput().GetFieldData().GetArray(arrayname)
        cardinality = an.GetComponent(0, 4)
        range = list()
        range.append(an.GetComponent(0, 0))
        range.append(an.GetComponent(0, 2))
        average = an.GetComponent(0, 1)
        stdDev = math.sqrt(math.fabs(an.GetComponent(0, 3)))
        outStr = '%s%g%s%g\n%s%g%s%g' % (
            '    range: ', range[0], '  -  ', range[1],
            '    average: ', average, '  , standard deviation: ', stdDev)
        return outStr

    # Here we define the various mesh types and labels for output.
    meshTypes = [
        ['Triangle', 'Triangle',
         [['QualityMeasureToArea', ' Area Ratio:'],
          ['QualityMeasureToEdgeRatio', ' Edge Ratio:'],
          ['QualityMeasureToAspectRatio', ' Aspect Ratio:'],
          ['QualityMeasureToRadiusRatio', ' Radius Ratio:'],
          ['QualityMeasureToAspectFrobenius', ' Frobenius Norm:'],
          ['QualityMeasureToMinAngle', ' Minimal Angle:']
          ]
         ],

        ['Quad', 'Quadrilateral',
         [['QualityMeasureToArea', ' Area Ratio:'],
          ['QualityMeasureToEdgeRatio', ' Edge Ratio:'],
          ['QualityMeasureToAspectRatio', ' Aspect Ratio:'],
          ['QualityMeasureToRadiusRatio', ' Radius Ratio:'],
          ['QualityMeasureToMedAspectFrobenius',
           ' Average Frobenius Norm:'],
          ['QualityMeasureToMaxAspectFrobenius',
           ' Maximal Frobenius Norm:'],
          ['QualityMeasureToMinAngle', ' Minimal Angle:']
          ]
         ]
    ]
    res = ''
    if polydata.GetNumberOfCells() > 0:
        for meshType in meshTypes:
            res += '\n%s%s' % (meshType[1], ' quality of the mesh ')
            quality.Update()
            an = quality.GetOutput().GetFieldData().GetArray('Mesh ' + meshType[1] + ' Quality')
            cardinality = an.GetComponent(0, 4)

            res = ''.join((res, '(%u elements):\n' % cardinality))

            # res += '('+str(cardinality) +meshType[1]+'):\n'

            for measure in meshType[2]:
                eval('quality.Set' + meshType[0] + measure[0] + '()')
                quality.Update()
                res += '\n%s\n%s' % (
                    measure[1],
                    DumpQualityStats(quality, 'Mesh ' + meshType[1] + ' Quality')
                )
            res += '\n'

    info = """\n\nDefinition of the different quality measures is given
    in the verdict library manual :
    http://www.vtk.org/Wiki/images/6/6b/VerdictManual-revA.pdf\n"""

    res += info
    print(res)
    return
