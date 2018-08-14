#!/usr/bin/env python
# coding: utf-8
"""Helper functions to study the connectivities of the mesh.

Based on Meshmagick by Francois Rongere (EC Nantes).
"""

def connectivity(mesh):
    """Compute the connectivities of the mesh.

    It concerns further connectivity than simple faces/vertices connectivities. It computes the vertices / vertices, vertices / faces and faces / faces connectivities.

    Note
    ----
    Note that if the mesh is not conformal, the algorithm may not perform correctly
    """

    nv = mesh.nb_vertices
    nf = mesh.nb_faces

    mesh_closed = True

    # Building connectivities

    # Establishing v_v and v_f connectivities
    v_v = dict([(i, set()) for i in range(nv)])
    v_f = dict([(i, set()) for i in range(nv)])
    for (iface, face) in enumerate(mesh._faces):
        if face[0] == face[-1]:
            face_w = face[:3]
        else:
            face_w = face
        for (index, iV) in enumerate(face_w):
            v_f[iV].add(iface)
            v_v[face_w[index - 1]].add(iV)
            v_v[iV].add(face_w[index - 1])

    # Connectivity f_f
    boundary_edges = dict()

    f_f = dict([(i, set()) for i in range(nf)])
    for ivertex in range(nv):
        set1 = v_f[ivertex]
        for iadj_v in v_v[ivertex]:
            set2 = v_f[iadj_v]
            intersection = list(set1 & set2)
            if len(intersection) == 2:
                f_f[intersection[0]].add(intersection[1])
                f_f[intersection[1]].add(intersection[0])

            elif len(intersection) == 1:
                boundary_face = mesh._faces[intersection[0]]

                if boundary_face[0] == boundary_face[-1]:
                    boundary_face = boundary_face[:3]
                ids = np.where((boundary_face == ivertex) + (boundary_face == iadj_v))[0]

                if ids[1] != ids[0]+1:
                    i_v_orig, i_v_target = boundary_face[ids]
                else:
                    i_v_target, i_v_orig = boundary_face[ids]

                boundary_edges[i_v_orig] = i_v_target
            else:
                raise RuntimeError('Unexpected error while computing mesh connectivities')

    # Computing boundaries
    boundaries = list()
    # TODO: calculer des boundaries fermees et ouvertes (closed_boundaries et open_boundaries) et mettre dans dict
    while True:
        try:
            boundary = list()
            i_v0_init, i_v1 = boundary_edges.popitem()
            boundary.append(i_v0_init)
            boundary.append(i_v1)
            i_v0 = i_v1

            while True:
                try:
                    i_v1 = boundary_edges.pop(i_v0)
                    boundary.append(i_v1)
                    i_v0 = i_v1
                except KeyError:
                    if boundary[0] != boundary[-1]:
                        print('Boundary is not closed !!!')
                    else:
                        boundaries.append(boundary)
                    break
        except KeyError:
            break

    return {'v_v': v_v,
            'v_f': v_f,
            'f_f': f_f,
            'boundaries': boundaries}
