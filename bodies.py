#!/usr/bin/env python
# coding: utf-8

import numpy as np
from itertools import product

import sys
sys.path.append("/home/ancellin/meshmagick/")

from meshmagick.mesh import Mesh

class FloattingBody(Mesh):
    pass


class Sphere(FloattingBody):
    def __init__(self, radius=1.0, ntheta=11, nphi=11, z0=0.0):

        # if z0 < -radius: # fully immerged
        theta_max = np.pi
        # elif z0 < radius:
        #     theta_max = np.arccos(z0/radius)
        # else:
        #     raise Exception("Sphere out of the water")

        theta = np.linspace(-theta_max, theta_max, ntheta)
        phi = np.linspace(-np.pi/2, np.pi/2, nphi)

        # Nodes
        nodes  = np.zeros((ntheta*nphi, 3), dtype=np.float32)

        for i, (t, p) in enumerate(product(theta, phi)):
            # The sign of theta below is a trick to get the correct orientation of the normal vectors...
            x =      radius * np.sin(t) * np.sin(np.sign(t)*p)
            y =      radius * np.sin(t) * np.cos(np.sign(t)*p)
            z = z0 - radius * np.cos(t)
            nodes[i, :] = (x, y, z)

        # Connectivities
        panels = np.zeros(((ntheta-1)*(nphi-1), 4), dtype=np.int)

        for k, (i, j) in enumerate(product(range(0, ntheta-1), range(0, nphi-1))):
            panels[k, :] = (j+i*nphi, j+(i+1)*nphi, j+1+(i+1)*nphi, j+1+i*nphi)

        Mesh.__init__(self, nodes, panels)
        self.merge_duplicates()
        self.heal_triangles()


class HorizontalCylinder(FloattingBody):
    def __init__(self, length=1.0, radius=1.0, nx=11, nr=3, ntheta=11, z0=0.0):

        # if z0 < -radius: # fully immerged
        theta_max = np.pi
        # elif z0 < radius:
        #     theta_max = np.arccos(z0/radius)
        # else:
        #     raise Exception("Cylinder out of the water")

        X = np.linspace(0.0, length, nx)
        R = np.linspace(0.0, radius, nr)
        theta = np.linspace(-theta_max, theta_max, ntheta)

        # Nodes
        nnodes = ntheta*(nx+2*nr)
        nodes  = np.zeros((nnodes, 3), dtype=np.float32)

        for i, (t, x) in enumerate(product(theta, X)):
            y = radius * np.sin(t)
            z = z0 - radius * np.cos(t)
            nodes[i, :] = (x, y, z)

        for i, (x, r, t) in enumerate(product([0, length], R, theta)):
            y = r * np.sin(t)
            z = z0 - r * np.cos(t)
            nodes[ntheta*nx+i, :] = (x, y, z)

        # Connectivities
        npanels = (ntheta-1)*((nx-1)+2*(nr-1))
        panels = np.zeros((npanels, 4), dtype=np.int)

        for k, (i, j) in enumerate(product(range(0, ntheta-1), range(0, nx-1))):
            panels[k, :] = (j+i*nx, j+(i+1)*nx, j+1+(i+1)*nx, j+1+i*nx)

        for k, (i, j) in enumerate(product(range(0, nr-1), range(ntheta*nx, ntheta*nx+ntheta-1))):
            panels[(ntheta-1)*(nx-1)+k, :] = (j+i*ntheta,  j+1+i*ntheta, j+1+(i+1)*ntheta,j+(i+1)*ntheta)

        for k, (i, j) in enumerate(product(range(0, nr-1), range(ntheta*(nx+nr), ntheta*(nx+nr)+ntheta-1))):
            panels[(ntheta-1)*((nx-1)+(nr-1))+k, :] = (j+i*ntheta, j+(i+1)*ntheta, j+1+(i+1)*ntheta, j+1+i*ntheta)

        Mesh.__init__(self, nodes, panels)
        self.merge_duplicates()
        self.heal_triangles()


class OneSidedRectangle(FloattingBody):
    def __init__(self, height=2.0, length=10.0, nh=5, nl=5, z0=0.0):

        X = np.linspace(-length/2, length/2, nl)
        Z = np.linspace(z0, z0+height, nh)

        nodes   = np.zeros((nl*nh, 3), dtype=np.float32)
        panels  = np.zeros(((nl-1)*(nh-1), 4), dtype=np.int)

        for i, (x, y, z) in enumerate(product(X, [0.0], Z)):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, nl-1), range(0, nh-1))):
            panels[k, :] = (j+i*nh, j+1+i*nh, j+1+(i+1)*nh, j+(i+1)*nh)

        Mesh.__init__(self, nodes, panels)


class TwoSidedRectangle(FloattingBody):
    def __init__(self, height=2.0, length=10.0, nh=5, nl=5, z0=0.0):

        X = np.linspace(-length/2, length/2, nl)
        Z = np.linspace(z0, z0+height, nh)

        nodes   = np.zeros((nl*nh, 3), dtype=np.float32)
        panels  = np.zeros((2*(nl-1)*(nh-1), 4), dtype=np.int)

        for i, (x, y, z) in enumerate(product(X, [0.0], Z)):
            nodes[i, :] = x, y, z

        for k, (i, j) in enumerate(product(range(0, nl-1), range(0, nh-1))):
            panels[k, :] = (j+i*nh, j+1+i*nh, j+1+(i+1)*nh, j+(i+1)*nh)
            panels[(nl-1)*(nh-1)+k, :] = (j+i*nh, j+(i+1)*nh, j+1+(i+1)*nh,  j+1+i*nh, )

        Mesh.__init__(self, nodes, panels)


# class OpenRectangularParallelepiped(FloattingBody):
#     def __init__(self, height=2.0, length=10.0, thickness=0.5, nh=5, nl=5, nth=3, z0=0.0):
#     def generate_nodes_and_panels(self, height=2, length=10, thickness=0.5, z0=-1.0, nh=3, nl=11, nth=7):
#         if z0 < 0:
#             height = min(height, -z0)
#         else:
#             raise Exception("Parallelepiped out of the water")

#         X = np.linspace(-length/2, length/2, nl)
#         Z = np.linspace(z0, z0+height, nh)

#         if self.y_sym:
#             if not nth % 2 == 1:
#                 raise Exception("Need an odd number of points to use y-symmetry")
#             nth = nth//2+1
#             Y = np.linspace(0, thickness/2, nth)
#             # nnodes  = 2*(nth*nl + nl*nh + nh*nth)
#             # npanels = 2*((nth-1)*(nl-1) + (nl-1)*(nh-1) + (nh-1)*(nth-1))
#         else:
#             Y = np.linspace(-thickness/2, thickness/2, nth)
#             # nnodes  = nl*nh + 2*(nth*nl + nh*nth)
#             # npanels = 2*((nth-1)*(nl-1) + (nl-1)*(nh-1) + (nh-1)*(nth-1))

#         # Nodes and panels
#         nodes   = np.zeros((nnodes, 3), dtype=np.float32)
#         panels  = np.zeros((npanels, 5), dtype=np.int32)

#         # Top panel
#         i0 = 0
#         nodes[:, :] = product(X, Y, [z0+height])
#         panels[:, :] = 
#         for i in range(nl-1):
#             for j in range(1, nth):
#                 connectivites.append(f"{i0+j+i*nth}\t{i0+j+(i+1)*nth}\t{i0+j+1+(i+1)*nth}\t{i0+j+1+i*nth}\n")

#         # Side panel
#         points.extend(product([-length/2], Y, Z))
#         for i in range(nth-1):
#             for j in range(1, nh):
#                 connectivites.append(f"{i0+j+i*nh}\t{i0+j+1+i*nh}\t{i0+j+1+(i+1)*nh}\t{i0+j+(i+1)*nh}\n")

#         # Side panel
#         i0 = len(points)
#         points.extend(product([length/2], Y, Z))
#         for i in range(nth-1):
#             for j in range(1, nh):
#                 connectivites.append(f"{i0+j+i*nh}\t{i0+j+(i+1)*nh}\t{i0+j+1+(i+1)*nh}\t{i0+j+1+i*nh}\n")

#         # Front panel
#         i0 = len(points)
#         points.extend(product(X, [thickness/2], Z))
#         for i in range(nl-1):
#             for j in range(1, nh):
#                 connectivites.append(f"{i0+j+i*nh}\t{i0+j+1+i*nh}\t{i0+j+1+(i+1)*nh}\t{i0+j+(i+1)*nh}\n")

#         if not self.y_sym:
#             # Back panel
#             i0 = len(points)
#             points.extend(product(X, [-thickness/2], Z))
#             for i in range(nl-1):
#                 for j in range(1, nh):
#                     connectivites.append(f"{i0+j+i*nh}\t{i0+j+(i+1)*nh}\t{i0+j+1+(i+1)*nh}\t{i0+j+1+i*nh}\n")


if __name__ == "__main__":
    # hc = HorizontalCylinder()
    # hc.show_matplotlib()

    # sp = Sphere(radius=1.0, ntheta=11, nphi=11)
    # sp.show_matplotlib()

    bd = TwoSidedRectangle()
    bd.translate_y(3.0)
    bd = bd + TwoSidedRectangle()
    bd.show_matplotlib()
