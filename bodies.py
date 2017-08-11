#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy.linalg import norm
from itertools import product

class FloattingBody:
    def __init__(self, y_sym=False, **kwargs):
        self.y_sym = y_sym
        self.generate_nodes_and_panels(**kwargs)
        self.panels[:, 4] = self.panels[:, 0]

        self.compute_normal_and_area()
        self.dof = {}

    def compute_normal_and_area(self):
        self.normal         = np.zeros((self.npanels, 3), dtype=np.float32)
        self.center_of_mass = np.zeros((self.npanels, 3), dtype=np.float32)
        self.area           = np.zeros((self.npanels),    dtype=np.float32)
        self.radius         = np.zeros((self.npanels),    dtype=np.float32)

        for i, panel in enumerate(self.panels):
            # Area of 1-2-4 triangle
            w0 = np.cross(
                    self.nodes[panel[1], :] - self.nodes[panel[0], :],
                    self.nodes[panel[3], :] - self.nodes[panel[1], :],
                    )

            # Area of 1-2-3 triangle
            w1 = np.cross(
                    self.nodes[panel[3], :] - self.nodes[panel[2], :],
                    self.nodes[panel[1], :] - self.nodes[panel[2], :],
                    )

            self.area[i] = (norm(w0) + norm(w1))/2
            self.normal[i, :] = (w0+w1)/norm(w0+w1)
            self.center_of_mass[i, :] = (
                    (self.nodes[panel[0], :] + self.nodes[panel[1], :] + self.nodes[panel[3], :])*norm(w0)/(6*self.area[i]) +
                    (self.nodes[panel[1], :] + self.nodes[panel[2], :] + self.nodes[panel[3], :])*norm(w1)/(6*self.area[i])
                    )

            self.radius[i] = max([norm(node - self.center_of_mass[i, :]) for node in self.nodes[panel[0:4], :]])

    def plot_mesh(self, normal_vectors=False):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")

        faces = []
        for face in self.panels:
            vertices = []
            for index_vertex in face:
                vertices.append(self.nodes[int(index_vertex), :])
            faces.append(vertices)
        ax.add_collection3d(Poly3DCollection(faces, facecolor=(0, 0, 0, 0), edgecolor='k'))

        if normal_vectors:
            ax.quiver(*zip(*self.center_of_mass), *zip(*self.normal), length=0.2)

        plt.xlabel("x")
        plt.ylabel("y")
        plt.xlim(min(self.nodes[:, 0]), max(self.nodes[:, 0]))
        plt.ylim(min(self.nodes[:, 1]), max(self.nodes[:, 1]))
        plt.gca().set_zlim(min(self.nodes[:, 2]), max(self.nodes[:, 2]))
        plt.show()

    def as_Nemoh_file(self):
        if self.y_sym:
            s = "\t2\t1\n"
        else:
            s = "\t2\t0\n"

        for i, node in enumerate(self.nodes):
            s += f"{i+1}\t{node[0]:f}\t{node[1]:f}\t{node[2]:f}\n"

        s += "0\t0\t0\t0\n"

        for panel in self.panels:
            s += f"{panel[0]+1}\t{panel[1]+1}\t{panel[2]+1}\t{panel[3]+1}\n"

        s += "0\t0\t0\t0\n"
        return s
    
    def save_as_file(self, filename):
        with open(filename, 'w') as f:
            f.write(self.as_Nemoh_file())


class Sphere(FloattingBody):
    def generate_nodes_and_panels(self, radius=1.0, z0=-3.0, ntheta=11, nphi=11):
        if z0 < -radius: # fully immerged
            theta_max = np.pi
        elif z0 < radius:
            theta_max = np.arccos(z0/radius)
        else:
            raise Exception("Sphere out of the water")

        if self.y_sym:
            if not ntheta % 2 == 1:
                raise Exception("Need an odd number of points to use y-symmetry")
            ntheta = ntheta//2+1
            theta = np.linspace(0, theta_max, ntheta)
        else:
            theta = np.linspace(-theta_max, theta_max, ntheta)
        phi = np.linspace(-np.pi/2, np.pi/2, nphi)

        # Nodes
        self.nnodes = ntheta*nphi
        self.nodes  = np.zeros((self.nnodes, 3), dtype=np.float32)

        for i, (t, p) in enumerate(product(theta, phi)):
            # The sign of theta is a trick to get the correct orientation of the normal vectors...
            x =      radius * np.sin(t) * np.sin(np.sign(t)*p)
            y =      radius * np.sin(t) * np.cos(np.sign(t)*p)
            z = z0 - radius * np.cos(t)
            self.nodes[i, :] = (x, y, z)

        # Connectivities
        self.npanels = (ntheta-1)*(nphi-1)
        self.panels = np.zeros((self.npanels, 5), dtype=np.int32)

        for k, (i, j) in enumerate(product(range(0, ntheta-1), range(0, nphi-1))):
            self.panels[k, 0:4] = (j+i*nphi, j+(i+1)*nphi, j+1+(i+1)*nphi, j+1+i*nphi)


class HorizontalCylinder(FloattingBody):
    def generate_nodes_and_panels(self, length=1.0, radius=1.0, z0=-2.0, nx=11, nr=3, ntheta=11):

        if z0 < -radius: # fully immerged
            theta_max = np.pi
        elif z0 < radius:
            theta_max = np.arccos(z0/radius)
        else:
            raise Exception("Cylinder out of the water")

        if self.y_sym:
            if not ntheta % 2 == 1:
                raise Exception("Need an odd number of points to use y-symmetry")
            ntheta = ntheta//2+1
            theta = np.linspace(0, theta_max, ntheta)
        else:
            theta = np.linspace(-theta_max, theta_max, ntheta)

        X = np.linspace(0.0, length, nx)
        R = np.linspace(0.0, radius, nr)

        # Nodes
        self.nnodes = ntheta*(nx+2*nr)
        self.nodes  = np.zeros((self.nnodes, 3), dtype=np.float32)

        for i, (t, x) in enumerate(product(theta, X)):
            y = radius * np.sin(t)
            z = z0 - radius * np.cos(t)
            self.nodes[i, :] = (x, y, z)

        for i, (x, r, t) in enumerate(product([0, length], R, theta)):
            y = r * np.sin(t)
            z = z0 - r * np.cos(t)
            self.nodes[ntheta*nx+i, :] = (x, y, z)

        # Connectivities
        self.npanels = (ntheta-1)*((nx-1)+2*(nr-1))
        self.panels = np.zeros((self.npanels, 5), dtype=np.int32)

        for k, (i, j) in enumerate(product(range(0, ntheta-1), range(0, nx-1))):
            self.panels[k, 0:4] = (j+i*nx, j+(i+1)*nx, j+1+(i+1)*nx, j+1+i*nx)

        for k, (i, j) in enumerate(product(range(0, nr-1), range(ntheta*nx, ntheta*nx+ntheta-1))):
            self.panels[(ntheta-1)*(nx-1)+k, 0:4] = (j+i*ntheta,  j+1+i*ntheta, j+1+(i+1)*ntheta,j+(i+1)*ntheta)

        for k, (i, j) in enumerate(product(range(0, nr-1), range(ntheta*(nx+nr), ntheta*(nx+nr)+ntheta-1))):
            self.panels[(ntheta-1)*((nx-1)+(nr-1))+k, 0:4] = (j+i*ntheta, j+(i+1)*ntheta, j+1+(i+1)*ntheta, j+1+i*ntheta)


# class Rectangle(FloattingBody):
#     def generate_nodes_and_panels(self, height=2, length=10, y0=0.0, z0=-1.0, nh=3, nl=11):
#         if z0 < 0:
#             height = min(height, -z0)
#         else:
#             raise Exception("Rectangle out of the water")

#         X = np.linspace(-length/2, length/2, nl)
#         Z = np.linspace(z0, z0+height, nh)

#         self.nodes   = np.zeros((self.nnodes, 3), dtype=np.float32)
#         self.panels  = np.zeros((self.npanels, 5), dtype=np.int32)

#         (product(X, [thickness/2], Z))

# class RectangularParallelepiped(FloattingBody):
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
#             # self.nnodes  = 2*(nth*nl + nl*nh + nh*nth)
#             # self.npanels = 2*((nth-1)*(nl-1) + (nl-1)*(nh-1) + (nh-1)*(nth-1))
#         else:
#             Y = np.linspace(-thickness/2, thickness/2, nth)
#             # self.nnodes  = nl*nh + 2*(nth*nl + nh*nth)
#             # self.npanels = 2*((nth-1)*(nl-1) + (nl-1)*(nh-1) + (nh-1)*(nth-1))

#         # Nodes and panels
#         self.nodes   = np.zeros((self.nnodes, 3), dtype=np.float32)
#         self.panels  = np.zeros((self.npanels, 5), dtype=np.int32)

#         # Top panel
#         i0 = 0
#         self.nodes[:, :] = product(X, Y, [z0+height])
#         self.panels[:, :] = 
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
    hc = HorizontalCylinder(y_sym=True)
    hc.plot_mesh(normal_vectors=True)
    # hc.save_as_file("test.dat")

    # sp = Sphere(radius=1.0, z0=-3.0, ntheta=11, nphi=11)
    # sp.plot_mesh(normal_vectors=True)
