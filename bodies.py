#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy.linalg import norm
from itertools import product

class FloattingBody:
    def __init__(self, **kwargs):
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


class HorizontalCylinder(FloattingBody):
    def generate_nodes_and_panels(self, length=1.0, radius=1.0, z0=-2.0, nx=10, nr=3, ntheta=10):
        if z0 < -radius: # fully immerged
            theta_max = np.pi
        elif z0 < radius:
            theta_max = np.arccos(z0/radius)
        else:
            raise Exception("Cylinder out of the water")

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


class Sphere(FloattingBody):
    def generate_nodes_and_panels(self, radius=1.0, z0=-3.0, ntheta=11, nphi=11):
        if z0 < -radius: # fully immerged
            theta_max = np.pi
        elif z0 < radius:
            theta_max = np.arccos(z0/radius)
        else:
            raise Exception("Sphere out of the water")

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


if __name__ == "__main__":
    # hc = HorizontalCylinder()
    # hc.plot_mesh()

    sp = Sphere(radius=1.0, z0=-3.0, ntheta=11, nphi=11)
    sp.plot_mesh()
