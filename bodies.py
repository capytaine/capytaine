#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy.linalg import norm
from itertools import product

class FloattingBody():

    def compute_normal_and_area(self):
        self.normal         = np.zeros((self.npanels, 3), dtype=np.float32)
        self.centre_of_mass = np.zeros((self.npanels, 3), dtype=np.float32)
        self.area           = np.zeros((self.npanels),    dtype=np.float32)

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
            self.normal[i, :] = w0+w1/norm(w0+w1)
            self.centre_of_mass[i, :] = (
                    (self.nodes[panel[0], :] + self.nodes[panel[1], :] + self.nodes[panel[3], :])*norm(w0)/(6*self.area[i]) +
                    (self.nodes[panel[1], :] + self.nodes[panel[2], :] + self.nodes[panel[3], :])*norm(w1)/(6*self.area[i])
                    )

    def plot_mesh(self):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        faces = []
        cdgs = []
        normal_vectors = []
        for face in self.panels:
            vertices = []
            for index_vertex in face:
                vertices.append(self.nodes[int(index_vertex), :])
            faces.append(vertices)
        ax.add_collection3d(Poly3DCollection(faces, facecolor=(0, 0, 0, 0), edgecolor='k'))
        ax.quiver(*zip(*self.centre_of_mass), *zip(*self.normal), length=0.2)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.xlim(min(self.nodes[:, 0]), max(self.nodes[:, 0]))
        plt.ylim(min(self.nodes[:, 1]), max(self.nodes[:, 1]))
        plt.gca().set_zlim(min(self.nodes[:, 2]), max(self.nodes[:, 2]))
        plt.show()

class HorizontalCylinder(FloattingBody):

    def __init__(self, length=1.0, radius=1.0, z0=-2.0, nx=10, nr=3, ntheta=10):

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
        self.panels = np.zeros((self.npanels, 4), dtype=np.int32)

        for k, (i, j) in enumerate(product(range(0, ntheta-1), range(0, nx-1))):
            self.panels[k, :] = (j+i*nx, j+(i+1)*nx, j+1+(i+1)*nx, j+1+i*nx)

        for k, (i, j) in enumerate(product(range(0, nr-1), range(ntheta*nx, ntheta*nx+ntheta-1))):
            self.panels[(ntheta-1)*(nx-1)+k, :] = (j+i*ntheta,  j+1+i*ntheta, j+1+(i+1)*ntheta,j+(i+1)*ntheta)

        for k, (i, j) in enumerate(product(range(0, nr-1), range(ntheta*(nx+nr), ntheta*(nx+nr)+ntheta-1))):
            self.panels[(ntheta-1)*((nx-1)+(nr-1))+k, :] = (j+i*ntheta, j+(i+1)*ntheta, j+1+(i+1)*ntheta, j+1+i*ntheta)

        # Normal vectors and areas
        self.compute_normal_and_area()

        self.dof = {}

if __name__ == "__main__":
    hc = HorizontalCylinder()
    hc.plot_mesh()

