#!/usr/bin/env python
# coding: utf-8
"""Matplotlib animation for the free surface elevation."""
# This file is part of "capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import animation

from capytaine.results import LinearPotentialFlowResult
from capytaine.geometric_bodies.free_surface import FreeSurface


def animation_matplotlib(result: LinearPotentialFlowResult, fs: FreeSurface,
                         complex_face_elevation=None, body_shape=None, fps=24):

    X = fs.mesh.faces_centers[:, 0].reshape(fs.nx, fs.ny)
    Y = fs.mesh.faces_centers[:, 1].reshape(fs.nx, fs.ny)
    fs_elevation = complex_face_elevation.reshape(fs.nx, fs.ny)
    scale = np.abs(fs_elevation).max()

    fig = plt.figure()
    ax = plt.gca()

    frames_per_period = int(result.period * fps)

    def animate(i_frame):
        ax.clear()

        # Draw body shape
        if body_shape is not None:
            ax.add_artist(body_shape)

        # Draw free surface elevation
        ax.contourf(X, Y, np.abs(fs_elevation)*np.cos(np.angle(fs_elevation) - 2*np.pi*i_frame/frames_per_period), 30,
                    vmin=-scale, vmax=scale)

        return ax

    plt.contourf(X, Y, np.real(fs_elevation), 30, vmin=-scale, vmax=scale)
    plt.colorbar()

    plt.axis('equal')
    ani = animation.FuncAnimation(fig, animate, frames_per_period, interval=1000/fps, blit=False)

    plt.show()
