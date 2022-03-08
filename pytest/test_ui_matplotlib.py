#!/usr/bin/env python
# coding: utf-8
"""Test of matplotlib user interface methods."""

import pytest
import numpy as np
import capytaine as cpt

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

@pytest.mark.skipif(plt is None,
                     reason='matplotlib is not installed')
def test_showmatplotlib_with_colors(tmp_path):
    cylinder = cpt.HorizontalCylinder(
        length=5.0, radius=1.0, center=(0, 0, -2), nr=5, nx=20, ntheta=20,
    )

    problem = cpt.DiffractionProblem(body=cylinder, wave_direction=0.0, omega=1.0)
    solver = cpt.BEMSolver()
    results = solver.solve(problem)

    cylinder.show_matplotlib(saveas=tmp_path / "showmpl.png")

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    import matplotlib.cm as cm
    colormap = cm.get_cmap('cividis')
    cylinder.show_matplotlib(color_field=np.abs(results.potential),
                             ax=ax, cbar_label=r'$\phi (m^2/s)$')
    ax.set_title('Potential distribution')
    plt.savefig(tmp_path / "showmpl_with_color_field.png")
