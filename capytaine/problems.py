#!/usr/bin/env python
# coding: utf-8
"""
Definition of the problems to solve with the BEM.
"""

from warnings import warn

import numpy as np

from capytaine._Wavenumber import invert_xtanhx
from capytaine.symmetries import PlanarSymmetry

class RadiationProblem:
    """A radiation problem to be solved by the BEM solver."""

    def __init__(self, body, free_surface=0.0, sea_bottom=-np.infty, omega=1.0, rho=1000.0, g=9.81):
        self.rho = rho
        self.g = g
        self.omega = omega

        if free_surface < sea_bottom:
            raise Exception("Sea bottom is above the free surface.")

        self.free_surface = free_surface
        self.sea_bottom = sea_bottom

        if self.depth == np.infty or omega**2*self.depth/g > 20:
            self.wavenumber = omega**2/g
        else:
            self.wavenumber = invert_xtanhx(omega**2*self.depth/g)/self.depth

        if any(body.vertices[:, 2] > free_surface) or any(body.vertices[:, 2] < sea_bottom):
            warn(f"""The mesh of the body {body.name} is not inside the domain.\nUse body.get_immersed_part() to clip the mesh.""")
        self.body = body

    @property
    def depth(self):
        return self.free_surface - self.sea_bottom

    def export_as_Nemoh_directory(self, directory_name, omega_range=None):
        import os
        from meshmagick.mmio import write_MAR

        if os.path.isdir(directory_name):
            warn(f"Exporting problem in existing directory: {directory_name}")
        else:
            os.makedirs(directory_name)

        write_MAR(
            os.path.join(directory_name, f'{self.body.name}.dat'),
            self.body.vertices,
            self.body.faces,
            xOz_symmetry=isinstance(self.body, PlanarSymmetry)
        )


        if omega_range is None:
            omega_nb_steps = 1
            omega_start = omega
            omega_stop = omega
        else:
            omega_nb_steps = len(omega_range)
            omega_start = min(omega_range)
            omega_stop = max(omega_range)

        with open(os.path.join(directory_name, "Nemoh.cal"), "w") as nemoh_cal:
            nemoh_cal.write(
                    DEFAULT_NEMOH_CAL.format(
                        rho=self.rho,
                        g=self.g,
                        depth=self.depth,
                        mesh_filename=f'{self.body.name}.dat',
                        mesh_vertices=self.body.nb_vertices,
                        mesh_faces=self.body.nb_faces,
                        omega_nb_steps=omega_nb_steps,
                        omega_start=omega_start,
                        omega_stop=omega_stop,
                        )
                    )

        with open(os.path.join(directory_name, "input.txt"), "w") as input_txt:
            input_txt.write(DEFAULT_INPUT_TXT)

        with open(os.path.join(directory_name, "ID.dat"), "w") as id_dat:
            id_dat.write(f"1\n.")


DEFAULT_NEMOH_CAL = """--- Environment ------------------------------------------------------------------------------------------------------------------
{rho}			! RHO 			! KG/M**3 	! Fluid specific volume
{g}				! G			! M/S**2	! Gravity
{depth}			! DEPTH			! M		! Water depth
0.	0.			! XEFF YEFF		! M		! Wave measurement point
--- Description of floating bodies -----------------------------------------------------------------------------------------------
1				! Number of bodies
--- Body 1 -----------------------------------------------------------------------------------------------------------------------
{mesh_filename}
{mesh_vertices} {mesh_faces}
1				! Number of degrees of freedom
1 0. 0. 1. 0. 0. 0.		! Heave
1				! Number of resulting generalised forces
1 0. 0. 1. 0. 0. 0.		! Heave
0				! Number of lines of additional information
--- Load cases to be solved -------------------------------------------------------------------------------------------------------
{omega_nb_steps} {omega_start} {omega_stop} ! Frequencies range
0	0.	0.		! Number of wave directions, Min and Max (degrees)
--- Post processing ---------------------------------------------------------------------------------------------------------------
0	0.1	10.		! IRF 				! IRF calculation (0 for no calculation), time step and duration
0				! Show pressure
0 	0.	180.		! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)
0	0	100.	100.	! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction
"""

DEFAULT_INPUT_TXT="""--- Calculation parameters ------------------------------------------------------------------------------------
1				! Indiq_solver		! - 		! Solver (0) Direct Gauss (1) GMRES (2) GMRES with FMM acceleration (2 not implemented yet)
20				! IRES			! - 		! Restart parameter for GMRES
5.E-07				! TOL_GMRES		! -		! Stopping criterion for GMRES
100				! MAXIT			! - 		! Maximum iterations for GMRES
1				! Sav_potential		! -		! Save potential for visualization
"""
