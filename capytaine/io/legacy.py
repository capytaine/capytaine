#!/usr/bin/env python
# coding: utf-8
"""Import or export Nemoh.cal files for backward compatibility with Nemoh 2."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import os
import logging

import numpy as np

from capytaine.io.mesh_writers import write_MAR
from capytaine.bodies.bodies import FloatingBody
from capytaine.bem.problems_and_results import DiffractionProblem, RadiationProblem
from capytaine.meshes.geometry import Axis

LOG = logging.getLogger(__name__)


def import_cal_file(filepath):
    """Read a Nemoh.cal file and return a list of problems."""

    with open(filepath, 'r') as cal_file:

        cal_file.readline()  # Unused line.
        rho = float(cal_file.readline().split()[0])
        g = float(cal_file.readline().split()[0])
        depth = float(cal_file.readline().split()[0])
        if depth == 0.0:
            sea_bottom = -np.infty
        else:
            sea_bottom = -depth
        xeff, yeff = (float(x) for x in cal_file.readline().split()[0:2])

        bodies = []

        cal_file.readline()  # Unused line.
        nb_bodies = int(cal_file.readline().split()[0])
        for _ in range(nb_bodies):
            cal_file.readline()  # Unused line.
            mesh_file = cal_file.readline().split()[0].strip()
            mesh_file = os.path.join(os.path.dirname(filepath), mesh_file)  # mesh path are relative to Nemoh.cal
            cal_file.readline()  # Number of points, number of panels (unused)

            if os.path.splitext(mesh_file)[1] == '.py':
                from importlib.util import spec_from_file_location, module_from_spec
                spec = spec_from_file_location("body_initialization_module", mesh_file)
                body_initialization = module_from_spec(spec)
                spec.loader.exec_module(body_initialization)
                body = body_initialization.body
            else:
                body = FloatingBody.from_file(mesh_file)

            nb_dofs = int(cal_file.readline().split()[0])
            for i_dof in range(nb_dofs):
                dof_data = cal_file.readline().split()
                if int(dof_data[0]) == 1:
                    direction = np.array([float(x) for x in dof_data[1:4]])
                    body.add_translation_dof(direction=direction)
                elif int(dof_data[0]) == 2:
                    direction = np.array([float(x) for x in dof_data[1:4]])
                    center_of_mass = np.array([float(x) for x in dof_data[4:7]])
                    body.add_rotation_dof(Axis(vector=direction, point=center_of_mass))

            nb_forces = int(cal_file.readline().split()[0])
            for i_force in range(nb_forces):
                force_data = cal_file.readline().split()
                if int(force_data[0]) == 1:
                    direction = np.array([float(x) for x in force_data[1:4]])
                elif int(force_data[0]) == 2:
                    direction = np.array([float(x) for x in force_data[1:4]])
                    center_of_mass = np.array([float(x) for x in force_data[4:7]])
            # TODO: use the generalized forces.

            nb_additional_lines = int(cal_file.readline().split()[0])
            for _ in range(nb_additional_lines):
                cal_file.readline()  # The additional lines are just ignored.

            bodies.append(body)

        if nb_bodies > 1:
            bodies = FloatingBody.join_bodies(*bodies)
        else:
            bodies = bodies[0]

        cal_file.readline()  # Unused line.
        frequency_data = cal_file.readline().split()
        omega_range = np.linspace(float(frequency_data[1]), float(frequency_data[2]), int(frequency_data[0]))

        direction_data = cal_file.readline().split()
        direction_range = np.linspace(float(direction_data[1]), float(direction_data[2]), int(direction_data[0]))
        direction_range = np.pi/180*direction_range  # conversion from degrees to radians.

        # The options below are not implemented yet.

        cal_file.readline()  # Unused line.
        irf_data = cal_file.readline()
        show_pressure = cal_file.readline().split()[0] == "1"
        kochin_data = cal_file.readline().split()
        kochin_range = np.linspace(float(kochin_data[1]), float(kochin_data[2]), int(kochin_data[0]))
        free_surface_data = cal_file.readline().split()

    # Generate Capytaine's problem objects
    env_args = dict(body=bodies, rho=rho, sea_bottom=sea_bottom, g=g)
    problems = []
    for omega in omega_range:
        for direction in direction_range:
            problems.append(DiffractionProblem(wave_direction=direction, omega=omega, **env_args))
        for dof in bodies.dofs:
            problems.append(RadiationProblem(radiating_dof=dof, omega=omega, **env_args))

    return problems


def export_as_Nemoh_directory(problem, directory_name, omega_range=None):
    """Export radiation problems as Nemoh 2.0 directory (experimental).

    TODO: Diffraction problem.

    Parameters
    ----------
    problem : RadiationProblem
        the problem that should be exported
    directory_name : string
        path to the directory
    omega_range : list of float or array of float, optional
        the exported problem will be set up with the following linear range:
        linspace(min(omega_range), max(omega_range), len(omega_range))
    """

    if os.path.isdir(directory_name):
        LOG.warning(f"""Exporting problem in already existing directory: {directory_name}
             You might be overwriting existing files!""")
    else:
        os.makedirs(directory_name)

    # Export the mesh
    write_MAR(
        os.path.join(directory_name, f'{problem.body.name}.dat'),
        problem.body.mesh.vertices,
        problem.body.mesh.faces,
        # xOz_symmetry=isinstance(problem.body, ReflectionSymmetry)
    )

    # Set range of frequencies
    if omega_range is None:
        omega_nb_steps = 1
        omega_start = problem.omega
        omega_stop = problem.omega
    else:
        omega_nb_steps = len(omega_range)
        omega_start = min(omega_range)
        omega_stop = max(omega_range)

    # Write Nemoh.cal
    with open(os.path.join(directory_name, "Nemoh.cal"), "w") as nemoh_cal:
        nemoh_cal.write(
                DEFAULT_NEMOH_CAL.format(
                    rho=problem.rho,
                    g=problem.g,
                    depth=problem.depth if problem.depth < np.infty else 0,
                    mesh_filename=f'{problem.body.name}.dat',
                    mesh_vertices=problem.body.mesh.nb_vertices,
                    mesh_faces=problem.body.mesh.nb_faces,
                    omega_nb_steps=omega_nb_steps,
                    omega_start=omega_start,
                    omega_stop=omega_stop,
                    )
                )

    # Write input.txt
    with open(os.path.join(directory_name, "input.txt"), "w") as input_txt:
        input_txt.write(DEFAULT_INPUT_TXT)

    # Write ID.dat
    with open(os.path.join(directory_name, "ID.dat"), "w") as id_dat:
        id_dat.write(f"1\n.")


DEFAULT_NEMOH_CAL = """--- Environment ------------------------------------------------------------------------------------------------------------------
{rho}			! RHO			! KG/M**3	! Fluid specific volume
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
0	0.1	10.		! IRF				! IRF calculation (0 for no calculation), time step and duration
0				! Show pressure
0	0.	180.		! Kochin function		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)
0	0	100.	100.	! Free surface elevation	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction
"""

DEFAULT_INPUT_TXT = """--- Calculation parameters ------------------------------------------------------------------------------------
1				! Indiq_solver		! -		! Solver (0) Direct Gauss (1) GMRES (2) GMRES with FMM acceleration (2 not implemented yet)
20				! IRES			! -		! Restart parameter for GMRES
5.E-07				! TOL_GMRES		! -		! Stopping criterion for GMRES
100				! MAXIT			! -		! Maximum iterations for GMRES
1				! Sav_potential		! -		! Save potential for visualization
"""


def write_dataset_as_tecplot_files(results_directory, data):
    """Write some of the data from a xarray dataset into legacy tecplot file outputs."""

    if 'added_mass' in data:
        with open(os.path.join(results_directory, 'RadiationCoefficients.tec'), 'w') as fi:
            for i in range(len(data['radiating_dof'])+1):
                fi.write(f'...\n')
            for dof in data.radiating_dof:
                fi.write(f'{dof.values}\n')
                for o in data.omega:
                    fi.write(f'  {o.values:e}  ')
                    for dof2 in data.influenced_dof:
                        fi.write(f"{data['added_mass'].sel(omega=o, radiating_dof=dof, influenced_dof=dof2).values:e}")
                        fi.write('  ')
                        fi.write(f"{data['radiation_damping'].sel(omega=o, radiating_dof=dof, influenced_dof=dof2).values:e}")
                        fi.write('  ')
                    fi.write('\n')

    if 'diffraction_force' in data:
        data['excitation_force'] = data['Froude_Krylov_force'] + data['diffraction_force']
        with open(os.path.join(results_directory, 'ExcitationForce.tec'), 'w') as fi:
            for i in range(len(data.influenced_dof)+1):
                fi.write(f'...\n')
            for wave_direction in data.wave_direction.values:
                fi.write(f'angle={wave_direction}\n')
                for o in data.omega.values:
                    fi.write(f'  {o:e}  ')
                    for dof in data.influenced_dof:
                        val = data['excitation_force'].sel(omega=o, wave_direction=wave_direction, influenced_dof=dof).values
                        fi.write(f'{np.abs(val):e}')
                        fi.write('  ')
                        fi.write(f'{np.angle(val):e}')
                        fi.write('  ')
                    fi.write('\n')

