# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Import or export Nemoh data files for backward compatibility with Nemoh 2."""

import os
import logging
from pathlib import Path
from typing import Union, List

import numpy as np

from capytaine.bem.solver import BEMSolver
from capytaine.io.xarray import assemble_dataset
from capytaine.meshes.io import load_mesh
from capytaine.bodies.bodies import FloatingBody
from capytaine.bodies.multibodies import Multibody
from capytaine.bem.problems_and_results import DiffractionProblem, RadiationProblem

LOG = logging.getLogger(__name__)


def import_cal_file(filepath):
    """Read a Nemoh.cal file and return a list of problems."""

    with open(filepath, 'r') as cal_file:

        cal_file.readline()  # Unused line.
        rho = float(cal_file.readline().split()[0])
        g = float(cal_file.readline().split()[0])
        water_depth = float(cal_file.readline().split()[0])
        if water_depth == 0.0:
            water_depth = np.inf
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
                body = FloatingBody(mesh=load_mesh(mesh_file, file_format='nemoh'))

            nb_dofs = int(cal_file.readline().split()[0])
            for i_dof in range(nb_dofs):
                dof_data = cal_file.readline().split()
                if int(dof_data[0]) == 1:
                    direction = np.array([float(x) for x in dof_data[1:4]])
                    body.add_translation_dof(direction=direction)
                elif int(dof_data[0]) == 2:
                    direction = np.array([float(x) for x in dof_data[1:4]])
                    center_of_mass = np.array([float(x) for x in dof_data[4:7]])
                    body.add_rotation_dof(direction=direction, rotation_center=center_of_mass)

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
        frequency_data_string_without_comment = cal_file.readline().split('!')[0]
        frequency_data = frequency_data_string_without_comment.split()
        if len(frequency_data) == 3:  # Nemoh v2 format
            omega_range = np.linspace(float(frequency_data[1]), float(frequency_data[2]), int(frequency_data[0]))
        else:
            type_of_frequency_data = int(frequency_data[0])
            if type_of_frequency_data == 1:  # angular frequency
                omega_range = np.linspace(float(frequency_data[2]), float(frequency_data[3]), int(frequency_data[1]))
            elif type_of_frequency_data == 2:  # frequency
                omega_range = 2*np.pi*np.linspace(float(frequency_data[2]), float(frequency_data[3]), int(frequency_data[1]))
            elif type_of_frequency_data == 3:  # period
                omega_range = 2*np.pi/np.linspace(float(frequency_data[2]), float(frequency_data[3]), int(frequency_data[1]))
            else:
                raise ValueError(f"Cannot parse the frequency data \"{frequency_data_string_without_comment}\" in {filepath}.")


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
    env_args = dict(body=bodies, rho=rho, water_depth=water_depth, g=g)
    problems = []
    for omega in omega_range:
        for direction in direction_range:
            problems.append(DiffractionProblem(wave_direction=direction, omega=omega, **env_args))
        for dof in bodies.dofs:
            problems.append(RadiationProblem(radiating_dof=dof, omega=omega, **env_args))

    return problems


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


def _hydrostatics_writer(
    hydrostatics_file_path,
    kh_file_path,
    center_of_buoyancy,
    center_of_mass,
    disp_volume,
    hydrostatic_stiffness
):
    """Write the Hydrostatics.dat and KH.dat files"""
    with open(hydrostatics_file_path, 'w') as hf:
        for j in range(3):
            line =  f'XF = {center_of_buoyancy[j]:7.4f} - XG = {center_of_mass[j]:7.4f} \n'
            hf.write(line)
        line = f'Displacement = {disp_volume:1.6E}\n'
        hf.write(line)
        hf.close()
    np.savetxt(kh_file_path, hydrostatic_stiffness.values, fmt='%1.6E')


def export_hydrostatics(
        hydrostatics_directory: str,
        bodies: Union[FloatingBody, Multibody, List[FloatingBody]],
) -> None:
    """Export rigid body hydrostatics in Nemoh's format (KH.dat and Hydrostatics.dat).
    If the bodies have hydrostatics matrices already defined, uses them.
    Otherwise, tires to compute the hydrostatics stiffness.

    Parameters
    ----------
    hydrostatics_directory: string
        Path to the directory in which the data will be written (two files per body)
    bodies: FloatingBody or Multibody
        The body or bodies, which are all expected to be rigid bodies with 6 dofs.
        A list of FloatingBody is also accepted for backward compatibility.

    Return
    ------
    None
    """

    if os.path.isdir(hydrostatics_directory):
        LOG.warning(f"""Exporting problem in already existing directory: {hydrostatics_directory}
             You might be overwriting existing files!""")
    else:
        os.makedirs(hydrostatics_directory)

    if isinstance(bodies, FloatingBody):
        bodies = [bodies]
    elif isinstance(bodies, Multibody):
        if not all(hasattr(b, 'hydrostatic_stiffness') for b in bodies.bodies):
            raise ValueError("All the bodies in the multibody should have a hydrostatic stiffness matrix defined separately to be able to export the hydrostatics with `io.legacy.export_hydrostatics`.")
        bodies = bodies.bodies

    body_count = len(bodies)
    if body_count == 1:
        body = bodies[0]
        hydrostatics_file_path = os.path.join(hydrostatics_directory, "Hydrostatics.dat")
        kh_file_path = os.path.join(hydrostatics_directory, "KH.dat")
        _hydrostatics_writer(
            hydrostatics_file_path,
            kh_file_path,
            body.center_of_buoyancy,
            body.center_of_mass,
            body.disp_volume,
            body.hydrostatic_stiffness
        )
    else:
        for (i, body) in enumerate(bodies):
            hydrostatics_file_path = os.path.join(hydrostatics_directory, f"Hydrostatics_{i}.dat")
            kh_file_path = os.path.join(hydrostatics_directory, f"KH_{i}.dat")
            _hydrostatics_writer(
                hydrostatics_file_path,
                kh_file_path,
                body.center_of_buoyancy,
                body.center_of_mass,
                body.disp_volume,
                body.hydrostatic_stiffness
            )


def export_hydrostatics_from_dataset(hydrostatics_directory: Union[str, Path], dataset):
    """
    """
    if os.path.isdir(hydrostatics_directory):
        LOG.warning(f"""Exporting problem in already existing directory: {hydrostatics_directory}
             You might be overwriting existing files!""")
    else:
        os.makedirs(hydrostatics_directory)

    bodies = dataset.coords["body"]
    if bodies.shape == () or bodies.shape == (1,):
        hydrostatics_file_path = os.path.join(hydrostatics_directory, "Hydrostatics.dat")
        kh_file_path = os.path.join(hydrostatics_directory, "KH.dat")
        _hydrostatics_writer(
            hydrostatics_file_path,
            kh_file_path,
            dataset["center_of_buoyancy"].values,
            dataset["center_of_mass"].values,
            dataset["disp_mass"].values/dataset.coords["rho"].values,
            dataset["hydrostatic_stiffness"],
        )
    else:
        for (i, body_name) in enumerate(bodies.values):
            hydrostatics_file_path = os.path.join(hydrostatics_directory, f"Hydrostatics_{i}.dat")
            kh_file_path = os.path.join(hydrostatics_directory, f"KH_{i}.dat")
            body_influenced_dofs = list(
                dof for dof in dataset.coords["influenced_dof"].values if dof.startswith(body_name)
            )
            body_radiating_dofs = list(
                dof for dof in dataset.coords["radiating_dof"].values if dof.startswith(body_name)
            )
            _hydrostatics_writer(
                hydrostatics_file_path,
                kh_file_path,
                dataset["center_of_buoyancy"].sel(body=body_name).values,
                dataset["center_of_mass"].sel(body=body_name).values,
                dataset["disp_mass"].sel(body=body_name).values/dataset.coords["rho"].values,
                dataset["hydrostatic_stiffness"].sel(influenced_dof=body_influenced_dofs,
                                                     radiating_dof=body_radiating_dofs),
            )


def run_cal_file(paramfile):
    problems = import_cal_file(paramfile)
    solver = BEMSolver()
    results = solver.solve_all(problems)
    data = assemble_dataset(results)

    results_directory = os.path.join(os.path.dirname(paramfile), 'results')
    try:
        os.mkdir(results_directory)
    except FileExistsError:
        LOG.warning(f"The output directory ({results_directory}) already exists. You might be overwriting existing data.")

    LOG.info("Write results in legacy tecplot format.")
    write_dataset_as_tecplot_files(results_directory, data)
