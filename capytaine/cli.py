#!/usr/bin/env python
# coding: utf-8
"""
Command-line interface for capytaine.

This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.
"""

import argparse
import logging
import os

import numpy as np

from capytaine.tools.import_export import import_cal_file
from capytaine.results import assemble_dataset
from capytaine.Nemoh import Nemoh

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")

LOG = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="Command-line interface for the BEM solver Nemoh.")
parser.add_argument('paramfiles',
                    default=['./Nemoh.cal'],
                    nargs='*',
                    help='path of parameters files (default: ./Nemoh.cal)')


def main():
    args = parser.parse_args()
    for paramfile in args.paramfiles:
        problems = import_cal_file(paramfile)
        solver = Nemoh()
        results = [solver.solve(pb) for pb in problems]
        data = assemble_dataset(results)
        print(data)

        results_directory = os.path.join(os.path.dirname(paramfile), 'results')
        try:
            os.mkdir(results_directory)
        except FileExistsError:
            LOG.warning("The 'results' directory already exists. You might be overwriting existing data.")

        LOG.info("Write radiation coefficients in legacy tecplot format.")
        if 'added_mass' in data:
            with open(os.path.join(results_directory, 'RadiationCoefficients.tec'), 'w') as fi:
                for i in range(len(data['radiating_dof'])+1):
                    fi.write(f'...\n')
                for dof in data.influenced_dof:
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
            LOG.info("Write excitation forces in legacy tecplot format.")
            with open(os.path.join(results_directory, 'ExcitationForce.tec'), 'w') as fi:
                for i in range(len(data.influenced_dof)+1):
                    fi.write(f'...\n')
                for angle in data.angle.values:
                    fi.write(f'angle={angle}\n')
                    for o in data.omega.values:
                        fi.write(f'  {o:e}  ')
                        for dof in data.influenced_dof:
                            val = data['excitation_force'].sel(omega=o, angle=angle, influenced_dof=dof).values
                            fi.write(f'{np.abs(val):e}')
                            fi.write('  ')
                            fi.write(f'{np.angle(val):e}')
                            fi.write('  ')
                        fi.write('\n')


if __name__ == '__main__':
    main()

