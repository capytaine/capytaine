#!/usr/bin/env python
# coding: utf-8
"""
Command-line interface for capytaine.
"""

import argparse
from itertools import chain
import logging
import os

import numpy as np
import matplotlib.pyplot as plt

from capytaine.import_export import import_cal_file
from capytaine.results import assemble_radiation_results_matrices, assemble_diffraction_results
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
        added_mass, radiation_damping = assemble_radiation_results_matrices(results)
        forces = assemble_diffraction_results(results)

        results_directory = os.path.join(os.path.dirname(paramfile), 'results')
        try:
            os.mkdir(results_directory)
        except FileExistsError:
            LOG.warning("The 'results' directory already exists. You might be overwriting existing data.")

        LOG.info("Write radiation coefficients in legacy tecplot format.")
        with open(os.path.join(results_directory, 'RadiationCoefficients.tec'), 'w') as fi:
            for i in range(len(added_mass.radiating_dof)+1):
                fi.write(f'...\n')
            for dof in added_mass.influenced_dof:
                fi.write(f'{dof.values}\n')
                for o in added_mass.omega:
                    fi.write(f'  {o.values:e}  ')
                    for dof2 in added_mass.influenced_dof:
                        fi.write(f'{added_mass.sel(omega=o, radiating_dof=dof, influenced_dof=dof2).values:e}')
                        fi.write('  ')
                        fi.write(f'{radiation_damping.sel(omega=o, radiating_dof=dof, influenced_dof=dof2).values:e}')
                        fi.write('  ')
                    fi.write('\n')

        LOG.info("Write excitation forces in legacy tecplot format.")
        with open(os.path.join(results_directory, 'ExcitationForce.tec'), 'w') as fi:
            for i in range(len(forces.influenced_dof)+1):
                fi.write(f'...\n')
            for angle in forces.angle.values:
                fi.write(f'angle={angle}\n')
                for o in forces.omega.values:
                    fi.write(f'  {o:e}  ')
                    for dof in forces.influenced_dof.values:
                        val = forces.sel(omega=o, angle=angle, influenced_dof=dof).values
                        fi.write(f'{o*abs(val):e}')
                        fi.write('  ')
                        fi.write(f'{np.angle(val):e}')
                        fi.write('  ')
                    fi.write('\n')

if __name__ == '__main__':
    main()

