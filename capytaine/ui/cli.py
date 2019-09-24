#!/usr/bin/env python
# coding: utf-8
"""Experimental command-line interface for Capytaine."""
# Copyright (C) 2017-2019 Matthieu Ancellin
# See LICENSE file at <https://github.com/mancellin/capytaine>

import argparse
import logging
import os

from capytaine.io.legacy import import_cal_file, write_dataset_as_tecplot_files
from capytaine.io.xarray import assemble_dataset
from capytaine.bem.solver import BEMSolver

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
        solver = BEMSolver()
        results = [solver.solve(pb) for pb in problems]
        data = assemble_dataset(results)
        print(data)

        results_directory = os.path.join(os.path.dirname(paramfile), 'results')
        try:
            os.mkdir(results_directory)
        except FileExistsError:
            LOG.warning("The 'results' directory already exists. You might be overwriting existing data.")

        LOG.info("Write results in legacy tecplot format.")
        write_dataset_as_tecplot_files(results_directory, data)


if __name__ == '__main__':
    main()

