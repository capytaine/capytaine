#!/usr/bin/env python
# coding: utf-8
"""
Command-line interface for capytaine.
"""

import argparse
import logging
from itertools import chain

import numpy as np
import matplotlib.pyplot as plt

from capytaine.import_export import import_cal_file
from capytaine.results import assemble_radiation_results_matrices
from capytaine.Nemoh import Nemoh

logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")

parser = argparse.ArgumentParser(description="Command-line interface for the BEM solver Nemoh.")
parser.add_argument('paramfiles',
                    default=['./Nemoh.cal'],
                    nargs='*',
                    help='path of parameters files (default: ./Nemoh.cal)')


def main():
    args = parser.parse_args()
    problems = chain(*(import_cal_file(paramfile) for paramfile in args.paramfiles))
    solver = Nemoh()
    results = [solver.solve(pb) for pb in problems]
    print(assemble_radiation_results_matrices(results))


if __name__ == '__main__':
    main()

