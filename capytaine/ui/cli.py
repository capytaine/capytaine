#!/usr/bin/env python
# coding: utf-8
"""Experimental command-line interface for Capytaine."""
# Copyright (C) 2017-2023 Matthieu Ancellin
# See LICENSE file at <https://github.com/capytaine/capytaine>

import argparse

import capytaine as cpt
from capytaine.io.legacy import run_cal_file

cpt.set_logging()

parser = argparse.ArgumentParser(description="Command-line interface for Capytaine taking Nemoh.cal files as input and returning Tecplots files.")
parser.add_argument('paramfiles',
                    default=['./Nemoh.cal'],
                    nargs='*',
                    help='path of parameters files (default: ./Nemoh.cal)')


def main():
    args = parser.parse_args()
    for paramfile in args.paramfiles:
        run_cal_file(paramfile)


if __name__ == '__main__':
    main()
