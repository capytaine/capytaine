#!/usr/bin/env python
# coding: utf-8
"""Experimental command-line interface for Capytaine."""
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
