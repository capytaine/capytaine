#!/usr/bin/env python
# coding: utf-8

import os
import logging
import subprocess
import cProfile, pstats, io

import numpy as np
import pandas as pd

from capytaine.Nemoh import Nemoh
from capytaine.problems import RadiationProblem
from capytaine.results import assemble_radiation_results_matrices
from capytaine.import_export import export_as_Nemoh_directory


def profile_capytaine(body, omega_range, result_dir, **problem_kwargs):
    if not os.path.isdir(result_dir):
        os.makedirs(result_dir)

    os.environ["MKL_NUM_THREADS"] = "1"

    if logging.root:
        del logging.root.handlers[:]

    logging.basicConfig(
        filename=f"{result_dir}/capytaine.log",
        level=logging.DEBUG,
        format="%(levelname)s:\t%(message)s"
    )

    pr = cProfile.Profile()
    pr.enable() #==Start profiler==

    problems = [RadiationProblem(body=body, omega=omega, **problem_kwargs) for omega in omega_range]
    solver = Nemoh()
    results = [solver.solve(pb) for pb in problems]

    pr.disable() #=================

    added_mass, dampings = assemble_radiation_results_matrices(results)
    np.savetxt(f'{result_dir}/results.csv', np.c_[added_mass.values[:, 0, 0], dampings.values[:, 0, 0]])

    s = io.StringIO()
    sortby = 'time'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    profiler_results = s.getvalue()
    with open(f'{result_dir}/profile.log', 'w') as log:
        log.write(profiler_results)

    os.environ["MKL_NUM_THREADS"] = "4"

    return float(profiler_results.split('\n')[0].split('in')[1].strip('seconds\n'))


def profile_Nemoh(body, omega_range, result_dir, nemoh_bin_dir="~/work/code/nemoh/bin", **problem_args):
    """Use Nemoh 2.0 to solve a problem and mesure computation time."""
    problem = RadiationProblem(body=body, omega=0.0, **problem_args)
    export_as_Nemoh_directory(problem, result_dir, omega_range)

    subprocess.run(
        f'cd {result_dir} && ' + os.path.join(nemoh_bin_dir, 'preProc'),
        shell=True,
        stdout=subprocess.PIPE,
    )
    solver_return = subprocess.run(
        f'cd {result_dir} && /usr/bin/time -p ' + os.path.join(nemoh_bin_dir, 'solver'),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding='utf8'
    )
    subprocess.run(
        f'cd {result_dir} && ' + os.path.join(nemoh_bin_dir, 'postProc'),
        shell=True,
        stdout=subprocess.PIPE,
    )

    with open(f'{result_dir}/profile.log', 'w') as log:
        log.write(solver_return.stdout)

    return float(solver_return.stderr.split('\n')[0].strip('real'))

