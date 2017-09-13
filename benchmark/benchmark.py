#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import cProfile, pstats, io

import numpy as np
import pandas as pd

from capytaine.Nemoh import Nemoh
from capytaine.problems import RadiationProblem


def profile_capytaine(body, omega_range, result_dir):
    if not os.path.isdir(result_dir):
        os.makedirs(result_dir)

    pr = cProfile.Profile()
    pr.enable() #========================

    problems = [RadiationProblem(body=body, rho=1000, omega=omega) for omega in omega_range]
    solver = Nemoh()
    results = [solver.solve(pb) for pb in problems]

    pr.disable() #========================

    results = np.asarray(results)
    np.savetxt(f'{result_dir}/results.csv', results)

    s = io.StringIO()
    sortby = 'time'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    profiler_results = s.getvalue()
    with open(f'{result_dir}/profile.log', 'w') as log:
        log.write(profiler_results)

    return float(profiler_results.split('\n')[0].split('in')[1].strip('seconds\n'))


def profile_Nemoh(body, omega_range, result_dir, nemoh_bin_dir="~/nemoh/bin"):
    problem = RadiationProblem(body=body, rho=1000, omega=0.0)
    problem.export_as_Nemoh_directory(result_dir, omega_range)

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

    with open(f'{result_dir}/profile.log', 'w') as log:
        log.write(solver_return.stdout)

    return float(solver_return.stderr.split('\n')[0].strip('real'))

