#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':

    # look for csv files in dir
    results = []
    for res_file in os.listdir('benchmarks/openmp'):
        if len(res_file) < 5 or res_file[-4:] != '.csv':
            continue
        df = pd.read_csv(os.path.join('benchmarks/openmp', res_file))
        df["filename"] = res_file
        results.append(df)

    # create dataframe
    df = pd.concat(results, axis="index")
    df.columns = df.columns.str.strip()
    df["kind"] = df["kind"].str.strip()
    df.reset_index(drop=True, inplace=True)

    df["branch"] = ''
    df["hash"] = ''
    for i in df.index:
        tmp1 = df.loc[i, "filename"][14:-4]
        tmp2 = tmp1.find('_')
        df.loc[i, "branch"] = tmp1[tmp2+1:]
        df.loc[i, "hash"] = tmp1[:tmp2]

    # plot

    mk = {'full': 'o', 'wave_only': 's', 'half_wave_only': '^'}
    hs = df['hash'].unique()
    cl = {}
    for i in range(len(hs)):
        cl[hs[i]] = 'C%1d' % (i % 10,)
    del(hs)


    plt.figure()
    for (commit, group_df) in df.groupby(["kind", "hash"]):
        plt.plot(group_df["n_threads"], group_df["elapsed_time"],
                 label='%s, %s' % (commit[1], commit[0]), ls='-',
                 marker=mk[commit[0]], c=cl[commit[1]])
    plt.xlabel("Number of threads")
    plt.ylabel("Computation time")
    plt.yscale("log")
    plt.title("OpenMP parallelization")
    plt.grid(True)
    plt.legend()

    plt.figure()
    n_threads_range = np.arange(df.n_threads.min(), df.n_threads.max()+1)
    plt.plot(n_threads_range, n_threads_range, '--', c='gray')
    for (commit, group_df) in df.groupby(["kind", "hash"]):
        ref_time = group_df.loc[group_df["n_threads"] == 1]["elapsed_time"].values[0]
        plt.plot(group_df["n_threads"], ref_time / group_df["elapsed_time"],
                 label='%s, %s' % (commit[1], commit[0]), ls='-',
                 marker=mk[commit[0]], c=cl[commit[1]])
    plt.xlabel("Number of threads")
    plt.ylabel("Speedup")
    plt.title("OpenMP parallelization (2)")
    plt.grid(True)
    plt.legend()

    plt.figure()
    for (commit, group_df) in df.groupby(["kind", "hash"]):
        ref_time = group_df.loc[group_df["n_threads"] == 1]["elapsed_time"].values[0]
        plt.plot(group_df["n_threads"], ref_time / group_df["elapsed_time"] / group_df["n_threads"],
                 label='%s, %s' % (commit[1], commit[0]), ls='-',
                 marker=mk[commit[0]], c=cl[commit[1]])
    plt.xlabel("Number of threads")
    plt.ylabel("Efficiency")
    plt.title("OpenMP parallelization (3)")
    plt.grid(True)
    plt.legend()

    plt.show()
