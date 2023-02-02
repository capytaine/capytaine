from glob import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

results = []
for res_file in glob("results/*/omp.csv"):
    df = pd.read_csv(res_file)
    df["filename"] = res_file
    results.append(df)

df = pd.concat(results, axis="index")
df.columns = df.columns.str.strip()
df["kind"] = df["kind"].str.strip()
df["commit_hash"] = df["filename"].str.extract("/(.*)/")

ref_time = df.loc[df["n_threads"] == 1]["elapsed_time"].loc[0]

# df["speedup"] = ref_time/df["elapsed_time"]
# df["efficiency"] = df["speedup"]/df["n_threads"]

plt.figure()
for (commit, group_df) in df.groupby(["kind", "commit_hash"]):
    plt.plot(group_df["n_threads"], group_df["elapsed_time"], label=commit)
n_threads_range = np.arange(df.n_threads.min(), df.n_threads.max()+1)
# plt.plot(n_threads_range, 1/n_threads_range, linestyle="--", color="grey")
plt.xlabel("Number of threads")
plt.ylabel("Computation time")
plt.yscale("log")
plt.title("OpenMP parallelization")
plt.legend()
plt.show()
