import timeit
from matplotlib.pylab import *
import torch
import numpy


### Define the custom solver function
def solve_numpy_cpu(a, b):
    numpy.linalg.solve(a, b)


## This function is to test baseline torch performance
def solve_torch_cpu(a, b):
    a = torch.from_numpy(a).cfloat()
    b = torch.from_numpy(b).cfloat()
    torch.linalg.solve(a, b)


## This function is to test GPU performance
def solve_torch_gpu(a, b):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    a = torch.from_numpy(a).cfloat().to(device)
    b = torch.from_numpy(b).cfloat().to(device)

    res = torch.linalg.solve(a, b)


results = {}

# Define the size of the matrix to test
dims = [5, 10, 50, 100, 500, 1000, 5000, 10000, 20000, 30000, 40000]
# our test functions
funcs = (solve_numpy_cpu, solve_torch_cpu, solve_torch_gpu)
func_dict = {f.__name__: f for f in funcs}

times = 5


def main():
    """run though each dimension and each function, comparing the time taken."""
    for dim in dims:
        a = numpy.random.rand(dim, dim)
        b = numpy.random.rand(dim)

        for fname, f in func_dict.items():
            time = timeit.timeit(lambda: f(a, b), number=times)

            if fname in results:
                results[fname].append(time)
            else:
                results[fname] = [time]

            print(f"{f.__name__}: {dim:<10}|{time:f}")


main()


fig, ax = subplots()
for fname, res in results.items():
    ax.plot(dims[: len(res)], res, label=fname)
ax.grid()
ax.legend()
ax.set(yscale="log", xscale="log", xlabel="Dimension", ylabel="Time (s)")
