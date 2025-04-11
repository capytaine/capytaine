import time
import numpy as np
import capytaine as cpt
import xarray as xr
import torch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

reference_cpu_solver = cpt.BEMSolver()


# A simple implementation of a linear solver on GPU using PyTorch
def linear_solver_on_GPU(A, b):
    A_on_GPU = torch.tensor(np.array(A), dtype=torch.complex64, device=device)
    b_on_GPU = torch.tensor(b, dtype=torch.complex64, device=device)
    res = torch.linalg.solve(A_on_GPU, b_on_GPU)
    return res.cpu().numpy()

naive_gpu_solver = cpt.BEMSolver(
    engine=cpt.BasicMatrixEngine(linear_solver=linear_solver_on_GPU)
)


# Alternatively, an optimization of the above to avoid doing the same GPU
# transfer and LU factorization twice
latest_A = None
latest_LU_decomp = None
def lu_linear_solver_with_cache_on_GPU(A, b):
    global latest_A, latest_LU_decomp
    if A is not latest_A:
        latest_A = A
        A_on_GPU = torch.tensor(np.array(A), dtype=torch.complex64, device=device)
        latest_LU_decomp = torch.linalg.lu_factor(A_on_GPU)
    b_on_GPU = torch.tensor(np.array(b), dtype=torch.complex64, device=device).reshape(-1, 1)
    # Reshaping because `torch.linalg.lu_solve` wants a column vector on the
    # right-hand-side (unlike `torch.linalg.solve`).
    res = torch.linalg.lu_solve(*latest_LU_decomp, b_on_GPU).reshape(-1)
    return res.cpu().numpy()

gpu_solver = cpt.BEMSolver(
    engine=cpt.BasicMatrixEngine(linear_solver=lu_linear_solver_with_cache_on_GPU)
)


###############
#  Benchmark  #
###############

mesh = cpt.mesh_sphere(resolution=(50, 50)).immersed_part()
print(mesh.nb_faces, "faces")
body = cpt.FloatingBody(mesh, dofs=cpt.rigid_body_dofs())

test_matrix = xr.Dataset(
    coords={
        "omega": np.linspace(0.0, 4.0, 20),
        "radiating_dof": list(body.dofs),
        "water_depth": [np.inf],
    }
)

start = time.perf_counter()
ds1 = reference_cpu_solver.fill_dataset(test_matrix, body)
print("CPU:", time.perf_counter() - start)

start = time.perf_counter()
ds2 = naive_gpu_solver.fill_dataset(test_matrix, body)
print("GPU:", time.perf_counter() - start)

start = time.perf_counter()
ds3 = gpu_solver.fill_dataset(test_matrix, body)
print("Slightly faster GPU:", time.perf_counter() - start)


# Check outputs
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(
    ds1.omega.values,
    ds1.added_mass.sel(radiating_dof="Roll", influenced_dof="Roll").values,
)
ax.plot(
    ds2.omega.values,
    ds2.added_mass.sel(radiating_dof="Roll", influenced_dof="Roll").values,
)
ax.plot(
    ds3.omega.values,
    ds3.added_mass.sel(radiating_dof="Roll", influenced_dof="Roll").values,
)
plt.show()
