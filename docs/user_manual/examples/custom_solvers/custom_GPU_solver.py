import numpy as np
import capytaine as cpt
import xarray as xr

## Custom linear solver setup
### In this example, we use Pytorch to compute the LU decomposition
### and solve the linear system on the GPU.
### You can use any other library to do the same.
import torch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def lu_linear_solver_with_cache_on_GPU(A, b):
    """our custom solver function takes a matrix A and a vector b
    and returns the solution of the linear system Ax = b."""

    A_on_GPU = torch.from_numpy(np.array(A)).cfloat().to(device)
    # Compute the LU decomposition, which is handy for recalculating the solution, with various B values for the same A matrix
    latest_LU_decomp = torch.linalg.lu_factor(A_on_GPU)
    print("Computing LU decomposition with Pytorch")

    # Its important to use the appropriate mechnism to transfer data to the GPU
    b_on_GPU = torch.from_numpy(b.reshape(-1, 1)).cfloat().to(device)
    res = torch.linalg.lu_solve(*latest_LU_decomp, b_on_GPU)
    return res.cpu().numpy()  # transfer the result back to CPU


# define the linear solver, with your custom solver function
solver = cpt.BEMSolver(
    engine=cpt.BasicMatrixEngine(linear_solver=lu_linear_solver_with_cache_on_GPU)
)

# Your test case
mesh = cpt.mesh_sphere().immersed_part()
body = cpt.FloatingBody(mesh, dofs=cpt.rigid_body_dofs())

# define the solution grid
test_matrix = xr.Dataset(
    coords={
        "omega": np.linspace(0.1, 4.0, 10),
        "radiating_dof": list(body.dofs),
    }
)

# your solution!
ds = solver.fill_dataset(test_matrix, body)
