
import os
import logging
from capytaine.io.mesh_writers import write_MAR

LOG = logging.getLogger(__name__)


def export_as_Nemoh_directory(problem, directory_name, omega_range=None):
    """Export radiation problems as Nemoh 2.0 directory (experimental).

    TODO: Diffraction problem.

    Parameters
    ----------
    problem : RadiationProblem
        the problem that should be exported
    directory_name : string
        path to the directory
    omega_range : list of float or array of float, optional
        the exported problem will be set up with the following linear range:
        linspace(min(omega_range), max(omega_range), len(omega_range))
    """

    if os.path.isdir(directory_name):
        LOG.warning(f"""Exporting problem in already existing directory: {directory_name}
             You might be overwriting existing files!""")
    else:
        os.makedirs(directory_name)

    # Export the mesh
    write_MAR(
        os.path.join(directory_name, f'{problem.body.name}.dat'),
        problem.body.mesh.vertices,
        problem.body.mesh.faces,
        # xOz_symmetry=isinstance(problem.body, ReflectionSymmetry)
    )

    # Set range of frequencies
    if omega_range is None:
        omega_nb_steps = 1
        omega_start = problem.omega
        omega_stop = problem.omega
    else:
        omega_nb_steps = len(omega_range)
        omega_start = min(omega_range)
        omega_stop = max(omega_range)

    # Write Nemoh.cal
    with open(os.path.join(directory_name, "Nemoh.cal"), "w") as nemoh_cal:
        nemoh_cal.write(
                DEFAULT_NEMOH_CAL.format(
                    rho=problem.rho,
                    g=problem.g,
                    depth=problem.water_depth if problem.water_depth < np.inf else 0,
                    mesh_filename=f'{problem.body.name}.dat',
                    mesh_vertices=problem.body.mesh.nb_vertices,
                    mesh_faces=problem.body.mesh.nb_faces,
                    omega_nb_steps=omega_nb_steps,
                    omega_start=omega_start,
                    omega_stop=omega_stop,
                    )
                )

    # Write input.txt
    with open(os.path.join(directory_name, "input.txt"), "w") as input_txt:
        input_txt.write(DEFAULT_INPUT_TXT)

    # Write ID.dat
    with open(os.path.join(directory_name, "ID.dat"), "w") as id_dat:
        id_dat.write(f"1\n.")


DEFAULT_NEMOH_CAL = """--- Environment ------------------------------------------------------------------------------------------------------------------
{rho}			! RHO			! KG/M**3	! Fluid specific volume
{g}				! G			! M/S**2	! Gravity
{depth}			! DEPTH			! M		! Water depth
0.	0.			! XEFF YEFF		! M		! Wave measurement point
--- Description of floating bodies -----------------------------------------------------------------------------------------------
1				! Number of bodies
--- Body 1 -----------------------------------------------------------------------------------------------------------------------
{mesh_filename}
{mesh_vertices} {mesh_faces}
1				! Number of degrees of freedom
1 0. 0. 1. 0. 0. 0.		! Heave
1				! Number of resulting generalised forces
1 0. 0. 1. 0. 0. 0.		! Heave
0				! Number of lines of additional information
--- Load cases to be solved -------------------------------------------------------------------------------------------------------
{omega_nb_steps} {omega_start} {omega_stop} ! Frequencies range
0	0.	0.		! Number of wave directions, Min and Max (degrees)
--- Post processing ---------------------------------------------------------------------------------------------------------------
0	0.1	10.		! IRF				! IRF calculation (0 for no calculation), time step and duration
0				! Show pressure
0	0.	180.		! Kochin function		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)
0	0	100.	100.	! Free surface elevation	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction
"""

DEFAULT_INPUT_TXT = """--- Calculation parameters ------------------------------------------------------------------------------------
1				! Indiq_solver		! -		! Solver (0) Direct Gauss (1) GMRES (2) GMRES with FMM acceleration (2 not implemented yet)
20				! IRES			! -		! Restart parameter for GMRES
5.E-07				! TOL_GMRES		! -		! Stopping criterion for GMRES
100				! MAXIT			! -		! Maximum iterations for GMRES
1				! Sav_potential		! -		! Save potential for visualization
"""
