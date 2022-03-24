import logging
from pathlib import Path

import numpy as np
from numpy import pi

import capytaine as cpt
from capytaine.ui.vtk import Animation

try:
    import meshmagick.mesh as mm
    import meshmagick.hydrostatics as hs
except:
    hs = None
    mm = None

from scipy.linalg import block_diag

logging.basicConfig(level=logging.INFO, format='%(levelname)-8s: %(message)s')

bem_solver = cpt.BEMSolver()


def generate_boat() -> cpt.FloatingBody:
    boat = cpt.FloatingBody.from_file("boat_200.mar", file_format="mar", name="pirate ship")
    boat.center_of_mass = boat.center_of_buoyancy if hasattr(boat, 'center_of_buoyancy') else np.zeros(3)
    boat.rotate_z(pi)
    boat.add_all_rigid_body_dofs()
    boat.keep_immersed_part()

    # The computation of the RAO requires the values of the inertia matrix and the hydrostatic stiffness matrix.
    
    try:
        boat.hydrostatic_stiffness = boat.hydrostatic_stiffness_xr()
        boat.mass = boat.rigid_dof_mass()
    except:
        if hs is not None and mm is not None:
            # You can use Meshmagick to compute the hydrostatic stiffness matrix.
            mesh = mm.Mesh(boat.mesh.vertices, boat.mesh.faces, name=boat.mesh.name)
            hsd = hs.compute_hydrostatics(mesh, boat.center_of_mass, 1000, 9.80665)
    
            m = hsd['disp_mass']
            I = np.array([[hsd['Ixx'], -1*hsd['Ixy'], -1*hsd['Ixz']],
                          [-1*hsd['Ixy'], hsd['Iyy'], -1*hsd['Iyz']],
                          [-1*hsd['Ixz'], -1*hsd['Iyz'], hsd['Izz']]])
            M = block_diag(m, m, m, I)
            boat.mass = boat.add_dofs_labels_to_matrix(M)
    
            kHS = block_diag(0,0,hsd['stiffness_matrix'],0)
            boat.hydrostatic_stiffness = boat.add_dofs_labels_to_matrix(kHS)
    
        else:
            # Alternatively, you can define these by hand
            boat.mass = boat.add_dofs_labels_to_matrix(
            [[1e6, 0,   0,   0,   0,   0],
              [0,   1e6, 0,   0,   0,   0],
              [0,   0,   1e6, 0,   0,   0],
              [0,   0,   0,   1e7, 0,   2e5],
              [0,   0,   0,   0,   4e7, 0],
              [0,   0,   0,   2e5, 0,   5e7]]
            )
    
            boat.hydrostatic_stiffness = boat.add_dofs_labels_to_matrix(
                [[0, 0, 0,    0,   0,    0],
                  [0, 0, 0,    0,   0,    0],
                  [0, 0, 3e6,  0,   -7e6, 0],
                  [0, 0, 0,    2e7, 0,    0],
                  [0, 0, -7e6, 0,   1e8,  0],
                  [0, 0, 0,    0,   0,    0]]
            )
    return boat


def setup_animation(body, fs, omega, wave_amplitude, wave_direction) -> Animation:
    # SOLVE BEM PROBLEMS
    problems = [cpt.RadiationProblem(omega=omega, body=body, radiating_dof=dof) for dof in body.dofs]
    problems += [cpt.DiffractionProblem(omega=omega, body=body, wave_direction=wave_direction)]
    results = bem_solver.solve_all(problems)
    *radiation_results, diffraction_result = results
    dataset = cpt.assemble_dataset(results)

    # COMPUTE RAO
    dataset['RAO'] = cpt.post_pro.rao(dataset, wave_direction=wave_direction)

    # Compute the motion of each face of the mesh for the animation
    rao_faces_motion = sum(dataset['RAO'].sel(omega=omega, radiating_dof=dof).data * body.full_body.dofs[dof] for dof in body.dofs)

    # COMPUTE FREE SURFACE ELEVATION
    # Compute the diffracted wave pattern
    diffraction_elevation = bem_solver.get_free_surface_elevation(diffraction_result, fs)
    incoming_waves_elevation = fs.incoming_waves(diffraction_result)

    # Compute the wave pattern radiated by the RAO
    radiation_elevations_per_dof = {res.radiating_dof: bem_solver.get_free_surface_elevation(res, fs) for res in radiation_results}
    rao_radiation_elevation = sum(dataset['RAO'].sel(omega=omega, radiating_dof=dof).data * radiation_elevations_per_dof[dof] for dof in body.dofs)

    # SET UP ANIMATION
    animation = Animation(loop_duration=2*pi/omega)
    animation.add_body(body.full_body, faces_motion=wave_amplitude*rao_faces_motion)
    animation.add_free_surface(fs, wave_amplitude * (incoming_waves_elevation + diffraction_elevation + rao_radiation_elevation))
    return animation


if __name__ == '__main__':
    body = generate_boat()
    fs = cpt.FreeSurface(x_range=(-100, 100), y_range=(-100, 100), nx=100, ny=100)

    wave_direction = 0.0
    omega = 1.5

    anim = setup_animation(body, fs, omega=omega, wave_amplitude=0.5, wave_direction=wave_direction)
    anim.run(camera_position=(-60, -60, 90), resolution=(800, 600))

    filename = f"{body.name}__omega_{omega:.2f}__beta_{wave_direction:.2f}.ogv"
    filepath = Path.cwd() / filename
    # anim.save(str(filepath), camera_position=(60, 60, 90), resolution=(800, 600))

