import logging
from pathlib import Path

import numpy as np
from numpy import pi

import capytaine as cpt
from capytaine.ui.vtk import Animation

logging.basicConfig(level=logging.INFO, format='%(levelname)-8s: %(message)s')

bem_solver = cpt.BEMSolver()


def generate_boat() -> cpt.FloatingBody:
    boat = cpt.FloatingBody.from_file("boat_200.mar", file_format="mar", name="pirate ship")
    boat.rotate_z(pi)
    boat.center_of_mass = boat.center_of_buoyancy if hasattr(boat, 'center_of_buoyancy') else np.zeros(3)
    boat.rotation_center = boat.center_of_mass
    boat.add_all_rigid_body_dofs()
    boat.keep_immersed_part()

    # Compute hydrostatics
    boat.hydrostatic_stiffness = boat.compute_hydrostatic_stiffness()
    boat.inertia_matrix = boat.compute_rigid_body_inertia()
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

