from numpy import pi

import capytaine as cpt
from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
from capytaine.ui.vtk import Animation

cpt.set_logging('INFO')

bem_solver = cpt.BEMSolver()


def generate_boat():
    boat_mesh = cpt.load_mesh("boat_200.mar", file_format="mar")
    boat = cpt.FloatingBody(
            mesh=boat_mesh,
            dofs=cpt.rigid_body_dofs(rotation_center=boat_mesh.center_of_buoyancy),
            center_of_mass = boat_mesh.center_of_buoyancy,
            name="pirate ship"
            )
    boat.inertia_matrix = boat.compute_rigid_body_inertia() / 10 # Artificially lower to have a more appealing animation
    boat.hydrostatic_stiffness = boat.immersed_part().compute_hydrostatic_stiffness()
    return boat


def setup_animation(body, fs, omega, wave_amplitude, wave_direction):
    # SOLVE BEM PROBLEMS
    radiation_problems = [cpt.RadiationProblem(omega=omega, body=body.immersed_part(), radiating_dof=dof) for dof in body.dofs]
    radiation_results = bem_solver.solve_all(radiation_problems)
    diffraction_problem = cpt.DiffractionProblem(omega=omega, body=body.immersed_part(), wave_direction=wave_direction)
    diffraction_result = bem_solver.solve(diffraction_problem)

    dataset = cpt.assemble_dataset(radiation_results + [diffraction_result])
    rao = cpt.post_pro.rao(dataset, wave_direction=wave_direction)

    # COMPUTE FREE SURFACE ELEVATION
    # Compute the diffracted wave pattern
    incoming_waves_elevation = airy_waves_free_surface_elevation(fs, diffraction_result)
    diffraction_elevation = bem_solver.compute_free_surface_elevation(fs, diffraction_result)

    # Compute the wave pattern radiated by the RAO
    radiation_elevations_per_dof = {res.radiating_dof: bem_solver.compute_free_surface_elevation(fs, res) for res in radiation_results}
    radiation_elevation = sum(rao.sel(omega=omega, radiating_dof=dof).data * radiation_elevations_per_dof[dof] for dof in body.dofs)

    # SET UP ANIMATION
    # Compute the motion of each face of the mesh for the animation
    rao_faces_motion = sum(rao.sel(omega=omega, radiating_dof=dof).data * body.dofs[dof] for dof in body.dofs)

    # Set up scene
    animation = Animation(loop_duration=2*pi/omega)
    animation.add_body(body, faces_motion=wave_amplitude*rao_faces_motion)
    animation.add_free_surface(fs, wave_amplitude * (incoming_waves_elevation + diffraction_elevation + radiation_elevation))
    return animation


if __name__ == '__main__':
    body = generate_boat()
    fs = cpt.FreeSurface(x_range=(-100, 75), y_range=(-100, 75), nx=100, ny=100)

    anim = setup_animation(body, fs, omega=1.5, wave_amplitude=0.5, wave_direction=pi)
    anim.run(camera_position=(70, 70, 100), resolution=(800, 600))
    anim.save("animated_boat.ogv", camera_position=(70, 70, 100), resolution=(800, 600))
