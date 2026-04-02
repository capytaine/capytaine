from numpy import pi

import capytaine as cpt
from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
from capytaine.ui.vedo_animations import Animation


bem_solver = cpt.BEMSolver()


def generate_boat():
    boat_mesh = cpt.load_mesh("docs/examples/src/boat_200.mar", file_format="mar")
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
    dataset["inertia_matrix"] = body.inertia_matrix
    dataset["hydrostatic_stiffness"] = body.hydrostatic_stiffness
    rao = cpt.post_pro.rao(dataset, wave_direction=wave_direction)

    # COMPUTE FREE SURFACE ELEVATION
    # Compute the diffracted wave pattern
    incoming_waves_elevation = airy_waves_free_surface_elevation(fs.vertices, diffraction_result)
    diffraction_elevation = bem_solver.compute_free_surface_elevation(fs.vertices, diffraction_result)

    # Compute the wave pattern radiated by the RAO
    radiation_elevations_per_dof = {res.radiating_dof: bem_solver.compute_free_surface_elevation(fs.vertices, res) for res in radiation_results}
    radiation_elevation = sum(rao.sel(omega=omega, radiating_dof=dof).data * radiation_elevations_per_dof[dof] for dof in body.dofs)

    # SET UP ANIMATION
    # Compute the motion of each face of the mesh for the animation
    dofs_vertices_motions = {dof: body.dofs[dof].evaluate_motion_at_points(body.mesh.vertices) for dof in body.dofs}
    rao_vertices_motion = sum(rao.sel(omega=omega, radiating_dof=dof).data * dofs_vertices_motions[dof] for dof in body.dofs)

    # Set up scene
    animation = Animation(loop_duration=2*pi/omega)
    animation.add_body(body, vertices_motion=wave_amplitude*rao_vertices_motion)
    animation.add_free_surface(fs, wave_amplitude * (incoming_waves_elevation + diffraction_elevation + radiation_elevation))
    return animation


if __name__ == '__main__':
    body = generate_boat()
    fs = cpt.mesh_rectangle(size=(200.0, 200.0), resolution=(100, 100))

    anim = setup_animation(body, fs, omega=1.5, wave_amplitude=0.5, wave_direction=pi)
    # anim.run(camera_position=(70, 70, 100), resolution=(800, 600))
    anim.save("animated_boat.mp4", camera={"pos": (70, 70, 100)}, resolution=(800, 600))
