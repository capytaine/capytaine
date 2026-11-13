from numpy import pi

import capytaine as cpt
from capytaine.bem.airy_waves import airy_waves_free_surface_elevation
from capytaine.ui.vedo_animations import Animation

from vedo import Light
from matplotlib.colors import LinearSegmentedColormap

# Colors of the Capytaine logo
NAVY, GOLD = "#21355d", "#fbc02d"
# Lighter colormap for the waves, inspired by the logo colors, such that the
# dark boat stands out. The warm grey avoids going through green between blue and gold.
WAVES_CMAP = LinearSegmentedColormap.from_list("capytaine_waves", ["#5a7aa6", "#4dbde0", "#c4bfb6", "#f2cf6a"])


bem_solver = cpt.BEMSolver()


def generate_boat():
    boat_mesh = cpt.load_mesh("docs/examples/src/boat_200.mar", file_format="mar")
    boat = cpt.FloatingBody(
            mesh=boat_mesh,
            dofs=cpt.rigid_body_dofs(rotation_center=boat_mesh.center_of_buoyancy),
            center_of_mass = boat_mesh.center_of_buoyancy,
            )
    boat.inertia_matrix = boat.compute_rigid_body_inertia() / 40 # Artificially lower to have a more appealing animation
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
    animation.add_body(body, vertices_motion=wave_amplitude*rao_vertices_motion, color=NAVY, edge_color=GOLD)
    animation.add_free_surface(fs, wave_amplitude * (incoming_waves_elevation + diffraction_elevation + radiation_elevation),
                               cmap=WAVES_CMAP)
    return animation


if __name__ == '__main__':
    body = generate_boat()
    fs = cpt.mesh_rectangle(size=(200.0, 200.0), resolution=(100, 100))

    anim = setup_animation(body, fs, omega=1.5, wave_amplitude=1.0, wave_direction=pi)

    params = dict(
        camera={"pos": (70, 70, 100)},
        lights=[Light([0, 0, 100], intensity=0.8)],
        resolution=(960, 720),  # Multiple of 16 to avoid resizing by the video encoder
        axes=0,
    )
    # anim.run(**params)
    anim.save("animated_boat.mp4", **params)  # Used as docs/_static/front_page_animation.mp4
