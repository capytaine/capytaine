import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.spatial import cKDTree

import capytaine as cpt
from capytaine.post_pro.mean_drift_force import near_field_mean_drift_force
from capytaine.tools.perturbation_expansion import DifferenceFrequencyExpansion

full_mesh = cpt.mesh_parallelepiped(resolution=(30, 30, 30))
mesh = full_mesh.immersed_part()
body = cpt.FloatingBody(
    mesh=mesh,
    lid_mesh=mesh.generate_lid(),
    dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
    center_of_mass=(0,0,0)
)

wavelength = 2.0
wave_directions = [0.0, np.pi/3]

solver = cpt.BEMSolver()
pbs = [cpt.DiffractionProblem(body=body, wavelength=wavelength, wave_direction=beta) for beta in wave_directions]
pbs += [cpt.RadiationProblem(body=body, wavelength=wavelength, radiating_dof=dof) for dof in body.dofs]
results = solver.solve_all(pbs)
dataset = cpt.assemble_dataset(results)
rao = cpt.post_pro.rao(dataset)
dataset.update(near_field_mean_drift_force(rao, results, solver, output_pressure=True))

ds = dataset.isel(wavelength=0)
rao = rao.isel(wavelength=0)
omega = float(ds["omega"].values)
period = 2*np.pi/omega

# The sea is made of the two waves, with the same frequency.
# Their complex amplitudes set their heights and their phases (at the origin, at t=0).
amplitude = np.array([0.1, 0.05*np.exp(1j*np.pi/3)])


def remap_waterline_pressure(amplitude, waterline_pressure, waterline_relative_elevation, waterline_mesh, box_mesh):
    """Mean waterline pressure of the sea made of all the waves, spread on the faces of a box around the mean waterline.

    The waterline term of the mean drift force is a lineic force <p h> per unit length of waterline, where p is the
    pressure of the waterline strip and h the relative elevation of the free surface. It is spread in the vertical
    coordinate z around the mean waterline along a semi-elliptic profile of half-width |h| (the amplitude of h):

        <p>(z) = 2 <p h> / (pi |h|) * sqrt(1 - (z/|h|)^2)    for |z| < |h|,

    whose integral over z is <p h>. Each face of the box gets the value of this profile at its center, using the
    waterline segment with the same horizontal position.

    The waves have the same frequency: the mean value in the sea made of all of them is the sum over all pairs of
    waves (m, n) of the difference-frequency QTF, which includes the cross terms between different waves.

    Parameters
    ----------
    amplitude: complex or array of shape (nb_waves,), complex
        Complex amplitude of each wave.
    waterline_pressure, waterline_relative_elevation: arrays of shape (nb_waves, nb_waterline_segments), complex
        Outputs of `near_field_mean_drift_force(..., output_pressure=True)`, for a unit amplitude.
    waterline_mesh: array of shape (nb_waterline_segments, 2, 3)
        Output of `near_field_mean_drift_force(..., output_pressure=True)`.
    box_mesh: Mesh
        Mesh of the whole box, including its part above the waterline, with the same horizontal resolution as the
        immersed mesh.

    Returns
    -------
    array of shape (box_mesh.nb_faces,), real
    """
    amplitude = np.broadcast_to(amplitude, waterline_pressure.shape[:1])[:, None]
    p = DifferenceFrequencyExpansion([0.0, amplitude*waterline_pressure, 0.0])
    h = DifferenceFrequencyExpansion([0.0, amplitude*waterline_relative_elevation, 0.0])
    mean_ph = (p * h).coefficients[2].sum(axis=(0, 1)).real  # (segment,)
    abs_h = np.sqrt(2*(h * h).coefficients[2].sum(axis=(0, 1)).real)  # (segment,), since <h h> = |h|^2/2

    _, segment = cKDTree(waterline_mesh.mean(axis=1)[:, :2]).query(box_mesh.faces_centers[:, :2])
    u = box_mesh.faces_centers[:, 2] / abs_h[segment]
    # Value at the center of the faces. When |h| is not much larger than the faces, their average
    # (with the antiderivative u sqrt(1 - u^2) + arcsin(u) of 2 sqrt(1 - u^2)) is more accurate.
    return 2 * mean_ph[segment] / (np.pi * abs_h[segment]) * np.sqrt(np.clip(1 - u**2, 0.0, None))  # zero for |u| > 1


# First-order pressure on the hull of each wave (wave_direction, hull_face), and of the sea (hull_face,)
first_order_hull_pressure = (amplitude[:, None] * (
        ds["pressure_diffraction"] + ds["pressure_Froude_Krylov"]
        + (ds["pressure_radiation"] * rao).sum("radiating_dof")
        ).transpose("wave_direction", "hull_face").values).sum(axis=0)  # Sum on incoming waves

# Mean pressure on the hull: sum over all pairs of waves of the difference-frequency QTF
amplitude_kl = amplitude[:, None] * np.conjugate(amplitude[None, :])  # The QTF of the waves (k, l) scales as a_k conj(a_l)
second_order_pressure = amplitude_kl[:, :, None] * ds["second_order_pressure"].transpose("wave_direction_k", "wave_direction_l", "hull_face").values
assert np.max(np.abs(second_order_pressure.sum(axis=(0, 1)).imag)) < 1e-10
mean_hull_pressure = second_order_pressure.sum(axis=(0, 1)).real

mean_waterline_pressure = remap_waterline_pressure(
        amplitude,
        ds["waterline_pressure"].transpose("wave_direction", "waterline_segment").values,
        ds["waterline_relative_elevation"].transpose("wave_direction", "waterline_segment").values,
        ds["waterline_mesh"].values,
        full_mesh,
        )

# Pressures on the immersed hull, transferred to the faces of the full mesh (zero on the emerged part)
distance, immersed_face = cKDTree(mesh.faces_centers).query(full_mesh.faces_centers)
is_immersed = distance < 1e-8
first_order_pressure = np.where(is_immersed, first_order_hull_pressure[immersed_face], 0.0)
mean_pressure = np.where(is_immersed, mean_hull_pressure[immersed_face], 0.0)

cmap = plt.get_cmap('coolwarm')
first_order_vabs = np.abs(first_order_pressure).max()
first_order_norm = plt.Normalize(vmin=-first_order_vabs, vmax=first_order_vabs)
mean_vabs = max(np.abs(mean_pressure).max(), np.abs(mean_waterline_pressure).max())
mean_norm = plt.Normalize(vmin=-mean_vabs, vmax=mean_vabs)

fields = [
        ("First order pressure", first_order_pressure.real, first_order_norm),
        ("Mean pressure on hull", mean_pressure, mean_norm),
        ("Mean waterline pressure", mean_waterline_pressure, mean_norm),
        ]

fig = plt.figure(figsize=(15, 5))
axs = [fig.add_subplot(1, len(fields), i+1, projection='3d') for i in range(len(fields))]
for ax, (title, field, norm) in zip(axs, fields):
    full_mesh.show(
            backend="matplotlib",
            ax=ax,
            facecolors=cmap(norm(field)),
            linewidth=0.1,
            )
    ax.set_title(title)
    ax.view_init(elev=30, azim=-135)  # Show the faces x = -0.5 and y = -0.5, on which the incoming waves hit
first_order_collection = axs[0].collections[-1]
fig.suptitle(", ".join(f"wave from {np.degrees(beta):.0f}°: amplitude {np.abs(a):.2f}m, phase {np.degrees(np.angle(a)):.0f}°"
                       for beta, a in zip(wave_directions, amplitude)))
fig.colorbar(plt.cm.ScalarMappable(norm=first_order_norm, cmap=cmap), ax=axs[0], shrink=0.6)
fig.colorbar(plt.cm.ScalarMappable(norm=mean_norm, cmap=cmap), ax=axs[1:], shrink=0.6)


def update(frame):
    t = frame/n_frames * period
    phase = np.exp(-1j*omega*t)
    first_order_collection.set_facecolor(cmap(first_order_norm((first_order_pressure * phase).real)))
    return first_order_collection,


n_frames = 30
animation_duration = 3.0  # seconds, wall-clock time for one wave period
anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=1000*animation_duration/n_frames, blit=False)

anim.save("second_order_pressure_field.mp4", fps=n_frames/animation_duration)

plt.show()
