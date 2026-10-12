import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.spatial import cKDTree

import capytaine as cpt
from capytaine.post_pro.mean_drift_force import near_field_mean_drift_force

full_mesh = cpt.mesh_parallelepiped(resolution=(30, 30, 30))
mesh = full_mesh.immersed_part()
body = cpt.FloatingBody(
    mesh=mesh,
    lid_mesh=mesh.generate_lid(),
    dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)),
    center_of_mass=(0,0,0)
)

wavelength = 2.0

solver = cpt.BEMSolver()
pbs = [cpt.DiffractionProblem(body=body, wavelength=wavelength, wave_direction=0.0)]
pbs += [cpt.RadiationProblem(body=body, wavelength=wavelength, radiating_dof=dof) for dof in body.dofs]
results = solver.solve_all(pbs)
dataset = cpt.assemble_dataset(results)
rao = cpt.post_pro.rao(dataset)
dataset.update(near_field_mean_drift_force(rao, results, solver, output_pressure=True))

ds = dataset.isel(wavelength=0, wave_direction_k=0, wave_direction_l=0, wave_direction=0)

amplitude = 0.1

first_order_hull_pressure = amplitude*(
        ds["diffraction_pressure"] + ds["Froude_Krylov_pressure"]
        + (ds["radiation_pressure"] * rao.isel(wavelength=0, wave_direction=0)).sum("radiating_dof")
        ).values
omega = float(ds["omega"].values)
period = 2*np.pi/omega

second_order_field = amplitude**2*ds["second_order_pressure"].values
assert np.max(np.abs(second_order_field.imag)) < 1e-10
second_order_field = second_order_field.real

waterline_pressure = amplitude*ds["waterline_pressure"].values  # complex
waterline_elevation = amplitude*ds["waterline_relative_elevation"].values  # complex
waterline_middles = ds["waterline_mesh"].values.mean(axis=1)  # (segment, xyz)

# The time average of the waterline pressure has a semi-elliptic profile
#   <p*>(s) = Re(p conj(h)) / (pi |h|) * sqrt(1 - (s/|h|)^2)    for |s| < |h|
# in the vertical coordinate s around the mean waterline.
# It is remapped on the faces of the full mesh (including its emerged part)
# by averaging it over the vertical extent of each face, using the waterline
# segment closest to the face. The antiderivative of the profile is
#   Re(p conj(h)) / (2 pi) * G(s/|h|)    with G(u) = u sqrt(1 - u^2) + arcsin(u)  for |u| <= 1
def G(u):
    u = np.clip(u, -1.0, 1.0)
    return u*np.sqrt(1 - u**2) + np.arcsin(u)

_, closest_segment = cKDTree(waterline_middles[:, :2]).query(full_mesh.faces_centers[:, :2])
re_ph = (waterline_pressure * np.conjugate(waterline_elevation)).real[closest_segment]
abs_h = np.abs(waterline_elevation)[closest_segment]

faces_z = full_mesh.vertices[full_mesh.faces, 2]
z_min, z_max = faces_z.min(axis=1), faces_z.max(axis=1)
vertical_extent = z_max - z_min
assert np.all(vertical_extent[(z_min < abs_h) & (-abs_h < z_max)] > 0), "Horizontal faces crossing the waterline zone are not supported"
remapped_waterline_pressure = np.zeros(full_mesh.nb_faces)
nonzero = vertical_extent > 0
remapped_waterline_pressure[nonzero] = (
        re_ph[nonzero] / (2*np.pi)
        * (G(z_max[nonzero]/abs_h[nonzero]) - G(z_min[nonzero]/abs_h[nonzero]))
        / vertical_extent[nonzero]
        )

# Pressures on the immersed hull, transferred to the faces of the full mesh (zero on the emerged part)
distance, immersed_face = cKDTree(mesh.faces_centers).query(full_mesh.faces_centers)
is_immersed = distance < 1e-8
first_order_pressure = np.where(is_immersed, first_order_hull_pressure[immersed_face], 0.0)
hull_pressure = np.where(is_immersed, second_order_field[immersed_face], 0.0)

fields = {
        "First order pressure": first_order_pressure.real,
        "Mean drift pressure on hull": hull_pressure,
        "Remapped mean waterline pressure": remapped_waterline_pressure,
        }

cmap = plt.get_cmap('coolwarm')
vabs = max(np.abs(first_order_pressure).max(),
           np.abs(hull_pressure).max(),
           np.abs(remapped_waterline_pressure).max())
norm = plt.Normalize(vmin=-vabs, vmax=vabs)

fig = plt.figure(figsize=(15, 5))
axs = [fig.add_subplot(1, len(fields), i+1, projection='3d') for i in range(len(fields))]
for ax, (title, field) in zip(axs, fields.items()):
    full_mesh.show(
            backend="matplotlib",
            ax=ax,
            facecolors=cmap(norm(field)),
            linewidth=0.1,
            )
    ax.set_title(title)
    ax.view_init(elev=30, azim=-150)  # Show the face x = -0.5, on which the incoming waves hit
first_order_collection = axs[0].collections[-1]
fig.suptitle(f"Amplitude = {amplitude}m")
fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs, shrink=0.6)


def update(frame):
    t = frame/n_frames * period
    phase = np.exp(-1j*omega*t)
    first_order_collection.set_facecolor(cmap(norm((first_order_pressure * phase).real)))
    return first_order_collection,


n_frames = 30
animation_duration = 3.0  # seconds, wall-clock time for one wave period
anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=1000*animation_duration/n_frames, blit=False)

anim.save("waterline_pressure_remapping.mp4", fps=n_frames/animation_duration)

plt.show()
