import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import capytaine as cpt
from capytaine.post_pro.mean_drift_force import near_field_mean_drift_force

mesh = cpt.mesh_parallelepiped(resolution=(30, 30, 30)).immersed_part()
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
rao = cpt.post_pro.rao(dataset) * 0.0
dataset.update(near_field_mean_drift_force(rao, results, solver, output_pressure=True))

ds = dataset.isel(wavelength=0, wave_direction_k=0, wave_direction_l=0, wave_direction=0)

amplitude = 0.1

first_order_hull_pressure = (
        ds["diffraction_pressure"] + ds["Froude_Krylov_pressure"]
        + (ds["radiation_pressure"] * rao.isel(wavelength=0, wave_direction=0)).sum("radiating_dof")
        ).values

first_order_field_real = amplitude*first_order_hull_pressure.real
first_order_field_imag = amplitude*first_order_hull_pressure.imag
omega = float(ds["omega"].values)
period = 2*np.pi/omega

second_order_field = amplitude**2*ds["second_order_pressure"].values
assert np.max(np.abs(second_order_field.imag)) < 1e-13
second_order_field = second_order_field.real

waterline_pressure = amplitude*ds["waterline_pressure"].values  # complex, oscillates in time
waterline_elevation = amplitude*ds["waterline_relative_elevation"].values  # complex, oscillates in time

cmap = plt.get_cmap('coolwarm')
vabs = max(np.abs(amplitude*first_order_hull_pressure).max(),
           np.abs(second_order_field).max(),
           np.abs(waterline_pressure).max())
norm = plt.Normalize(vmin=-vabs, vmax=vabs)
mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

fig = plt.figure(figsize=(15, 5))
gs = fig.add_gridspec(1, 3)

ax1 = fig.add_subplot(gs[0, 0], projection='3d')
body.mesh.show(
        backend="matplotlib",
        ax=ax1,
        facecolors=cmap(norm(first_order_field_real)),
        linewidth=0.1,
        )
hull_collection = ax1.collections[-1]
ax1.set_title("First order pressure")

ax4 = fig.add_subplot(gs[0, 1], projection='3d')
body.mesh.show(
    backend="matplotlib",
    ax=ax4,
    facecolors="none",
    edgecolor="k",
    linewidth=0.1,
)

ax3 = fig.add_subplot(gs[0, 2], projection='3d')
body.mesh.show(
        backend="matplotlib",
        ax=ax3,
        facecolors=cmap(norm(second_order_field)),
        linewidth=0.1,
        )
ax3.set_title("Mean drift pressure")

p0 = ds["waterline_mesh"].values[:, 0, :]
p1 = ds["waterline_mesh"].values[:, 1, :]


def waterline_quads(elevation):
    return [
            [(p0[i, 0], p0[i, 1], 0), (p1[i, 0], p1[i, 1], 0),
             (p1[i, 0], p1[i, 1], elevation[i]), (p0[i, 0], p0[i, 1], elevation[i])]
            for i in range(len(p0))
            ]


curtain = Poly3DCollection(waterline_quads(waterline_elevation.real), facecolors=cmap(norm(waterline_pressure.real)))
ax4.add_collection3d(curtain)
ax4.set_title("Waterline")
fig.suptitle(f"Amplitude = {amplitude}m")

fig.colorbar(mappable, ax=[ax1, ax3, ax4], shrink=0.6)


def update(frame):
    t = frame/n_frames * period
    phase = np.exp(-1j*omega*t)

    hull_pressure_t = (amplitude*first_order_hull_pressure * phase).real
    hull_collection.set_facecolor(cmap(norm(hull_pressure_t)))

    waterline_pressure_t = np.abs((waterline_pressure * phase).real)
    waterline_elevation_t = (waterline_elevation * phase).real
    curtain.set_verts(waterline_quads(waterline_elevation_t))
    curtain.set_facecolor(cmap(norm(waterline_pressure_t)))

    return hull_collection, curtain


n_frames = 30
animation_duration = 3.0  # seconds, wall-clock time for one wave period
anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=1000*animation_duration/n_frames, blit=False)

anim.save("second_order_pressure_field.mp4", fps=n_frames/animation_duration)

plt.show()
