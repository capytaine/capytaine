"""Example of extracting mean drift force pressure field from mean drift computation in Capytaine.
The whole file is heavily using DataArray, such that the meaning of each array dimension are explicit, for teaching purpose.
"""
import numpy as np
import xarray as xr
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

wavelengths = [2.0]
wave_directions = [0.0, np.pi/3]

solver = cpt.BEMSolver()
pbs = [cpt.DiffractionProblem(body=body, wavelength=wavelength, wave_direction=beta)
       for beta in wave_directions for wavelength in wavelengths]
pbs += [cpt.RadiationProblem(body=body, wavelength=wavelength, radiating_dof=dof)
        for dof in body.dofs for wavelength in wavelengths]
results = solver.solve_all(pbs)
ds = cpt.assemble_dataset(results)
rao = cpt.post_pro.rao(ds)
ds.update(near_field_mean_drift_force(rao, results, solver, output_pressure=True))

# The sea is made of the two waves, with the same frequency.
# Their complex amplitudes set their heights and their phases (at the origin, at t=0).
amplitude = xr.DataArray(
    np.array([[0.1, 0.05*np.exp(1j*np.pi/3)]]),
    dims=["wavelength", "wave_direction"],
    coords={"wavelength": ds.coords["wavelength"], "wave_direction": ds.coords["wave_direction"]}
)

# Pressures on the immersed hull, transferred to the faces of the full mesh
# without interpolation (nearest neighbors), setting zero on the emerged part.
distance, immersed_face = cKDTree(mesh.faces_centers).query(full_mesh.faces_centers)
immersed_face = xr.DataArray(immersed_face, dims=["full_mesh_face"])
is_immersed = xr.DataArray(distance < 1e-8, dims=["full_mesh_face"])
def remap_from_mesh_to_full_mesh(data):
    return data.isel(hull_face=immersed_face).where(is_immersed, 0.0)


# ============================ FIRST ORDER PRESSURE ============================
first_order_pressure = (
    amplitude * (
        ds["pressure_diffraction"] + ds["pressure_Froude_Krylov"]
        + (ds["pressure_radiation"] * rao).sum("radiating_dof")
        )).sum(dim=["wave_direction"])
# dims = (wavelength, hull_face)
first_order_pressure_full_mesh = remap_from_mesh_to_full_mesh(first_order_pressure)
# dims = (wavelength, full_mesh_face)


# ========================= SECOND ORDER HULL PRESSURE =========================
# Sum over all pairs of waves of the difference-frequency QTF
amplitude_kl = (amplitude.rename({"wave_direction": "wave_direction_k"})
                    * np.conjugate(amplitude.rename({"wave_direction": "wave_direction_l"})))
mean_drift_hull_pressure_per_pair_of_waves = amplitude_kl * ds["second_order_pressure"]
mean_drift_hull_pressure = mean_drift_hull_pressure_per_pair_of_waves.sum(dim=["wave_direction_k", "wave_direction_l"])
# dims = (wavelength, hull_face)
assert np.max(np.abs(mean_drift_hull_pressure.values.imag)) < 1e-10
mean_drift_hull_pressure = mean_drift_hull_pressure.real  # Time independent pressure
mean_drift_hull_pressure_full_mesh = remap_from_mesh_to_full_mesh(mean_drift_hull_pressure)
# dims = (wavelength, full_mesh_face)


# ====================== SECOND ORDER WATERLINE PRESSURE =======================
def mean_of_product(a, b):
    # Difference-frequency QTF of the product of the first order fields a and b,
    # 1/4 (a_k conj(b_l) + conj(a_l) b_k), summed over all pairs of waves (k, l).
    a_k = a.rename({"wave_direction": "wave_direction_k"})
    a_l = a.rename({"wave_direction": "wave_direction_l"})
    b_k = b.rename({"wave_direction": "wave_direction_k"})
    b_l = b.rename({"wave_direction": "wave_direction_l"})
    qtf = (a_k * np.conjugate(b_l) + np.conjugate(a_l) * b_k) / 4
    return qtf.sum(["wave_direction_k", "wave_direction_l"]).real  # The imaginary part is zero.


def remap_waterline_pressure(amplitude, waterline_pressure, waterline_relative_elevation, waterline_mesh, box_mesh):
    """Mean waterline pressure of the sea made of all the waves, spread on the faces of a box around the mean waterline.

    The waterline term of the mean drift force is a lineic force <p h> per unit length of waterline, where p is the
    pressure of the waterline strip and h the relative elevation of the free surface. It is spread in the vertical
    coordinate z around the mean waterline along a semi-elliptic profile of half-width |h| (the amplitude of h):

        <p>(z) = 2 <p h> / (pi |h|) * sqrt(1 - (z/|h|)^2)    for |z| < |h|,

    whose integral over z is <p h>. Each face of the box gets the value of this profile at its center, using the
    waterline segment with the same horizontal position.

    The waves have the same frequency: the mean value in the sea made of all of them is the sum over all pairs of
    waves (k, l) of the difference-frequency QTF, which includes the cross terms between different waves.

    Parameters
    ----------
    amplitude: DataArray with dims (..., wave_direction), complex
        Complex amplitude of each wave, possibly with other dimensions such as `wavelength`.
    waterline_pressure, waterline_relative_elevation: DataArrays with dims (..., wave_direction, waterline_segment), complex
        Outputs of `near_field_mean_drift_force(..., output_pressure=True)`, for a unit amplitude.
    waterline_mesh: DataArray with dims (waterline_segment, 2, 3)
        Output of `near_field_mean_drift_force(..., output_pressure=True)`.
    box_mesh: Mesh
        Mesh of the whole box, including its part above the waterline, with the same horizontal resolution as the
        immersed mesh.

    Returns
    -------
    DataArray with dims (..., full_mesh_face), real
        The other dimensions (such as `wavelength`) of the inputs are kept.
    """
    p = amplitude * waterline_pressure
    h = amplitude * waterline_relative_elevation
    ph = mean_of_product(p, h)  # dims = (..., waterline_segment)
    abs_h = np.sqrt(2*mean_of_product(h, h))  # dims = (..., waterline_segment), since <h h> = |h|^2/2

    _, segment = cKDTree(waterline_mesh.mean(axis=1)[:, :2]).query(box_mesh.faces_centers[:, :2])
    segment = xr.DataArray(segment, dims=["full_mesh_face"])
    u = xr.DataArray(box_mesh.faces_centers[:, 2], dims=["full_mesh_face"]) / abs_h.isel(waterline_segment=segment)

    # Value at the center of the faces. When |h| is not much larger than the faces, their average
    # (with the antiderivative u sqrt(1 - u^2) + arcsin(u) of 2 sqrt(1 - u^2)) is more accurate.
    return (2 * ph.isel(waterline_segment=segment) / (np.pi * abs_h.isel(waterline_segment=segment))
              * np.sqrt(np.clip(1 - u**2, 0.0, None)))  # zero for |u| > 1


mean_drift_waterline_pressure_full_mesh = remap_waterline_pressure(
        amplitude,
        ds["waterline_pressure"],
        ds["waterline_relative_elevation"],
        ds["waterline_mesh"],
        full_mesh,
        )
# dims = (wavelength, full_mesh_face)


# ================================== DISPLAY ==================================
# Display only a single frequency
freq_display = wavelengths[0]

first_order_pressure_displayed = first_order_pressure_full_mesh.sel(wavelength=freq_display)
mean_drift_hull_pressure_displayed = mean_drift_hull_pressure_full_mesh.sel(wavelength=freq_display)
mean_drift_waterline_pressure_displayed = mean_drift_waterline_pressure_full_mesh.sel(wavelength=freq_display)

omega = float(ds["omega"].sel(wavelength=freq_display))
period = 2*np.pi/omega

cmap = plt.get_cmap('coolwarm')
first_order_vabs = np.abs(first_order_pressure_displayed).max()
first_order_norm = plt.Normalize(vmin=-first_order_vabs, vmax=first_order_vabs)
mean_vabs = max(np.abs(mean_drift_hull_pressure_displayed).max(), np.abs(mean_drift_waterline_pressure_displayed).max())
mean_norm = plt.Normalize(vmin=-mean_vabs, vmax=mean_vabs)

fields = [
        ("First order pressure", first_order_pressure_displayed.real, first_order_norm),
        ("Mean pressure on hull", mean_drift_hull_pressure_displayed, mean_norm),
        ("Remapped mean waterline pressure", mean_drift_waterline_pressure_displayed, mean_norm),
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
                       for beta, a in zip(wave_directions, amplitude.sel(wavelength=freq_display).values)))
fig.colorbar(plt.cm.ScalarMappable(norm=first_order_norm, cmap=cmap), ax=axs[0], shrink=0.6)
fig.colorbar(plt.cm.ScalarMappable(norm=mean_norm, cmap=cmap), ax=axs[1:], shrink=0.6)

n_frames = 30
animation_duration = 3.0  # seconds, wall-clock time for one wave period


def update(frame):
    t = frame/n_frames * period
    phase = np.exp(-1j*omega*t)
    first_order_collection.set_facecolor(cmap(first_order_norm((first_order_pressure_displayed * phase).real)))
    return first_order_collection,


anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=1000*animation_duration/n_frames, blit=False)

anim.save("second_order_pressure_field.mp4", fps=n_frames/animation_duration)

plt.show()
