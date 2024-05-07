"""
Example file based on Malenica 1995:
Š. Malenica, P.J. Clark, B. Molin,
Wave and current forces on a vertical cylinder free to surge and sway,
Applied Ocean Research, Volume 17, Issue 2, 1995, Pages 79-90,
https://doi.org/10.1016/0141-1187(95)00002-I.

This script is based on the work of the Maritime Technology Division (MTD) of Ghent University, Belgium.
Please refer to following papers:

I. Herdayanditya, L. Donatini, G. Verao Fernandez, A. B. K. Pribadi, and P. Rauwoens,
“Waves-current effect investigation on monopile excitation force employing approximate forward speed approach,”
in 6th MASHCON : international conference on ship manoeuvring in shallow and confined water
with special focus on port manoeuvres, Glasgow, UK, 2022, pp. 73–80.
http://hdl.handle.net/1854/LU-8756871

L. Donatini, I. Herdayanditya, G. Verao Fernandez, A. B. K. Pribadi, E. Lataire, and G. Delefortrie,
“Implementation of forward speed effects on an open source seakeeping solver,”
in 6th MASHCON : international conference on ship manoeuvring in shallow and confined water
with special focus on port manoeuvres, Glasgow, UK, 2022, pp. 20–33.
http://hdl.handle.net/1854/LU-8756868

"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import capytaine as cpt

g = 9.81
rho = 1025.0

water_depth = 1.0
radius = 1.0
center = (0, 0, -0.5)
resolution = (2, 20, 20)

mesh = cpt.mesh_vertical_cylinder(length=1.001*water_depth, radius=radius, center=center, resolution=resolution)
cylinder = cpt.FloatingBody(mesh=mesh, dofs=cpt.rigid_body_dofs(rotation_center=center), center_of_mass=center)
cylinder = cylinder.immersed_part(water_depth=water_depth)

wavenumber_range = np.linspace(0.2,2,25)
froude_number_range = np.array([0, 0.05, -0.05])
forward_speed_range = froude_number_range * np.sqrt(g*radius)

test_matrix = xr.Dataset(coords={
    'wavenumber': wavenumber_range,
    'forward_speed': forward_speed_range,
    'wave_direction': [0.0],
    'radiating_dof': ["Surge"],
    'water_depth': [water_depth],
    'rho': [rho],
    })
solver = cpt.BEMSolver()
dataset = solver.fill_dataset(test_matrix, cylinder)
dataset["excitation_force"] = dataset["Froude_Krylov_force"] + dataset["diffraction_force"]

dataset.coords["encounter_omega"] = dataset["omega"] - dataset["wavenumber"] * dataset["forward_speed"] * np.cos(dataset["wave_direction"])
dataset.coords["Froude_number"] = dataset.coords["forward_speed"]/np.sqrt(g*radius)

plt.figure()
dataset["norm_excitation"] = np.abs(dataset["excitation_force"]/(g*rho*radius**2))
dataset["norm_excitation"].attrs['long_name'] = "Normalized excitation force F/ρgr²"
dataset["norm_excitation"].sel(wave_direction=0.0, influenced_dof="Surge").plot(x="wavenumber", hue="Froude_number")

plt.figure()
dataset["norm_added_mass"] = dataset["added_mass"]/(rho*radius**3)
dataset["norm_added_mass"].attrs['long_name'] = "Normalized added mass A/ρr³"
dataset["norm_added_mass"].sel(influenced_dof="Surge", radiating_dof="Surge").plot(x="wavenumber", hue="Froude_number")

plt.figure()
dataset["norm_rad_damping"] = dataset["radiation_damping"]/(rho*dataset.encounter_omega*radius**3)
dataset["norm_rad_damping"].attrs['long_name'] = "Normalized radiation damping B/ρωr³"
dataset["norm_rad_damping"].sel(influenced_dof="Surge", radiating_dof="Surge").plot(x="wavenumber", hue="Froude_number")

print(dataset)
plt.show()
