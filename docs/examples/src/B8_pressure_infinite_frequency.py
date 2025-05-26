import numpy as np
import capytaine as cpt

# Initialize an example mesh and test case.
mesh = cpt.mesh_sphere().immersed_part()
body = cpt.FloatingBody(mesh, dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, 0)))

pb = cpt.RadiationProblem(body=body, omega=np.inf, radiating_dof="Heave")
# Zero and infinity are accepted values for "omega", as well as for "period", "wavelength" and "wavenumber"

solver = cpt.BEMSolver()
res = solver.solve(pb)

print(f"Heave-heave added mass: {res.added_mass['Heave']:.2e}")

p = np.real(res.pressure)
# p is a custom Python object storing the pressure on each panel

i_panel = 0  # Some panel for testing

# All pressures are actually infinite
print(f"Pressure on panel {i_panel}:", float(p[i_panel]))

# However p/omega^2 is finite (hence the finite added mass)
p_over_omega2 = p/pb.omega**2
print(f"Pressure/omega^2 on panel {i_panel}: {p_over_omega2[i_panel]:.2e}")
