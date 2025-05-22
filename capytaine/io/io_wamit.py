import numpy as np

def _format_excitation_line(period, beta_deg, i_dof, force):
    mod_f = np.abs(force)
    phi_f = np.degrees(np.angle(force))
    return "{:12.6e}\t{:12.6f}\t{:5d}\t{:12.6e}\t{:12.3f}\t{:12.6e}\t{:12.6e}\n".format(
        period, beta_deg, i_dof, mod_f, phi_f, force.real, -force.imag
    )


def export_wamit_1_from_dataset(dataset, filename):
    """
    Export added mass and radiation damping coefficients to a WAMIT .1 file.

    Format:
        PER     I     J     Aij     Bij

    Where:
        - PER is the wave period (s),
        - I and J are DOF indices (1-based),
        - Aij is the added mass coefficient,
        - Bij is the radiation damping coefficient.

    Parameters
    ----------
    dataset : xarray.Dataset
        Must contain "added_mass" and "radiation_damping" terms.
    filename : str
        Output file path for the .1 file.
    """
    if "added_mass" not in dataset or "radiation_damping" not in dataset:
        raise ValueError("Missing 'added_mass' or 'radiation_damping' in dataset")

    added_mass = dataset["added_mass"]
    radiation_damping = dataset["radiation_damping"]
    omegas = added_mass.coords["omega"].values
    dofs = list(added_mass.coords["influenced_dof"].values)
 
    with open(filename, "w") as f:
        for omega in omegas:
            period = 2 * np.pi / omega
            for i, dof_i in enumerate(dofs, 1):
                for j, dof_j in enumerate(dofs, 1):
                    A = added_mass.sel(
                        influenced_dof=dof_i, radiating_dof=dof_j, omega=omega
                    ).item()
                    B = radiation_damping.sel(
                        influenced_dof=dof_i, radiating_dof=dof_j, omega=omega
                    ).item()
                    f.write(f"{period:12.6e}\t{i:5d}\t{j:5d}\t{A:12.6e}\t{B:12.6e}\n")
    
def export_wamit_2_from_dataset(dataset, filename):
    """
    Export excitation forces to a WAMIT .2 file (body-fixed reference frame).

    Format:
        PER     BETA     I     |Fi|     Pha(Fi)     Re(Fi)     Im(Fi)

    Where:
        - PER is the wave period (s),
        - BETA is the wave direction (deg),
        - I is the DOF index (1-based),
        - |Fi| is the excitation force magnitude,
        - Pha(Fi) is the phase in degrees,
        - Re(Fi), Im(Fi) are the real and imaginary parts.

    Parameters
    ----------
    dataset : xarray.Dataset
        Must contain the complex-valued "excitation" field.
    filename : str
        Output file path for the .2 file.
    """
    if "excitation_force" not in dataset:
        raise ValueError("Missing 'excitation_force' in dataset.")

    omega_list = dataset["omega"].values
    beta_list = dataset["wave_direction"].values
    excitation_force = dataset["excitation_force"]
    dofs = list(excitation_force.coords["influenced_dof"].values)
        
    with open(filename, "w") as f:
        for omega in omega_list:
            period = 2 * np.pi / omega
            for beta in beta_list:
                beta_deg = np.degrees(beta)
                for i_dof, dof in enumerate(dofs, 1):
                    force = excitation_force.sel(omega=omega, wave_direction=beta, influenced_dof=dof).item()
                    line = _format_excitation_line(period, beta_deg, i_dof, force)
                    f.write(line)

def export_wamit_3_from_dataset(dataset, filename):
    """
    Export diffraction forces to a WAMIT .3 file.

    Format:
        PER     BETA     I     |Fi|     Pha(Fi)     Re(Fi)     Im(Fi)

    Where:
        - PER is the wave period (s),
        - BETA is the wave direction (deg),
        - I is the DOF index (1-based),
        - |Fi| is the excitation force magnitude,
        - Pha(Fi) is the phase in degrees,
        - Re(Fi), Im(Fi) are the real and imaginary parts.

    Parameters
    ----------
    dataset : xarray.Dataset
        Must contain the complex-valued "diffraction_force" field in global frame.
    filename : str
        Output file path for the .3 file.
    """
    if "diffraction_force" not in dataset:
        raise ValueError("Missing diffraction_force in dataset.")

    diffraction_force = dataset["diffraction_force"]
    omegas = diffraction_force.coords["omega"].values
    betas = diffraction_force.coords["wave_direction"].values
    dofs = list(diffraction_force.coords["influenced_dof"].values)

    with open(filename, "w") as f:
        for omega in omegas:
            period = 2 * np.pi / omega
            for beta in betas:
                beta_deg = np.degrees(beta)
                for i_dof, dof in enumerate(dofs, 1):
                    force = dataset["diffraction_force"].sel(
                        omega=omega, wave_direction=beta, influenced_dof=dof
                    ).item()
                    line = _format_excitation_line(period, beta_deg, i_dof, force)
                    f.write(line)

def export_to_wamit(dataset, problem_name="WAMIT_OUTPUT", exports=("1", "2", "3")):
    """
    Master function to export a Capytaine dataset to WAMIT-format files.

    Parameters
    ----------
    dataset : xarray.Dataset
        The Capytaine dataset containing hydrodynamic results (radiation, excitation, etc.).
    problem_name : str
        Base name for the output files (e.g., 'problem' → 'problem.1', 'problem.2', etc.).
    exports : tuple of str
        Which files to export. Valid values are any combination of "1", "2", and "3":
            - "1": Added mass and radiation damping (WAMIT .1)
            - "2": Excitation from Haskind relations (WAMIT .2)
            - "3": Excitation from diffraction potential (WAMIT .3)
    """
    export_status = {
        "1": ("radiation coefficients", export_wamit_1_from_dataset, ".1"),
        "2": ("excitation forces (haskind)", export_wamit_2_from_dataset, ".2"),
        "3": ("excitation forces (diffraction)", export_wamit_3_from_dataset, ".3"),
    }

    for key in exports:
        if key not in export_status:
            print(f"[!] Unknown export option: '{key}' – skipping.")
            continue

        description, export_func, ext = export_status[key]
        output_file = f"{problem_name}{ext}"

        try:
            export_func(dataset, output_file)
            print(f"[✓] Exported {output_file} ({description})")
        except Exception as e:
            print(f"[X] Failed to export {output_file}: {e}")