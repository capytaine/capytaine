import numpy as np

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
    dof_order : list of str, optional
        Custom DOF ordering. If None, the order from the dataset is used.
    """
    if "added_mass" not in dataset or "radiation_damping" not in dataset:
        raise ValueError("Missing 'added_mass' or 'radiation_damping' in dataset")

    added_mass = dataset["added_mass"]
    radiation_damping = dataset["radiation_damping"]
    omegas = added_mass.coords["omega"].values
    dofs = list(added_mass.coords["influenced_dof"].values)
 
    with open(filename, "w") as f:
        for omega in omegas:
            for i, dof_i in enumerate(dofs, 1):
                for j, dof_j in enumerate(dofs, 1):
                
                    period = 2 * np.pi / omega
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
    dof_order : list of str, optional
        Custom DOF ordering. If None, the order from the dataset is used.
    """
    if "excitation" not in dataset:
        raise ValueError("Missing excitation in dataset.")

    omega_list = dataset["omega"].values
    beta_list = dataset["wave_direction"].values
    excitation = dataset["excitation"]
    dofs = list(excitation.coords["influenced_dof"].values)
        
    with open(filename, "w") as f:
        for omega in omega_list:
            period = 2 * np.pi / omega
            for beta in beta_list:
                beta_deg = np.degrees(beta)
                for i_dof, dof in enumerate(dofs, 1):
                    F = excitation.sel(omega=omega, wave_direction=beta, influenced_dof=dof).item()
                    absF = np.abs(F)
                    phaF = np.degrees(np.angle(F))
                    line = "{:12.6e}\t{:12.6f}\t{:5d}\t{:12.6e}\t{:12.3f}\t{:12.6e}\t{:12.6e}\n".format(period, beta_deg, i_dof, absF, phaF, F.real, -F.imag)
                    f.write(line)

def export_wamit_3_from_dataset(dataset, filename):
    """
    Export diffraction forces to a WAMIT .3 file (earth-fixed reference frame).

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
    dof_order : list of str, optional
        Custom DOF ordering. If None, the order from the dataset is used.
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
            for i_dof, dof in enumerate(dofs, 1):
                for beta in betas:
                    beta_deg = np.degrees(beta)
                    F = diffraction_force.sel(influenced_dof=dof, omega=omega, wave_direction=beta).item()
                    absF = np.abs(F)
                    phaF = np.degrees(np.angle(F))
                    f.write("{:12.6e}\t{:12.6f}\t{:5d}\t{:12.6e}\t{:12.3f}\t{:12.6e}\t{:12.6e}\n".format(
                        period, beta_deg, i_dof, absF, phaF, F.real, -F.imag))

def export_to_wamit(dataset, problem_name="WAMIT_OUTPUT", exports=("1", "2", "3")):
    """
    Master function to export Capytaine dataset to WAMIT format files.

    Parameters
    ----------
    dataset : xarray.Dataset
        The Capytaine dataset containing radiation and/or diffraction results.
    problem_name : str
        Base name used for output files (e.g., 'problem' → 'problem.1').
    exports : tuple of str
        Choose which files to export. Any combination of ("1", "2", "3").
    """

    if "excitation" not in dataset and "diffraction_force" in dataset:
        dataset["excitation"] = dataset["diffraction_force"]


    if "1" in exports:
        try:
            export_wamit_1_from_dataset(dataset, f"{problem_name}.1")
            print(f"[✓] Exported {problem_name}.1 (radiation coefficients)")
        except Exception as e:
            print(f"[!] Failed to export {problem_name}.1: {e}")

    if "2" in exports:
        try:
            export_wamit_2_from_dataset(dataset, f"{problem_name}.2")
            print(f"[✓] Exported {problem_name}.2 (excitation forces, mod/phase)")
        except Exception as e:
            print(f"[!] Failed to export {problem_name}.2: {e}")

    if "3" in exports:
        try:
            export_wamit_3_from_dataset(dataset, f"{problem_name}.3")
            print(f"[✓] Exported {problem_name}.3 (excitation forces, Re/Im)")
        except Exception as e:
            print(f"[!] Failed to export {problem_name}.3: {e}")
