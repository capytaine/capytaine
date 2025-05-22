import numpy as np


# DOF string to index (WAMIT-style)
DOF_INDEX = {
    "Surge": 1,
    "Sway": 2,
    "Heave": 3,
    "Roll": 4,
    "Pitch": 5,
    "Yaw": 6,
}

# DOF string to type
DOF_TYPE = {
    "Surge": "trans", "Sway": "trans", "Heave": "trans",
    "Roll": "rot", "Pitch": "rot", "Yaw": "rot",
}

# Type-pair to exponent k
K_LOOKUP = {
    ("trans", "trans"): 3,
    ("trans", "rot"): 4,
    ("rot", "trans"): 4,
    ("rot", "rot"): 5,
}

def get_dof_index_and_k(dof_i, dof_j):
    i = DOF_INDEX[dof_i]
    j = DOF_INDEX[dof_j]
    k = K_LOOKUP[(DOF_TYPE[dof_i], DOF_TYPE[dof_j])]
    return i, j, k

def export_wamit_1_from_dataset(dataset, filename, L=1.0):
    """
    Export added mass and radiation damping coefficients to a WAMIT .1 file.

    Coefficients are normalized as:
        Aij = Aij / (rho * L^k)
        Bij = Bij / (omega * rho * L^k)

    Format:
        PER     I     J     Aij     Bij

    Parameters
    ----------
    dataset : xarray.Dataset
        Must contain 'added_mass' and 'radiation_damping'.
    filename : str
        Path to output .1 file.
    L : float
        Reference length scale for normalization (default is 1.0).
    """
    if "added_mass" not in dataset or "radiation_damping" not in dataset:
        raise ValueError("Missing 'added_mass' or 'radiation_damping' in dataset")

    rho = dataset["rho"].item()
    omegas = dataset["omega"].values
    periods = dataset["period"].values
    added_mass = dataset["added_mass"]
    damping = dataset["radiation_damping"]
    dofs = list(added_mass.coords["influenced_dof"].values)

    with open(filename, "w") as f:
        for omega, period in zip(omegas, periods):
            for dof_i in dofs:
                for dof_j in dofs:
                    A = added_mass.sel(omega=omega, influenced_dof=dof_i, radiating_dof=dof_j).item()
                    B = damping.sel(omega=omega, influenced_dof=dof_i, radiating_dof=dof_j).item()
                    i, j, k = get_dof_index_and_k(dof_i, dof_j)
                    norm = rho * (L ** k)
                    A_norm = A / norm
                    B_norm = B / (omega * norm)

                    f.write(f"{period:12.6e}\t{i:5d}\t{j:5d}\t{A_norm:12.6e}\t{B_norm:12.6e}\n")

def _format_excitation_line(period, beta_deg, i_dof, force):
    force_conj = np.conj(force)
    mod_f = np.abs(force_conj)
    phi_f = np.degrees(np.angle(force_conj))
    return "{:12.6e}\t{:12.6f}\t{:5d}\t{:12.6e}\t{:12.3f}\t{:12.6e}\t{:12.6e}\n".format(
        period, beta_deg, i_dof, mod_f, phi_f, force_conj.real, force_conj.imag
    )

def _export_wamit_excitation_force(dataset, field_name, filename):
    """
    Generic exporter for excitation-like forces in WAMIT .3-style format.

    Parameters
    ----------
    dataset : xarray.Dataset
        Dataset containing the desired complex-valued force field.
    field_name : str
        Name of the variable in dataset (e.g. "excitation_force", "Froude_Krylov_force").
    filename : str
        Output path for the .3/.3fk/.3sc file.
    """
    if field_name not in dataset:
        raise ValueError(f"Missing field '{field_name}' in dataset.")

    field = dataset[field_name]
    periods = dataset["period"].values
    omegas = dataset["omega"].values
    betas = field.coords["wave_direction"].values
    dofs = list(field.coords["influenced_dof"].values)

    with open(filename, "w") as f:
        for period, omega in zip(periods, omegas):
            for beta in betas:
                beta_deg = np.degrees(beta)
                for dof in dofs:
                    i_dof = DOF_INDEX.get(dof)
                    if i_dof is None:
                        raise KeyError(f"DOF '{dof}' is not recognized in DOF_INDEX mapping.")
                    force = field.sel(
                        omega=omega, wave_direction=beta, influenced_dof=dof
                    ).item()
                    line = _format_excitation_line(period, beta_deg, i_dof, force)
                    f.write(line)

def export_wamit_3_from_dataset(dataset, filename):
    """Export total excitation to WAMIT .3 file."""
    _export_wamit_excitation_force(dataset, "excitation_force", filename)

def export_wamit_3fk_from_dataset(dataset, filename):
    """Export Froude-Krylov contribution to WAMIT .3fk file."""
    _export_wamit_excitation_force(dataset, "Froude_Krylov_force", filename)

def export_wamit_3sc_from_dataset(dataset, filename):
    """Export scattered (diffraction) contribution to WAMIT .3sc file."""
    _export_wamit_excitation_force(dataset, "diffraction_force", filename)

def export_to_wamit(dataset, problem_name="WAMIT_OUTPUT", exports=("1", "3", "3fk", "3sc")):
    """
    Master function to export a Capytaine dataset to WAMIT-format files.

    Parameters
    ----------
    dataset : xarray.Dataset
        The Capytaine dataset containing hydrodynamic results.
    problem_name : str
        Base filename for WAMIT files (e.g. "output" → output.1, output.3fk, etc.).
    exports : tuple of str
        Which files to export: any combination of "1", "2", "3", "3fk", "3sc".
    """
    export_map = {
        "1":    ("radiation coefficients", export_wamit_1_from_dataset, ".1"),
        "3":    ("total excitation force", export_wamit_3_from_dataset, ".3"),
        "3fk":  ("Froude-Krylov force",    export_wamit_3fk_from_dataset, ".3fk"),
        "3sc":  ("diffraction force",      export_wamit_3sc_from_dataset, ".3sc"),
    }

    for key in exports:
        if key not in export_map:
            print(f"[!] Unknown export option '{key}' — skipping.")
            continue

        description, func, ext = export_map[key]
        filepath = f"{problem_name}{ext}"

        try:
            func(dataset, filepath)
            print(f"[✓] Exported {filepath} ({description})")
        except Exception as e:
            print(f"[X] Failed to export {filepath}: {e}")
