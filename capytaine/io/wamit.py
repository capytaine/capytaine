import numpy as np
import logging

logger = logging.getLogger(__name__)

DOF_INDEX = {"Surge": 1, "Sway": 2, "Heave": 3, "Roll": 4, "Pitch": 5, "Yaw": 6}

DOF_TYPE = {dof: "trans" if dof in {"Surge", "Sway", "Heave"} else "rot" for dof in DOF_INDEX}
K_LOOKUP = {
    ("trans", "trans"): 3,
    ("trans", "rot"): 4,
    ("rot", "trans"): 4,
    ("rot", "rot"): 5,
}

def get_dof_index_and_k(dof_i, dof_j):
    i = DOF_INDEX[dof_i]
    j = DOF_INDEX[dof_j]
    t_i = DOF_TYPE[dof_i]
    t_j = DOF_TYPE[dof_j]
    k = K_LOOKUP[(t_i, t_j)]
    return i, j, k

def export_wamit_hst(dataset, filename, length_scale=1.0):
    """
    Export the nondimensional hydrostatic stiffness matrix to a WAMIT .hst file.

    Format:
        I     J     C(I,J)

    Parameters
    ----------
    dataset : xarray.Dataset
        Must contain 'hydrostatics' field with 'hydrostatic_stiffness' (named-dict or labeled 6x6 array).
        Must also contain 'rho' and 'g' (either in dataset or within hydrostatics).
    filename : str
        Output path for the .hst file.
    length_scale : float
        Reference length scale L for nondimensionalization.
    """
    if "hydrostatics" not in dataset:
        raise ValueError("Dataset must contain a 'hydrostatics' field.")

    hydro = dataset["hydrostatics"].item()
    C = np.asarray(hydro.get("hydrostatic_stiffness", None))
    if C is None or C.shape != (6, 6):
        raise ValueError("'hydrostatic_stiffness' must be a 6x6 matrix.")

    rho = dataset.get("rho", hydro.get("rho", 1025.0))
    g = dataset.get("g", hydro.get("g", 9.81))
    L = length_scale

    # DOF order used in Capytaine
    dof_names = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

    with open(filename, "w") as f:
        for i_local, dof_i in enumerate(dof_names):
            for j_local, dof_j in enumerate(dof_names):
                cij = C[i_local, j_local]
                if np.isclose(cij, 0.0):
                    continue
                i, j, k = get_dof_index_and_k(dof_i, dof_j)
                norm = rho * g * (L ** k)
                cij_nd = cij / norm
                f.write(f"{i:5d} {j:5d} {cij_nd:12.6e}\n")


def export_wamit_1(dataset, filename, length_scale=1.0):
    """
    Export added mass and radiation damping coefficients to a WAMIT .1 file.

    Coefficients are normalized as:
        Aij = Aij / (rho * length_scale^k)
        Bij = Bij / (omega * rho * length_scale^k)

    Format:
        PER     I     J     Aij [Bij]

    Special handling:
        - For PER = -1 (omega = inf), only Aij is written.
        - For PER = 0  (omega = 0), only Aij is written.

    Parameters
    ----------
    dataset : xarray.Dataset
        Must contain 'added_mass', 'radiation_damping', 'omega', 'period'.
    filename : str
        Path to output .1 file.
    length_scale : float
        Reference length scale (L) for normalization.
    """
    if "added_mass" not in dataset or "radiation_damping" not in dataset:
        raise ValueError("Missing 'added_mass' or 'radiation_damping' in dataset")

    forward_speed = dataset['forward_speed'].values
    if not np.isclose(forward_speed, 0.0):
        raise ValueError("Forward speed must be zero for WAMIT export.")

    rho = dataset["rho"].item()
    omegas = dataset["omega"].values
    periods = dataset["period"].values
    added_mass = dataset["added_mass"]
    damping = dataset["radiation_damping"]
    dofs = list(added_mass.coords["influenced_dof"].values)

    omega_blocks = {
        "inf": [],
        "zero": [],
        "regular": []
    }

    for omega, period in zip(omegas, periods):
        for dof_i in dofs:
            for dof_j in dofs:
                A = added_mass.sel(omega=omega, influenced_dof=dof_i, radiating_dof=dof_j).item()
                B = damping.sel(omega=omega, influenced_dof=dof_i, radiating_dof=dof_j).item()
                i_dof, j_dof, k = get_dof_index_and_k(dof_i, dof_j)
                norm = rho * (length_scale ** k)
                A_norm = A / norm

                if np.isinf(omega):
                    line = f"{-1.0:12.6e}\t{i_dof:5d}\t{j_dof:5d}\t{A_norm:12.6e}\n"
                    omega_blocks["inf"].append(line)
                elif np.isclose(omega, 0.0):
                    line = f"{0.0:12.6e}\t{i_dof:5d}\t{j_dof:5d}\t{A_norm:12.6e}\n"
                    omega_blocks["zero"].append(line)
                else:
                    B_norm = B / (omega * norm)
                    line = f"{period:12.6e}\t{i_dof:5d}\t{j_dof:5d}\t{A_norm:12.6e}\t{B_norm:12.6e}\n"
                    omega_blocks["regular"].append(line)

    with open(filename, "w") as f:
        f.writelines(omega_blocks["inf"])
        f.writelines(omega_blocks["zero"])
        f.writelines(omega_blocks["regular"])


def _format_excitation_line(period, beta_deg, i_dof, force):
    force_conj = np.conj(force)
    mod_f = np.abs(force_conj)
    phi_f = np.degrees(np.angle(force_conj))
    return "{:12.6e}\t{:12.6f}\t{:5d}\t{:12.6e}\t{:12.3f}\t{:12.6e}\t{:12.6e}\n".format(
        period, beta_deg, i_dof, mod_f, phi_f, force_conj.real, force_conj.imag
    )

def _write_wamit_excitation_line(f, period, omega, beta, dof, field, rho, g, wave_amplitude, length_scale):
    beta_deg = np.degrees(beta)
    i_dof = DOF_INDEX.get(dof)
    if i_dof is None:
        raise KeyError(f"DOF '{dof}' is not recognized in DOF_INDEX mapping.")
    force = field.sel(
        omega=omega, wave_direction=beta, influenced_dof=dof
    ).item()
    dof_type = DOF_TYPE.get(dof, "trans")
    m = 2 if dof_type == "trans" else 3
    norm = rho * g * wave_amplitude * (length_scale ** m)
    force_normalized = force / norm
    line = _format_excitation_line(period, beta_deg, i_dof, force_normalized)
    f.write(line)

def _export_wamit_excitation_force(dataset, field_name, filename, length_scale=1.0, wave_amplitude=1.0):
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
    forward_speed = dataset['forward_speed'].values
    if not np.isclose(forward_speed, 0.0):
        raise ValueError("Forward speed must be zero for WAMIT export.")
        
    if field_name not in dataset:
        raise ValueError(f"Missing field '{field_name}' in dataset.")

    field = dataset[field_name]
    periods = dataset["period"].values
    omegas = dataset["omega"].values
    rho = dataset["rho"].values
    g = dataset["g"].values

    betas = field.coords["wave_direction"].values
    dofs = list(field.coords["influenced_dof"].values)

    with open(filename, "w") as f:
        for period, omega in zip(periods, omegas):
            for beta in betas:
                for dof in dofs:
                    _write_wamit_excitation_line(
                        f, period, omega, beta, dof, field, rho, g, wave_amplitude, length_scale
                    )

def export_wamit_3(dataset, filename):
    """Export total excitation to WAMIT .3 file."""
    _export_wamit_excitation_force(dataset, "excitation_force", filename)

def export_wamit_3fk(dataset, filename):
    """Export Froude-Krylov contribution to WAMIT .3fk file."""
    _export_wamit_excitation_force(dataset, "Froude_Krylov_force", filename)

def export_wamit_3sc(dataset, filename):
    """Export scattered (diffraction) contribution to WAMIT .3sc file."""
    _export_wamit_excitation_force(dataset, "diffraction_force", filename)

def export_to_wamit(dataset, problem_name, exports=("1", "3", "3fk", "3sc", "hst")):
    """
    Master function to export a Capytaine dataset to WAMIT-format files.

    Parameters
    ----------
    dataset : xarray.Dataset
        The Capytaine dataset containing hydrodynamic results.
    problem_name : str
        Base filename for WAMIT files (e.g. "output" → output.1, output.3fk, etc.).
    exports : tuple of str
        Which files to export: any combination of "1", "3", "3fk", "3sc", "hst".
    """
    export_map = {
        "1":    ("radiation coefficients", export_wamit_1, ".1"),
        "3":    ("total excitation force", export_wamit_3, ".3"),
        "3fk":  ("Froude-Krylov force",    export_wamit_3fk, ".3fk"),
        "3sc":  ("diffraction force",      export_wamit_3sc, ".3sc"),
        "hst":  ("hydrostatics",           export_wamit_hst, ".hst"),
    }

    for key in exports:
        if key not in export_map:
            logger.warning(f"Unknown export option '{key}' — skipping.")
            continue

        description, func, ext = export_map[key]
        filepath = f"{problem_name}{ext}"

        try:
            func(dataset, filepath)
            logger.info(f"Exported {filepath} ({description})")
        except Exception as e:
            logger.error(f"Failed to export {filepath}: {e}")