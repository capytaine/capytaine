import numpy as np
import logging
from typing import Iterable, Tuple, TextIO
import xarray

logger = logging.getLogger(__name__)

DOF_INDEX = {"Surge": 1, "Sway": 2, "Heave": 3, "Roll": 4, "Pitch": 5, "Yaw": 6}

DOF_TYPE = {
    dof: "trans" if dof in {"Surge", "Sway", "Heave"} else "rot" for dof in DOF_INDEX
}
K_LOOKUP = {
    ("trans", "trans"): 3,
    ("trans", "rot"): 4,
    ("rot", "trans"): 4,
    ("rot", "rot"): 5,
}


def get_dof_index_and_k(dof_i: str, dof_j: str) -> Tuple[int, int, int]:
    """Get the degree of freedom indices and the corresponding stiffness matrix index.

    Parameters
    ----------
    dof_i: str
        The name of the first degree of freedom.
    dof_j: str
        The name of the second degree of freedom.

    Returns
    -------
    tuple
        A tuple containing the indices (i, j) and the stiffness matrix index k.
    """
    i = DOF_INDEX[dof_i]
    j = DOF_INDEX[dof_j]
    t_i = DOF_TYPE[dof_i]
    t_j = DOF_TYPE[dof_j]
    k = K_LOOKUP[(t_i, t_j)]
    return i, j, k


def export_wamit_hst(
    dataset: xarray.Dataset, filename: str, length_scale: float = 1.0
) -> None:
    """
    Export the nondimensional hydrostatic stiffness matrix to a WAMIT .hst file.

    Format:
        I     J     C(I,J)

    Parameters
    ----------
    dataset: xarray.Dataset
        Must contain 'hydrostatics' field with 'hydrostatic_stiffness' (named-dict or labeled 6x6 array).
        Must also contain 'rho' and 'g' (either in dataset or within hydrostatics).
    filename: str
        Output path for the .hst file.
    length_scale: float
        Reference length scale L for nondimensionalization.
    """
    if "hydrostatic_stiffness" not in dataset:
        raise ValueError("Dataset must contain a 'hydrostatic_stiffness' field.")

    # Reduce all extra dimensions to their first value, except the last two (should be 6x6)
    hydrostatic = dataset["hydrostatic_stiffness"]
    C = np.asarray(hydrostatic)
    if C is None or C.shape != (6, 6):
        raise ValueError("'hydrostatic_stiffness' must be a 6x6 matrix.")

    rho = float(np.atleast_1d(dataset.get("rho")).item())
    g = float(np.atleast_1d(dataset.get("g")).item())

    # DOF order used in Capytaine
    dof_names = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

    with open(filename, "w") as f:
        for i_local, dof_i in enumerate(dof_names):
            for j_local, dof_j in enumerate(dof_names):
                cij = C[i_local, j_local]
                i, j, k = get_dof_index_and_k(dof_i, dof_j)
                norm = rho * g * (length_scale**k)
                cij_nd = cij / norm
                f.write(f"{i:5d} {j:5d} {cij_nd:12.6e}\n")


def export_wamit_1(
    dataset: xarray.Dataset, filename: str, length_scale: float = 1.0
) -> None:
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
    dataset: xarray.Dataset
        Must contain 'added_mass', 'radiation_damping', 'omega', 'period'.
    filename: str
        Path to output .1 file.
    length_scale: float
        Reference length scale (L) for normalization.

    Raises
    ------
    ValueError
        If the data are missing in the dataset or if forward speed is not zero
    """
    if "added_mass" not in dataset or "radiation_damping" not in dataset:
        raise ValueError("Missing 'added_mass' or 'radiation_damping' in dataset")

    forward_speed = dataset["forward_speed"].values
    if not np.isclose(forward_speed, 0.0):
        raise ValueError("Forward speed must be zero for WAMIT export.")

    rho = dataset["rho"].item()
    omegas = np.asarray(dataset["omega"].values)
    periods = np.asarray(dataset["period"].values)
    added_mass = dataset["added_mass"]
    damping = dataset["radiation_damping"]
    dofs = list(added_mass.coords["influenced_dof"].values)

    omega_blocks = {"inf": [], "zero": [], "regular": []}

    # Identify regular frequencies (finite, > 0)
    regular_mask = np.isfinite(omegas) & (omegas > 0)
    regular_freqs = freqs[regular_mask]
    regular_omegas = omegas[regular_mask]
    regular_periods = periods[regular_mask]

    # Sort by increasing period
    sorted_indices = np.argsort(regular_periods)
    sorted_freqs = regular_freqs[sorted_indices]
    sorted_omegas = regular_omegas[sorted_indices]
    sorted_periods = regular_periods[sorted_indices]

    for omega, period, freq in zip(sorted_omegas, sorted_periods, sorted_freqs):
        for dof_i in dofs:
            for dof_j in dofs:
                A = added_mass.sel(
                    omega=omega, influenced_dof=dof_i, radiating_dof=dof_j
                ).item()
                B = damping.sel(
                    omega=omega, influenced_dof=dof_i, radiating_dof=dof_j
                ).item()
                j_dof, i_dof, k = get_dof_index_and_k(dof_i, dof_j)
                norm = rho * (length_scale**k)
                A_norm = A / norm
                B_norm = B / (omega * norm)
                line = f"{period:12.6e}\t{i_dof:5d}\t{j_dof:5d}\t{A_norm:12.6e}\t{B_norm:12.6e}\n"
                omega_blocks["regular"].append(line)

    with open(filename, "w") as f:
        f.writelines(omega_blocks["inf"])
        f.writelines(omega_blocks["zero"])
        f.writelines(omega_blocks["regular"])


def _format_excitation_line(
    period: float, beta_deg: float, i_dof: int, force: complex
) -> str:
    """Format a WAMIT excitation line.

    Parameters
    ----------
    period: float
        Wave period.
    beta_deg: float
        Wave direction (degrees).
    i_dof: int
        Degree of freedom index.
    force: complex
        Excitation force (complex).

    Returns
    -------
    str
        Formatted excitation line.
    """
    force_conj = np.conj(force)
    mod_f = np.abs(force_conj)
    phi_f = np.degrees(np.angle(force_conj))
    return "{:12.6e}\t{:12.6f}\t{:5d}\t{:12.6e}\t{:12.3f}\t{:12.6e}\t{:12.6e}\n".format(
        period, beta_deg, i_dof, mod_f, phi_f, force_conj.real, force_conj.imag
    )


def _write_wamit_excitation_line(
    f: TextIO,
    period: float,
    omega: float,
    beta: float,
    dof: str,
    field: xarray.DataArray,
    rho: float,
    g: float,
    wave_amplitude: float,
    length_scale: float,
):
    """Write a WAMIT excitation line to the file.

    Parameters
    ----------
    f: TextIO
        File object to write to.
    period: float
        Wave period.
    omega: float
        Wave frequency.
    beta: float
        Wave direction (radians).
    dof: str
        Degree of freedom.
    field: xarray.DataArray
        Field containing the excitation forces.
    rho: float
        Water density.
    g: float
        Gravitational acceleration.
    wave_amplitude: float
        Wave amplitude.
    length_scale: float
        Length scale for normalization.

    Raises
    ------
    KeyError
        If the degree of freedom is not recognized.
    """
    beta_deg = np.degrees(beta)
    i_dof = DOF_INDEX.get(dof)
    if i_dof is None:
        raise KeyError(f"DOF '{dof}' is not recognized in DOF_INDEX mapping.")
    force = field.sel(omega=omega, wave_direction=beta, influenced_dof=dof).item()
    dof_type = DOF_TYPE.get(dof, "trans")
    m = 2 if dof_type == "trans" else 3
    norm = rho * g * wave_amplitude * (length_scale**m)
    force_normalized = force / norm
    line = _format_excitation_line(period, beta_deg, i_dof, force_normalized)
    f.write(line)


def _export_wamit_excitation_force(
    dataset: xarray.Dataset,
    field_name: str,
    filename: str,
    length_scale: float = 1.0,
    wave_amplitude: float = 1.0,
):
    """
    Generic exporter for excitation-like forces in WAMIT .3-style format.

    Format:
        PER     BETA    I     Fmagnitude  Fphase  Freal  Fimaginary

    Parameters
    ----------
    dataset: xarray.Dataset
        Dataset containing the desired complex-valued force field.
    field_name: str
        Name of the variable in dataset (e.g. "excitation_force", "Froude_Krylov_force").
    filename: str
        Output path for the .3/.3fk/.3sc file.
    """
    forward_speed = dataset["forward_speed"].values
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

    sorted_indices = np.argsort(periods)
    periods = periods[sorted_indices]
    omegas = omegas[sorted_indices]

    with open(filename, "w") as f:
        for period, omega in zip(periods, omegas):
            # TODO: So far, we skip omega=0 and omega=inf for .3-style exports
            # TODO: Check https://www.wamit.com/manual6.4/Chap4.pdf - CH4-13
            if np.isclose(omega, 0.0) or np.isinf(omega):
                continue
            for beta in betas:
                for dof in dofs:
                    _write_wamit_excitation_line(
                        f,
                        period,
                        omega,
                        beta,
                        dof,
                        field,
                        rho,
                        g,
                        wave_amplitude,
                        length_scale,
                    )


def export_wamit_3(dataset: xarray.Dataset, filename: str) -> None:
    """Export total excitation to WAMIT .3 file.

    Parameters
    ----------
    dataset: xarray.Dataset
        Dataset containing the desired complex-valued force field.
    filename: str
        Output path for the .3 file.
    """
    _export_wamit_excitation_force(dataset, "excitation_force", filename)


def export_wamit_3fk(dataset: xarray.Dataset, filename: str) -> None:
    """Export Froude-Krylov contribution to WAMIT .3fk file.

    Parameters
    ----------
    dataset: xarray.Dataset
        Dataset containing the desired complex-valued force field.
    filename: str
        Output path for the .3fk file.
    """
    _export_wamit_excitation_force(dataset, "Froude_Krylov_force", filename)


def export_wamit_3sc(dataset: xarray.Dataset, filename: str) -> None:
    """Export scattered (diffraction) contribution to WAMIT .3sc file.

    Parameters
    ----------
    dataset: xarray.Dataset
        Dataset containing the desired complex-valued force field.
    filename: str
        Output path for the .3sc file.
    """
    _export_wamit_excitation_force(dataset, "diffraction_force", filename)


def export_to_wamit(
    dataset: xarray.Dataset,
    problem_name: str,
    exports: Iterable[str] = ("1", "3", "3fk", "3sc", "hst"),
) -> None:
    """
    Master function to export a Capytaine dataset to WAMIT-format files.

    Parameters
    ----------
    dataset: xarray.Dataset
        Dataset containing the desired complex-valued force field.
    problem_name: str
        Base filename for WAMIT files (e.g. "output" → output.1, output.3fk, etc.).
    exports: iterable of str
        Which files to export: any combination of "1", "3", "3fk", "3sc", "hst".
    """
    export_map = {
        "1": ("radiation coefficients", export_wamit_1, ".1"),
        "3": ("total excitation force", export_wamit_3, ".3"),
        "3fk": ("Froude-Krylov force", export_wamit_3fk, ".3fk"),
        "3sc": ("diffraction force", export_wamit_3sc, ".3sc"),
        "hst": ("hydrostatics", export_wamit_hst, ".hst"),
    }

    for key in exports:
        if key not in export_map:
            logger.warning(
                f"Export to WAMIT format: unknown option '{key}' — skipping."
            )
            continue

        description, func, ext = export_map[key]
        filepath = f"{problem_name}{ext}"

        try:
            func(dataset, filepath)
            logger.info(f"Export to WAMIT format: exported {filepath} ({description})")
        except Exception as e:
            logger.warning(f"Export to WAMIT format: did not export {filepath}: {e}")
