import numpy as np
import logging
from typing import Union, Iterable, Tuple, TextIO
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


def check_dataset_ready_for_export(ds: xarray.Dataset) -> None:
    """
    Sanity checks to validate that the dataset is exportable to BEMIO/WAMIT-like formats.

    Parameters
    ----------
    ds : xarray.Dataset
        The dataset to be validated.

    Raises
    ------
    ValueError
        If any unsupported coordinate has multiple values or
        if non-rigid-body DOFs are present.
    """
    # 1. Check for singleton coordinates
    critical_coords = ["water_depth", "g", "rho"]
    coords_with_multiple_values = [
        k
        for k in critical_coords
        if k in ds.coords and len(ds.coords[k].dims) > 0 and ds.sizes[k] > 1
    ]

    if coords_with_multiple_values:
        msg = (
            "Export formats like WAMIT require only one value for each of the following coordinates: "
            f"{', '.join(critical_coords)}.\n"
            f"Problematic dimensions: {coords_with_multiple_values}.\n"
            "You can extract a subset using:\n"
            f"    ds_slice = ds.sel({', '.join([f'{k}={str(ds.coords[k].values[0])}' for k in coords_with_multiple_values])})"
        )
        raise ValueError(msg)

    # 2. Check for rigid-body DOFs only
    rigid_body_dofs = ("Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw")
    if "influenced_dof" in ds.coords:
        dofs = set(ds.influenced_dof.values)
        non_rigid_dofs = dofs.difference(set(rigid_body_dofs))
        if non_rigid_dofs:
            raise ValueError(
                "WAMIT Export is only supported for single rigid body.\n"
                f"Unexpected DOFs: {non_rigid_dofs}.\n"
                f"Allowed DOFs: {rigid_body_dofs}"
            )


def identify_frequency_axis(
    dataset: Union[xarray.Dataset, xarray.DataArray],
) -> Tuple[str, np.ndarray, np.ndarray]:
    """
    Identify the frequency axis in the dataset and return its values along with the periods.

    Parameters
    ----------
    dataset : xarray.Dataset or xarray.DataArray
        Dataset that must include 'period' coordinate and at least one of:
        'omega', 'freq', 'period', 'wavenumber', or 'wavelength' as dimension.

    Returns
    -------
    freq_key : str
        The name of the main frequency-like coordinate.
    freq_vals : np.ndarray
        The values from the frequency coordinate (as present in the dataset).
    period_vals : np.ndarray
        The values of the 'period' coordinate in seconds.

    Raises
    ------
    ValueError
        If 'period' is not a coordinate or if no frequency-like dimension is found.
    """
    allowed_keys = {"omega", "freq", "wavenumber", "wavelength", "period"}
    dataset_dims = set(dataset.dims)
    keys_in_dataset = dataset_dims & allowed_keys

    if "period" not in dataset.coords:
        raise ValueError("Dataset must contain 'period' as a coordinate.")

    # Prioritize 'period' if it is one of the dimensions
    if "period" in dataset_dims:
        freq_key = "period"
    elif len(keys_in_dataset) >= 1:
        freq_key = sorted(keys_in_dataset)[0]  # deterministic choice
    else:
        raise ValueError(
            "Dataset must contain at least one frequency-like dimension among: "
            "'omega', 'freq', 'wavenumber', 'wavelength', or 'period'."
        )

    freq_vals = np.asarray(dataset[freq_key].values)
    period_vals = np.asarray(dataset["period"].values)

    return freq_key, freq_vals, period_vals


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
        Must contain 'added_mass', 'radiation_damping', and either 'omega' or 'period'.
    filename: str
        Output .1 file.
    length_scale: float
        Reference length scale (L) used for normalization.

    Raises
    ------
    ValueError
        If required fields are missing or forward speed is not zero.
    """
    if "added_mass" not in dataset or "radiation_damping" not in dataset:
        raise ValueError("Missing 'added_mass' or 'radiation_damping' in dataset.")

    if not np.isclose(dataset["forward_speed"].item(), 0.0):
        raise ValueError("Forward speed must be zero for WAMIT export.")

    rho = dataset["rho"].item()
    added_mass = dataset["added_mass"]
    damping = dataset["radiation_damping"]
    omegas = dataset["omega"]
    dofs = list(added_mass.coords["influenced_dof"].values)

    # Determine main frequency coordinate
    freq_key, freq_vals, period_vals = identify_frequency_axis(dataset=added_mass)

    # Separate lines into blocks depending on period type
    period_blocks = {
        "T_zero": [],  # period == 0    → omega = inf  → PER = -1
        "T_inf": [],  # period == inf  → omega = 0    → PER = 0
        "T_regular": [],  # finite, non-zero periods
    }

    for omega, freq_val, period in zip(omegas, freq_vals, period_vals):
        for dof_i in dofs:
            for dof_j in dofs:
                j_dof, i_dof, k = get_dof_index_and_k(dof_i, dof_j)
                A = added_mass.sel(
                    {
                        freq_key: freq_val,
                        "influenced_dof": dof_i,
                        "radiating_dof": dof_j,
                    }
                ).item()
                norm = rho * (length_scale**k)
                A_norm = A / norm

                if np.isclose(period, 0.0):
                    # Case PER = -1 (omega = inf)
                    line = f"{-1.0:12.6e}\t{i_dof:5d}\t{j_dof:5d}\t{A_norm:12.6e}\n"
                    period_blocks["T_zero"].append(line)
                elif np.isinf(period):
                    # Case PER = 0 (omega = 0)
                    line = f"{0.0:12.6e}\t{i_dof:5d}\t{j_dof:5d}\t{A_norm:12.6e}\n"
                    period_blocks["T_inf"].append(line)
                else:
                    B = damping.sel(
                        {
                            freq_key: freq_val,
                            "influenced_dof": dof_i,
                            "radiating_dof": dof_j,
                        }
                    ).item()
                    B_norm = B / (omega * norm)
                    line = f"{period:12.6e}\t{i_dof:5d}\t{j_dof:5d}\t{A_norm:12.6e}\t{B_norm:12.6e}\n"
                    period_blocks["T_regular"].append((period, line))

    # Sort regular lines by increasing period
    sorted_regular = sorted(period_blocks["T_regular"], key=lambda t: t[0])
    sorted_lines = [line for _, line in sorted_regular]

    # Write to file
    with open(filename, "w") as f:
        f.writelines(period_blocks["T_zero"])
        f.writelines(period_blocks["T_inf"])
        f.writelines(sorted_lines)


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
    freq_key: str,
    freq_val: float,
    period: float,
    beta: float,
    dof: str,
    field: xarray.DataArray,
    rho: float,
    g: float,
    wave_amplitude: float,
    length_scale: float,
):
    """Write a single line for WAMIT .3 file format, using freq_key (omega or period)."""
    beta_deg = np.degrees(beta)
    i_dof = DOF_INDEX.get(dof)
    if i_dof is None:
        raise KeyError(f"DOF '{dof}' is not recognized in DOF_INDEX mapping.")

    dof_type = DOF_TYPE.get(dof, "trans")
    m = 2 if dof_type == "trans" else 3
    norm = rho * g * wave_amplitude * (length_scale**m)

    # Select value using appropriate key (omega or period)
    force = field.sel(
        {freq_key: freq_val, "wave_direction": beta, "influenced_dof": dof}
    ).item()
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
    Export excitation-like forces to a WAMIT .3-style file.

    Format:
        PER     BETA    I     Fmagnitude  Fphase  Freal  Fimaginary
    """
    forward_speed = dataset["forward_speed"].values
    if not np.isclose(forward_speed, 0.0):
        raise ValueError("Forward speed must be zero for WAMIT export.")

    if field_name not in dataset:
        raise ValueError(f"Missing field '{field_name}' in dataset.")

    field = dataset[field_name]
    rho = dataset["rho"].item()
    betas = field.coords["wave_direction"].values
    dofs = list(field.coords["influenced_dof"].values)
    g = dataset["g"].item()

    # Determine main frequency coordinate
    freq_key, freq_vals, period_vals = identify_frequency_axis(dataset=dataset)

    # Sort by increasing period
    sorted_indices = np.argsort(period_vals)
    sorted_periods = period_vals[sorted_indices]
    sorted_freqs = freq_vals[sorted_indices]

    with open(filename, "w") as f:
        for freq_val, period in zip(sorted_freqs, sorted_periods):
            # Skip WAMIT special cases
            if np.isclose(freq_val, 0.0) or np.isinf(freq_val):
                continue
            for beta in betas:
                for dof in dofs:
                    _write_wamit_excitation_line(
                        f=f,
                        freq_key=freq_key,
                        freq_val=freq_val,
                        period=period,
                        beta=beta,
                        dof=dof,
                        field=field,
                        rho=rho,
                        g=g,
                        wave_amplitude=wave_amplitude,
                        length_scale=length_scale,
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
    check_dataset_ready_for_export(dataset)

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
