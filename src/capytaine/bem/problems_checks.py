import os
import logging

import numpy as np

from capytaine.tools.optional_imports import silently_import_optional_dependency

LOG = logging.getLogger(__name__)


def _check_wavelength_and_mesh_resolution(problems):
    """Display a warning if some of the problems have a mesh resolution
    that might not be sufficient for the given wavelength."""
    LOG.debug("Check wavelength with mesh resolution.")
    risky_problems = [pb for pb in problems
                      if 0.0 < pb.wavelength < pb.body.minimal_computable_wavelength]
    nb_risky_problems = len(risky_problems)
    if any(pb.body.lid_mesh is not None for pb in problems):
        mesh_str = "mesh or lid_mesh"
    else:
        mesh_str = "mesh"
    if nb_risky_problems == 1:
        pb = risky_problems[0]
        freq_type = risky_problems[0].provided_freq_type
        freq = pb.__getattribute__(freq_type)
        LOG.warning(f"Mesh resolution for {pb}:\n"
                    f"The resolution of the {mesh_str} of body {pb.body.__short_str__()} might "
                    f"be insufficient for {freq_type}={freq}.\n"
                    f"This warning appears because the largest panel of the {mesh_str} "
                    f"has radius ({pb.body.mesh_including_lid.faces_radiuses.max():.3f} m) > wavelength/8 ({pb.wavelength / 8:.3f} m)."
                    )
    elif nb_risky_problems > 1:
        freq_type = risky_problems[0].provided_freq_type
        risky_freqs = np.array([float(pb.__getattribute__(freq_type)) for pb in risky_problems])
        risky_wavelengths = np.array([pb.wavelength for pb in risky_problems])
        max_radius = max(pb.body.mesh_including_lid.faces_radiuses.max() for pb in risky_problems)
        LOG.warning(f"Mesh resolution for {nb_risky_problems} problems:\n"
                    f"The resolution of the {mesh_str} might be insufficient "
                    f"for {freq_type} ranging from {risky_freqs.min():.3f} to {risky_freqs.max():.3f}.\n"
                    f"This warning appears because the largest panel of the {mesh_str} "
                    f"has radius ({max_radius:.3f} m) > wavelength/8 ({risky_wavelengths.min() / 8:.3f} to {risky_wavelengths.max() / 8:.3f} m)."
                    )


def _check_wavelength_and_water_depth(problems):
    """Check for problems with finite depth that could as well be computed in infinite depth."""
    LOG.debug("Check wavelength and water depth.")
    def check(pb):
        return ((pb.water_depth < np.inf)
                and (pb.water_depth > 5*pb.wavelength)
                and (pb.water_depth > -2*pb.body.mesh.vertices[:, 2].min()))
    filtered_pbs = [pb for pb in problems if check(pb)]
    nb_filtered_pbs = len(filtered_pbs)
    if nb_filtered_pbs > 1:
        wd = list(np.unique([pb.water_depth for pb in filtered_pbs]))
        if len(wd) == 1:
            wd = wd[0]

        LOG.warning(f"Water depth for {nb_filtered_pbs} problems:\n"
                    f"Very deep finite water depth have been set (`water_depth={wd}`).\n"
                    "Computation might be faster by using the infinite water depth instead: `water_depth=np.inf`."
                    )


def _check_wavelength_and_irregular_frequencies(problems):
    """Display a warning if some of the problems might encounter irregular frequencies."""
    LOG.debug("Check wavelength with estimated irregular frequency.")
    risky_problems = [pb for pb in problems
                      if pb.free_surface != np.inf and
                      pb.body.first_irregular_frequency_estimate(g=pb.g) < pb.omega < np.inf]
    nb_risky_problems = len(risky_problems)
    if nb_risky_problems >= 1:
        if any(pb.body.lid_mesh is None for pb in problems):
            recommendation = "Setting a lid for the floating body is recommended."
        else:
            recommendation = "The lid might need to be closer to the free surface."
        if nb_risky_problems == 1:
            pb = risky_problems[0]
            freq_type = risky_problems[0].provided_freq_type
            freq = pb.__getattribute__(freq_type)
            LOG.warning(f"Irregular frequencies for {pb}:\n"
                        f"The body {pb.body.__short_str__()} might display irregular frequencies "
                        f"for {freq_type}={freq}.\n"
                        + recommendation
                        )
        elif nb_risky_problems > 1:
            freq_type = risky_problems[0].provided_freq_type
            freqs = np.array([float(pb.__getattribute__(freq_type)) for pb in risky_problems])
            LOG.warning(f"Irregular frequencies for {nb_risky_problems} problems:\n"
                        "Irregular frequencies might be encountered "
                        f"for {freq_type} ranging from {freqs.min():.3f} to {freqs.max():.3f}.\n"
                        + recommendation
                        )


def _check_ram(problems, engine, n_jobs=1):
    """Display a warning if the RAM estimation is larger than a certain limit."""
    LOG.debug("Check RAM estimation.")
    psutil = silently_import_optional_dependency("psutil")
    if psutil is None:
        ram_limit = 8
    else :
        ram_limit = psutil.virtual_memory().total / (1024**3) * 0.3

    if n_jobs == - 1:
        n_jobs = os.cpu_count()

    estimated_peak_memory = n_jobs*max(engine.compute_ram_estimation(pb) for pb in problems)

    if estimated_peak_memory < 0.5:
        LOG.info("Estimated peak RAM usage: <1 GB.")

    elif estimated_peak_memory < ram_limit:
        LOG.info(f"Estimated peak RAM usage: {int(np.ceil(estimated_peak_memory))} GB.")

    else:
        LOG.warning(f"Estimated peak RAM usage: {int(np.ceil(estimated_peak_memory))} GB.")
