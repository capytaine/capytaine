import logging
from collections import defaultdict


LOG = logging.getLogger(__name__)


def display_grouped_errors(results):
    """Displays errors that occur during the solver execution and groups them according
    to the problem type and exception type for easier reading."""
    failed_results = defaultdict(list)
    for res in results:
        if hasattr(res, "exception") and hasattr(res, "problem"):
            key = (type(res.exception), str(res.exception), res.problem.omega, res.problem.water_depth, res.problem.forward_speed)
            failed_results[key].append(res.problem)

    for (exc_type, exc_msg, omega, water_depth, forward_speed), problems in failed_results.items():
        nb = len(problems)
        if nb > 1:
            if forward_speed != 0.0:
                LOG.warning("Skipped %d problems for body=%s, omega=%s, water_depth=%s, forward_speed=%s\nbecause of %s(%r)",
                            nb, problems[0].body.__short_str__(), omega, water_depth, forward_speed, exc_type.__name__, exc_msg)
            else:
                LOG.warning("Skipped %d problems for body=%s, omega=%s, water_depth=%s\nbecause of %s(%r)",
                            nb, problems[0].body.__short_str__(), omega, water_depth, exc_type.__name__, exc_msg)
        else:
            LOG.warning("Skipped %s\nbecause of %s(%r)", problems[0], exc_type.__name__, exc_msg)
