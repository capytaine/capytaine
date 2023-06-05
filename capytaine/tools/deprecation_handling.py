import logging

import numpy as np

LOG = logging.getLogger(__name__)

def _get_water_depth(free_surface, water_depth, sea_bottom, default_water_depth=np.infty):
    if water_depth is None and sea_bottom is None:
        return default_water_depth
    elif water_depth is not None and sea_bottom is None:
        return float(water_depth)
    elif water_depth is None and sea_bottom is not None:
        LOG.warning("Deprecation warning: please prefer giving a `water_depth` rather than a `sea_bottom`.")
        return float(free_surface - sea_bottom)
    else:
        raise ValueError("Cannot give both a `water_depth` and a `sea_bottom`.")
