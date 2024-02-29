import logging

import numpy as np

LOG = logging.getLogger(__name__)

def _get_water_depth(free_surface, water_depth, sea_bottom, default_water_depth=np.inf):
    if water_depth is None and sea_bottom is None:
        return default_water_depth
    elif water_depth is not None and sea_bottom is None:
        if water_depth <= 0.0:
            raise ValueError(f"`water_depth` should be strictly positive. Received value: {water_depth}")
        return float(water_depth)
    elif water_depth is None and sea_bottom is not None:
        LOG.warning("To uniformize notations througouth Capytaine, setting `water_depth` is preferred to `sea_bottom` since version 2.0.")
        return float(free_surface - sea_bottom)
    else:
        raise ValueError("Cannot give both a `water_depth` and a `sea_bottom`.")
