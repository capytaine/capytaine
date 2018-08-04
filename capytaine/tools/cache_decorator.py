# coding: utf-8
# This file is part of "Capytaine" (https://github.com/mancellin/capytaine).
# It has been written by Matthieu Ancellin and is released under the terms of the GPLv3 license.

import logging

import numpy as np

from capytaine.tools.max_length_dict import MaxLengthDict

LOG = logging.getLogger(__name__)


def keep_in_cache(cache_name="unnamed_cache"):
    def _decorator(matrix_building_function):
        def matrix_build_function_with_cache(self, mesh1, *args, _rec_depth=(1,)):
            if self.use_cache and mesh1 not in self.__cache__[cache_name]:
                self.__cache__[cache_name][mesh1] = MaxLengthDict({}, max_length=int(np.product(_rec_depth)))
                LOG.debug("\t" * len(_rec_depth) +
                          f"\tCreate {cache_name} cache (max_length={int(np.product(_rec_depth))}) for {mesh1.name}")

            mesh2 = args[0]
            if not self.use_cache or args not in self.__cache__[cache_name][mesh1]:
                LOG.debug("\t" * len(_rec_depth) +
                          f"\tComputing {cache_name} of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name}")
                matrices = matrix_building_function(self, mesh1, *args, _rec_depth)
                if self.use_cache:
                    self.__cache__[cache_name][mesh1][args] = matrices
                return matrices

            else:
                LOG.debug("\t" * len(_rec_depth) +
                          f"\tRetrieving stored matrix {cache_name} of {mesh1.name} on {'itself' if mesh2 is mesh1 else mesh2.name}")
                return self.__cache__[cache_name][mesh1][args]

        return matrix_build_function_with_cache
    return _decorator


