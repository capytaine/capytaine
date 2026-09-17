# Copyright 2026 Capytaine developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Adapted from https://github.com/platformdirs/platformdirs (MIT Licensed)
"""

import os
import sys
from pathlib import Path

from capytaine import __version__


def cache_directory():
    if "CAPYTAINE_CACHE_DIR" in os.environ:
        path = os.path.join(os.environ["CAPYTAINE_CACHE_DIR"], __version__)
    elif sys.platform == "win32":  # Windows
        path = os.path.normpath(os.environ.get("LOCALAPPDATA"))
        path = os.path.join(path, "capytaine", "Cache", __version__)
    elif sys.platform == "darwin":  # MacOS
        path = os.path.expanduser("~/Library/Caches")
        path = os.path.join(path, "capytaine", __version__)
    else:
        path = os.environ.get("XDG_CACHE_HOME", "")
        if path.strip() == "":
            path = os.path.expanduser("~/.cache")
        path = os.path.join(path, "capytaine", __version__)
    Path(path).mkdir(parents=True, exist_ok=True)
    return path
