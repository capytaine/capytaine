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
