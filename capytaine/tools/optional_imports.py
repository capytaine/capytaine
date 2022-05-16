#!/usr/bin/env python
# coding: utf-8
"""Tool to import optional dependencies. Inspired by similar code in pandas."""

import importlib

def import_optional_dependency(module_name: str, package_name: str = None):
    try:
        module = importlib.import_module(module_name)
    except ImportError:
        if package_name is None:
            package_name = module_name

        message = (
            f"Missing optional dependency '{module_name}'. "
            f"Use pip or conda to install {package_name}."
        )
        raise ImportError(message) from None

    return module

def silently_import_optional_dependency(module_name: str):
    # Same as above, except it does not raise a exception when the module is not found.
    # Instead, simply returns None.
    try:
        module = importlib.import_module(module_name)
    except ImportError:
        module = None
    return module
