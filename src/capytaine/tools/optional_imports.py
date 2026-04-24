"""Tool to import optional dependencies. Inspired by similar code in pandas."""

import importlib

def import_optional_dependency(
    module_name: str,
    package_name: str = None,
    error_message: str = None,
):
    """Return an imported module or raises an error.
    The error message can be customized either by passing the name of the package to install
    (if it's different from the module name),
    or by passing a full error message."""
    try:
        module = importlib.import_module(module_name)
    except ImportError:
        if error_message is None:
            if package_name is None:
                package_name = module_name
            error_message = (
                f"Missing optional dependency '{module_name}'. "
                f"Use pip or conda to install {package_name}."
            )
        raise ImportError(error_message) from None

    return module

def silently_import_optional_dependency(module_name: str):
    # Same as above, except it does not raise a exception when the module is not found.
    # Instead, simply returns None.
    try:
        module = importlib.import_module(module_name)
    except ImportError:
        module = None
    return module
