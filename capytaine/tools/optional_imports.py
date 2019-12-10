#!/usr/bin/env python
# coding: utf-8
"""Tool to import optional dependencies. Inspired by similar code in pandas."""

import importlib

def import_optional_dependency(name: str):
    try:
        module = importlib.import_module(name)
    except ImportError:
        message = (
            "Missing optional dependency '{name}'."
            "Use pip or conda to install {name}."
        )
        raise ImportError(message.format(name=name)) from None

    return module
