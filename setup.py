#!/usr/bin/env python
# coding: utf-8

from numpy.distutils.core import Extension, setup

def libDelhommeau_src(precision):
    return [
        "capytaine/green_functions/libDelhommeau/src/{}.f90".format(precision),
        "capytaine/green_functions/libDelhommeau/src/constants.f90",
        "capytaine/green_functions/libDelhommeau/src/Delhommeau_integrals.f90",
        "capytaine/green_functions/libDelhommeau/src/old_Prony_decomposition.f90",
        "capytaine/green_functions/libDelhommeau/src/Green_Rankine.f90",
        "capytaine/green_functions/libDelhommeau/src/Green_wave.f90",
        "capytaine/green_functions/libDelhommeau/src/matrices.f90",
    ]

extensions_names_and_extra_arguments = [
        ("capytaine.green_functions.libs.Delhommeau", []),
        ("capytaine.green_functions.libs.XieDelhommeau", ["-DXIE_CORRECTION"]),
        ]

extensions_modules = [
        Extension(
            name=name + "_" + precision,
            sources=libDelhommeau_src(precision),
            extra_f90_compile_args=['-fopenmp', '-cpp'] + extra_args,
            extra_link_args=['-fopenmp'],
            # # Uncomment the following lines to get more verbose output from f2py.
            # define_macros=[
            #     ('F2PY_REPORT_ATEXIT', 1),
            #     ('F2PY_REPORT_ON_ARRAY_COPY', 1),
            # ],
        )
        for (name, extra_args) in extensions_names_and_extra_arguments
        for precision in ["float32", "float64"]
        ]


if __name__ == "__main__":
    setup(ext_modules=extensions_modules)
