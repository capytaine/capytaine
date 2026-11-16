"""Check that the compiled Fortran core of Capytaine has been built with OpenMP.

Prints the shared libraries the extension module is linked to and compares the
time to evaluate the Green function with one and several OpenMP threads.
"""

import os
import platform
import subprocess
import sys
from pathlib import Path

import capytaine

so_files = list((Path(capytaine.__file__).parent / "green_functions").glob("Delhommeau_float64*"))
print("Extension module:", so_files)
if platform.system() == "Darwin":
    linked = subprocess.run(["otool", "-L", str(so_files[0])], capture_output=True, text=True).stdout
elif platform.system() == "Linux":
    linked = subprocess.run(["ldd", str(so_files[0])], capture_output=True, text=True).stdout
else:
    linked = ""
print(linked)

TIMING = """
import time
import numpy as np
import capytaine as cpt
mesh = cpt.mesh_sphere(radius=1.0, center=(0, 0, -2), resolution=(40, 40))
gf = cpt.Delhommeau()
gf.evaluate(mesh, mesh, wavenumber=1.0)  # Warm-up
t0 = time.perf_counter()
for _ in range(3):
    gf.evaluate(mesh, mesh, wavenumber=1.0)
print(time.perf_counter() - t0)
"""

def timing(nb_threads):
    env = {**os.environ, "OMP_NUM_THREADS": str(nb_threads)}
    out = subprocess.run([sys.executable, "-c", TIMING], env=env, capture_output=True, text=True, check=True)
    return float(out.stdout.strip().splitlines()[-1])

nb_cpus = os.cpu_count()
t1 = timing(1)
tn = timing(nb_cpus)
print(f"cpu_count={nb_cpus}  time with 1 thread: {t1:.3f}s  time with {nb_cpus} threads: {tn:.3f}s  speedup: {t1/tn:.2f}")

if platform.system() != "Windows":
    assert "gomp" in linked or "omp" in linked, "Extension module not linked to an OpenMP runtime"
assert t1 / tn > 1.3, "No speedup with OpenMP threads"
print("OpenMP OK")
