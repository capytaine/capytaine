project('capytaine',
  ['c', 'fortran'],
  version : run_command('capytaine/__about__.py', check: true).stdout().strip(),
  default_options : ['warning_level=2'])

py_mod = import('python')
py3 = py_mod.find_installation('python3')
py3_dep = py3.dependency()
message('Python executable: @0@'.format(py3.path()))
numpy_version = run_command(py3, ['-c', 'import numpy; print(numpy.__version__)']).stdout().strip()
message('Numpy version: @0@'.format(numpy_version))

incdir_numpy = run_command(py3,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include() if hasattr(numpy.f2py, "get_include") else os.path.join(os.path.dirname(numpy.f2py.__file__), "src"))'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

libDelhommeau_src = [
      'capytaine/green_functions/libDelhommeau/src/constants.f90',
      'capytaine/green_functions/libDelhommeau/src/Delhommeau_integrals.f90',
      'capytaine/green_functions/libDelhommeau/src/old_Prony_decomposition.f90',
      'capytaine/green_functions/libDelhommeau/src/Green_Rankine.f90',
      'capytaine/green_functions/libDelhommeau/src/Green_wave.f90',
      'capytaine/green_functions/libDelhommeau/src/matrices.f90',
]

float32 = ['capytaine/green_functions/libDelhommeau/src/float32.f90',]
float64 = ['capytaine/green_functions/libDelhommeau/src/float64.f90',]

Delhommeau_float32_f2py_wrapper = custom_target('Delhommeau_float32module.c',
                                           input : float32 + libDelhommeau_src,
                                           output : ['Delhommeau_float32module.c', 'Delhommeau_float32-f2pywrappers2.f90'],
                                           command : [ py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'Delhommeau_float32', '--lower'])

Delhommeau_float64_f2py_wrapper = custom_target('Delhommeau_float64module.c',
                                           input : float64 + libDelhommeau_src,
                                           output : ['Delhommeau_float64module.c', 'Delhommeau_float64-f2pywrappers2.f90'],
                                           command : [ py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'Delhommeau_float64', '--lower'])

install_subdir(
  'capytaine',
  exclude_directories: ['green_functions/libDelhommeau/'],  # The Fortran source files do not need to be distributed in the wheel
  install_dir: py3.get_install_dir(pure: false),
  )
# Alternatively: would need to list all the .py files and use py3.install_source()

if build_machine.system() == 'darwin'
  omp = []  # Waiting for https://github.com/mesonbuild/meson/issues/7435
else
  omp = dependency('openmp')
endif

args = ['-cpp']
link_args = ['-lquadmath']  # Required for the Windows build with cibuildwheel
deps = [py3_dep, omp]

py3.extension_module('Delhommeau_float32',
           float32, libDelhommeau_src, Delhommeau_float32_f2py_wrapper, incdir_f2py+'/fortranobject.c',
           include_directories: inc_np,
           dependencies : deps,
           fortran_args : args,
           link_args : link_args,
           install : true,
           subdir : 'capytaine/green_functions/libs/')

py3.extension_module('Delhommeau_float64',
           float64, libDelhommeau_src, Delhommeau_float64_f2py_wrapper, incdir_f2py+'/fortranobject.c',
           include_directories: inc_np,
           dependencies : deps,
           fortran_args : args,
           link_args : link_args,
           install : true,
           subdir : 'capytaine/green_functions/libs/')
