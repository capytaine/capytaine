project('libDelhommeau bindings for Capytaine',
  ['c', 'fortran'],
  version : '2.0.dev0',
  default_options : ['warning_level=2'])


py_mod = import('python')
py3 = py_mod.find_installation('python3')
py3_dep = py3.dependency()
message(py3.path())
message(py3.get_install_dir())

incdir_numpy = run_command(py3,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

libDelhommeau_src = [
      'capytaine/green_functions/libDelhommeau/src/float64.f90',
      'capytaine/green_functions/libDelhommeau/src/constants.f90',
      'capytaine/green_functions/libDelhommeau/src/Delhommeau_integrals.f90',
      'capytaine/green_functions/libDelhommeau/src/old_Prony_decomposition.f90',
      'capytaine/green_functions/libDelhommeau/src/Green_Rankine.f90',
      'capytaine/green_functions/libDelhommeau/src/Green_wave.f90',
      'capytaine/green_functions/libDelhommeau/src/matrices.f90'
]

libDelhommeau_f2py_wrapper = custom_target('Delhommeau_float64module.c',
                                           input : libDelhommeau_src,
                                           output : ['Delhommeau_float64module.c', 'Delhommeau_float64-f2pywrappers2.f90'],
                                           command : [ py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'Delhommeau_float64', '--lower'])

install_subdir('capytaine', install_dir: py3.get_install_dir(pure: false))
# Alternatively: would need to list all the .py files and use py3.install_source()

py3.extension_module('Delhommeau_float64',
           libDelhommeau_src,
           libDelhommeau_f2py_wrapper,
           incdir_f2py+'/fortranobject.c',
           include_directories: inc_np,
           dependencies : py3_dep,
           fortran_args : ['-cpp', '-fopenmp'],
           link_args : ['-fopenmp'],
           install : true,
           subdir : 'capytaine/green_functions/libs/')

