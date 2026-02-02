project('capytaine',
  ['c', 'fortran'],
  meson_version : '>= 1.0.0',
  version : '3.0.0a1',
  default_options : ['warning_level=2'])

py_mod = import('python')
py3 = py_mod.find_installation('python3')
py3_dep = py3.dependency()
message('Python executable: @0@'.format(py3.full_path()))
numpy_version = run_command(py3, ['-c', 'import numpy; print(numpy.__version__)'], check: true).stdout().strip()
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

float32 = ['libs/libDelhommeau/src/float32.f90',]
float64 = ['libs/libDelhommeau/src/float64.f90',]

src_to_interface_with_python = [
      'libs/libDelhommeau/src/constants.f90',
      'libs/libDelhommeau/src/Delhommeau_integrals.f90',
      'libs/libDelhommeau/src/old_Prony_decomposition.f90',
      'libs/LiangWuNoblesse/LiangWuNoblesseWaveTerm.f90',
      'libs/libDelhommeau/src/Green_Rankine.f90',
      'libs/libDelhommeau/src/Green_wave.f90',
      'libs/libDelhommeau/src/matrices.f90',
      'libs/libDelhommeau/src/interface.f90',
]

other_src = [
      'libs/FinGreen3D/FinGreen3D.f90',
  ]

Delhommeau_float32_f2py_wrapper = custom_target('Delhommeau_float32module.c',
                                           input : float32 + src_to_interface_with_python,
                                           output : ['Delhommeau_float32module.c', 'Delhommeau_float32-f2pywrappers2.f90'],
                                           command : [ py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'Delhommeau_float32', '--lower'])

Delhommeau_float64_f2py_wrapper = custom_target('Delhommeau_float64module.c',
                                           input : float64 + src_to_interface_with_python,
                                           output : ['Delhommeau_float64module.c', 'Delhommeau_float64-f2pywrappers2.f90'],
                                           command : [ py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'Delhommeau_float64', '--lower'])

install_subdir(
  'src/capytaine',
  exclude_files: ['__about__.py.in'],
  install_dir: py3.get_install_dir(pure: false),
  )
# Alternatively: would need to list all the .py files and use py3.install_source()

# Write version and build information to be accessible at runtime
fc = meson.get_compiler('fortran')
configure_file(
  input: 'src/capytaine/__about__.py.in',
  output: '__about__.py',
  configuration: {
    'version': meson.project_version(),
    'build_time': run_command('date', '+%Y-%m-%dT%H:%M:%SZ').stdout().strip(),
    'compiler': fc.get_id(),
    'compiler_version': fc.version(),
  },
  install: true,
  install_dir: py3.get_install_dir(pure: false) / 'capytaine',
)

if build_machine.system() == 'darwin'
  omp = []  # Waiting for https://github.com/mesonbuild/meson/issues/7435
else
  omp = dependency('openmp')
endif

args = ['-cpp']
optional_dependencies_flags = ['-DLIANGWUNOBLESSE_OPTIONAL_DEPENDENCY', '-DFINGREEN3D_OPTIONAL_DEPENDENCY']
link_args = ['-lquadmath']  # Required for the Windows build with cibuildwheel
deps = [py3_dep, omp]

py3.extension_module('Delhommeau_float32',
           float32, src_to_interface_with_python, Delhommeau_float32_f2py_wrapper, incdir_f2py+'/fortranobject.c',
           include_directories: inc_np,
           dependencies : deps,
           fortran_args : args,
           link_args : link_args,
           install : true,
           subdir : 'capytaine/green_functions/')

py3.extension_module('Delhommeau_float64',
           float64, other_src, src_to_interface_with_python, Delhommeau_float64_f2py_wrapper, incdir_f2py+'/fortranobject.c',
           include_directories: inc_np,
           dependencies : deps,
           fortran_args : args + optional_dependencies_flags,
           link_args : link_args,
           install : true,
           subdir : 'capytaine/green_functions/')
