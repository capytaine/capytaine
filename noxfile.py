import os
import tempfile
import nox

# Create virtual environments in a temporary directory somewhere else because
# meson-python does not like local virtualenvironments
# (https://github.com/capytaine/capytaine/issues/396)
nox.options.envdir = os.path.join(tempfile.gettempdir(), "nox-capytaine")

os.environ.update({"PDM_IGNORE_SAVED_PYTHON": "1"})
os.environ.update({"PDM_USE_VENV": "1"})

NOXFILE_DIR = os.path.dirname(__file__)

EXAMPLE_FILES = ["compare_Green_functions.py", "convergence_study.py", "custom_dofs.py", "custom_Green_function.py", "finite_depth_cylinder.py", "free_surface_elevation.py", "haskind.py", "kochin.py", "multibody.py", "plot_influence_matrix.py", "plot_velocity_in_domain.py", "radiation_cylinder.py"]

def run_tests(session):
    with session.chdir(session.create_tmp()):
        session.run("python", "-m", "pytest", os.path.join(NOXFILE_DIR, "pytest"))
        session.run('python', '-c', '"import capytaine; print(capytaine.__version__)"')
        session.run('capytaine', '--help')
        for example_file in EXAMPLE_FILES:
            session.run('python', os.path.join(NOXFILE_DIR, "examples", example_file), env={'MPLBACKEND': 'pdf'})


@nox.session
def build_and_test_on_latest_env(session):
    session.install(".[test,optional]")
    run_tests(session)


@nox.session
@nox.parametrize("env_file", ["2023-10-02"])
def build_and_test_on_locked_env(session, env_file):
    session.install("pdm")
    session.run_always('pdm', 'sync', '--no-self', '-L', f"pytest/envs/{env_file}.lock")
    # Lock file was created with
    # pdm lock -d -G build -G test -L pytest/envs/2023-10-02.lock
    session.install("--no-build-isolation", "--no-deps", "-e", ".")
    # Editable install using pip and not pdm because pdm's editable install
    # does not seems to work with meson-python.
    session.install("matplotlib")  # Need to be installed in the lock file in the future.
    run_tests(session)
