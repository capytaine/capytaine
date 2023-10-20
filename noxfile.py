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

def run_tests(session):
    with session.chdir(session.create_tmp()):
        session.run("python", "-m", "pytest", os.path.join(NOXFILE_DIR, "pytest"))
        session.run('python', '-c', '"import capytaine; print(capytaine.__version__)"')
        session.run('capytaine', '--help')
        for example_file in ["custom_dofs.py", "haskind.py", "multibody.py"]:
            session.run('python', os.path.join(NOXFILE_DIR, "examples", example_file))


@nox.session
def build_and_test_on_latest_env(session):
    session.install(".[test]")
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
    run_tests(session)
