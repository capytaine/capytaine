import os
import nox

os.environ.update({"PDM_IGNORE_SAVED_PYTHON": "1"})

@nox.session
def build_and_test_on_latest_env(session):
    session.install(".[ci]")
    noxfile_dir = os.path.dirname(__file__)
    with session.chdir(session.create_tmp()):
        session.run("pytest", os.path.join(noxfile_dir, "pytest"))
        session.run('python', '-c', '"import capytaine; print(capytaine.__version__)"')
        session.run('capytaine', '--help')
        for example_file in ["custom_dofs.py", "haskind.py", "multibody.py"]:
            session.run('python', os.path.join(noxfile_dir, "examples", example_file))

@nox.session
@nox.parametrize("env_file", ["2023-10-02"])
def build_and_test_on_locked_env(session, env_file):
    session.run_always('pdm', 'sync', '--no-self', '-L', f"pytest/envs/{env_file}.lock", external=True)
    session.install("--no-build-isolation", "--no-deps", ".")
    noxfile_dir = os.path.dirname(__file__)
    with session.chdir(session.create_tmp()):
        session.run("pytest", os.path.join(noxfile_dir, "pytest"))
        session.run('python', '-c', '"import capytaine; print(capytaine.__version__)"')
        session.run('capytaine', '--help')
        for example_file in ["custom_dofs.py", "haskind.py", "multibody.py"]:
            session.run('python', os.path.join(noxfile_dir, "examples", example_file))

