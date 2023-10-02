import os
import nox

os.environ.update({"PDM_IGNORE_SAVED_PYTHON": "1"})

@nox.session
def build_and_test_on_latest_env(session):
    session.install(".[ci]")
    session.run('python', '-c', '"import capytaine; print(capytaine.__version__)"')
    session.run("pytest")

@nox.session
@nox.parametrize("env_file", ["2023-10-02"])
def build_and_test_on_locked_env(session, env_file):
    session.run_always('pdm', 'sync', '--no-self', '-L', f"pytest/envs/{env_file}.lock", external=True)
    session.install(".")
    session.run('python', '-c', '"import capytaine; print(capytaine.__version__)"')
    session.run('pytest')

