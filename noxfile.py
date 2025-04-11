import os
import sys
import tempfile
import nox

# Create virtual environments in a temporary directory somewhere else because
# meson-python does not like local virtualenvironments
# (https://github.com/capytaine/capytaine/issues/396)
tempdir = tempfile.mkdtemp(prefix="nox-capytaine-")
nox.options.envdir = os.path.join(tempdir, "venvs")

ENV = {
        'MPLBACKEND': 'pdf',
        'CAPYTAINE_CACHE_DIR': os.path.join(tempdir, "capytaine_cache"),
        'CAPYTAINE_PROGRESS_BAR': 'False'
        }

EXAMPLE_FILES = [
    # "quickstart.py",
    "A1_radiation_cylinder.py",
    "A2_multibody.py",
    "A3_finite_depth_cylinder.py",
    "A4_custom_dofs.py",
    "A5_convergence_study.py",
    # "A6_irregular_frequency_removal.py",      # Slow
    "A7_elasticity_of_beam.py",
    "B1_pressure_on_hull.py",
    "B2_haskind.py",
    "B3_free_surface_elevation.py",
    "B4_kochin.py",
    "B5_plot_velocity_in_domain.py",
    # "B6_animate_free_surface.py",             # Requires VTK
    # "B7_boat_animation.py",                   # Requires VTK
    "C5_plot_influence_matrix.py",
    # "C6_axisymmetric_buoy.py",                # Requires VTK
    # "C7_h_matrices_with_preconditionner.py",  # Slow
    "C8_compare_Green_functions.py",
    "C9_custom_Green_function.py",
    # "C10_custom_linear_solver_on_gpu.py",     # Requires torch
        ]

NOXFILE_DIR = os.path.dirname(__file__)

NEMOH_CASES = os.path.join(NOXFILE_DIR, "pytest", "Nemoh_verification_cases", "Cylinder")

def run_tests(session):
    with session.chdir(session.create_tmp()):
        session.run("python", "-m", "pytest", os.path.join(NOXFILE_DIR, "pytest"), env=ENV)
        session.run('python', '-c', '"import capytaine; print(capytaine.__version__)"')
        session.run('capytaine', '--help')
        session.run('capytaine', os.path.join(NEMOH_CASES, "Nemoh.cal"), env=ENV)
        session.run('capytaine', os.path.join(NEMOH_CASES, "Nemoh_v3.cal"), env=ENV)
        for example_file in EXAMPLE_FILES:
            session.run('python', os.path.join(NOXFILE_DIR, "docs", "examples", "src", example_file), env=ENV)


@nox.session
def build_and_test_on_latest_env(session):
    # By default, pip will install the latest dependencies compatible with the
    # constraints in pyproject.toml.
    session.install(".[test,optional]")
    run_tests(session)


@nox.session
def editable_build_and_test_on_latest_env(session):
    session.install("-r", "editable_install_requirements.txt")
    session.install("--no-build-isolation", "--editable", ".[test,optional]")
    run_tests(session)


@nox.session(python=['3.8', '3.12'])
def build_and_test_on_locked_env(session):
    if session.python == '3.8':
        env_file = "2023-08-01-py3.8.txt"
        # Lock file was created with the following command
        # PY=3.8 DATE=2023-08-01 uv pip compile \
        # pyproject.toml editable_install_requirements.txt \
        # --python-version $PY --exclude-newer $DATE \
        # --extra optional --extra test \
        # -o pytest/envs/$DATE-py$PY.txt
        # Older date where not possible to reach because of the joblib>=1.3 requirement.
    elif session.python == '3.12':
        env_file = "2024-10-22-py3.12.txt"
        # PY=3.12 DATE=2024-10-22 uv pip compile \
        # pyproject.toml editable_install_requirements.txt \
        # --python-version $PY --exclude-newer $DATE \
        # --extra optional --extra test \
        # -o pytest/envs/$DATE-py$PY.txt
    else:
        # On CI, this session is only run on Python 3.8 and 3.12
        # (see .github/workflows/test_new_commits.yaml)
        # This fallback might be useful for local tests:
        env_file = "2024-04-08-py3.12.txt"

    session.install('-r', f"pytest/envs/{env_file}")

    # We install without build isolation in order to control also the build environment.
    # There might be other ways to do that.
    session.install("--no-deps", "--no-build-isolation", ".")

    run_tests(session)


@nox.session
def build_and_test_on_nightly_builds(session):
    # https://scientific-python.org/specs/spec-0004/
    session.install("--pre", "--upgrade",
                    "--extra-index-url", "https://pypi.anaconda.org/scientific-python-nightly-wheels/simple",
                    "numpy", "scipy", "pandas", "xarray", "rich")
    session.install("meson-python", "ninja", "charset-normalizer")
    session.install("--no-deps", "--no-build-isolation", ".")
    session.install("pytest", "matplotlib")
    run_tests(session)
