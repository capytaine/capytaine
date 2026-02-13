set shell := ["bash", "-c"]
set windows-shell := ["powershell.exe", "-c"]

default:
    just --list

# Run the recipe below in a clean virtual environment to install Capytaine in
# editable mode.
editable_install:
    uv pip install -r pyproject.toml \
        --extra optional \
        --group editable_install \
        --group dev
    uv pip install --no-build-isolation --editable .

# Define the temporary directory differently based on OS
TEMP_DIR := if os_family() == 'windows' {
    `$temp = Join-Path $env:TEMP ('capytainedev-' + [System.IO.Path]::GetRandomFileName()); New-Item -ItemType Directory -Path $temp -Force | Out-Null; Write-Output $temp`
} else {
    `mktemp --directory --tmpdir capytainedev-XXXXX`
}

ENV := "MPLBACKEND=pdf CAPYTAINE_PROGRESS_BAR=False CAPYTAINE_CACHE_DIR=" + TEMP_DIR + "/cache/"

TEST_DIR := justfile_directory() / 'pytest'
NEMOH_CASES := TEST_DIR / 'Nemoh_verification_cases' / 'Cylinder'

EXAMPLES_DIR := justfile_directory() / 'docs' / 'examples' / 'src'
EXAMPLES_FILES := ' \
    A1_single_body_hydrodynamics.py \
    A2_multibody.py \
    A3_finite_depth_flap.py \
    A4_forward_speed_on_vertical_cylinder.py \
    A7_custom_dofs.py \
    A8_convergence_study.py \
    A10_elasticity_of_beam.py \
    A11_parametric_study_depth.py \
    B1_pressure_on_hull.py \
    B2_haskind.py \
    B3_free_surface_elevation.py \
    B4_kochin.py \
    B5_plot_velocity_in_domain.py \
    B8_pressure_infinite_frequency.py \
    C5_plot_influence_matrix.py \
    C8_compare_Green_functions.py \
    C9_custom_Green_function.py \
'

## Reason for skipping some example files:
# A5_benchmark_plane_symmetries.py \
# A6_benchmark_axisymmetric_mesh.py \
# A9_test_irregular_frequency_removal.py \ # Slow
# B6_animate_free_surface.py  \           # Requires VTK
# B7_boat_animation.py  \                 # Requires VTK
# C6_axisymmetric_buoy.py  \              # Requires VTK
# C10_custom_linear_solver_on_gpu.py \    # Requires torch

# Run the test suite and the example files assuming a virtual environment with
# Capytaine installed has been activated
[unix]
_test:
    #!/usr/bin/env bash
    set -euxo pipefail
    uv pip list | grep -Ei "numpy|xarray|capytaine" || true
    cd {{TEMP_DIR}}
    export {{ENV}}
    python -c "import capytaine; print(capytaine.__version__)"
    python -m pytest {{TEST_DIR}}
    capytaine --help
    capytaine {{NEMOH_CASES}}/Nemoh.cal
    capytaine {{NEMOH_CASES}}/Nemoh_v3.cal
    set +x  # Stop displaying commands
    for f in {{EXAMPLES_FILES}}; do
        echo -e "\n---> Running $f";
        python {{EXAMPLES_DIR}}/$f;
    done

# Simplified version of the above, working on Windows
[windows]
_test:
    #! powershell
    cd {{TEMP_DIR}}
    python -c "import capytaine; print(capytaine.__version__)"
    python -m pytest {{TEST_DIR}}
    capytaine --help
    capytaine {{NEMOH_CASES}}/Nemoh.cal
    capytaine {{NEMOH_CASES}}/Nemoh_v3.cal

_install_and_test:
    uv pip install --no-deps --no-build-isolation .
    just _test


# In the recipes below, we have
#
#    --no-editable --with "capytaine @ ." -- just _test
#
# or
#
#    --no-editable -- just _install_and_test
#
# because the default behavior of uv is to install the local project in
# editable mode, but we don't want that because of an incompatibility with
# meson-python (https://github.com/astral-sh/uv/issues/10214)
#
# besides, we have --no-default-groups because uv loads by default the `dev`
# dependency-group from pyproject.toml, but we don't need it here and it can
# actually cause issue in CI.

test_in_latest_env:
    uv run \
        --isolated --no-default-groups \
        --no-editable --with "capytaine[optional] @ ." \
        --group test \
        -- \
        just _test

# In the recipe below,
# "--index-strategy unsafe-best-match" means uv should not ignore wheels
# from PyPI during universal resolution

test_in_nightly_env:
    uv run \
        --isolated --no-default-groups \
        --python 3.14 \
        --pre --extra-index-url https://pypi.anaconda.org/scientific-python-nightly-wheels/simple \
        --index-strategy unsafe-best-match \
        --no-editable --with "capytaine @ ." \
        --group test \
        -- \
        just _test


test_in_py38_reference_env:
    uv run \
        --isolated --no-default-groups \
        --python 3.8 \
        --no-editable \
        --with-requirements {{TEST_DIR}}/envs/2023-08-01-py3.8.txt \
        -- \
        just _install_and_test

test_in_py313_reference_env:
    uv run \
        --isolated --no-default-groups \
        --python 3.13 \
        --no-editable \
        --with-requirements {{TEST_DIR}}/envs/2025-11-25-py3.13.txt \
        -- \
        just _install_and_test

# How the requirements files from the above recipes where generated.
create_test_env_file python="3.8" date="2023-08-01":
    uv pip compile \
        pyproject.toml \
        --python-version {{python}} \
        --group test \
        --group editable_install \
        --extra optional \
        --exclude-newer {{date}} \
        -o {{TEST_DIR}}/envs/{{date}}-py{{python}}.txt


# Compile the Fortran code without parallelism for easier reading of the errors.
test_fortran_compilation:
    # It is assumed that meson and ninja are already installed (e.g. with editable_install).
    meson setup --wipe {{TEMP_DIR}}/build && meson compile -C {{TEMP_DIR}}/build -j 1


build_docs:
    uv run \
        --isolated \
        --no-editable --with "capytaine[optional] @ ." \
        --group docs \
        -- \
        make --directory="./docs/"


clean:
    rm -f src/capytaine/green_functions/libs/*.so
    rm -rf build/
    rm -rf dist/
    rm -rf src/capytaine.egg-info/
    rm -rf docs/_build
    rm -rf .pytest_cache/
    rm -rf .venv/
    rm -rf __pycache__ */__pycache__ */*/__pycache__ */*/*/__pycache__
    rm -rf $HOME/.cache/capytaine/*
    rm -rf /tmp/capytainedev-*/
