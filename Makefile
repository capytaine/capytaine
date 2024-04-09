install:
	pip install .

develop:
	pip install -r editable_install_requirements.txt
	pip install --no-build-isolation --editable .

TEMP_DIR := $(shell mktemp -d)
test_fortran_compilation:
	# Compile the Fortran code without parallelism for easier reading of the errors.
	# It is assumed that meson and ninja are already installed.
	meson setup --wipe $(TEMP_DIR) && meson compile -C $(TEMP_DIR) -j 1

test:
	# Build and test the current repository in a fixed environment.
	nox -s build_and_test_on_locked_env

clean:
	rm -f capytaine/green_functions/libs/*.so
	rm -rf build/
	rm -rf dist/
	rm -rf capytaine.egg-info/
	rm -rf docs/_build
	rm -rf .pytest_cache/
	rm -rf .nox/
	rm -rf .venv/
	rm -rf __pycache__ */__pycache__ */*/__pycache__ */*/*/__pycache__
	rm -rf ${HOME}/.cache/capytaine/*

.PHONY: install develop test clean test_fortran_compilation
