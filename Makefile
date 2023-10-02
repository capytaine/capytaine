install:
	pip install .

develop:
	pip install meson-python ninja numpy charset-normalizer # No installed from pyproject.toml in this case...
	pip install --no-build-isolation -e .

test_fortran_compilation:
	meson setup --wipe build_meson && meson compile -C build_meson -j 1

test: develop
	# TODO: use something like nox instead.
	# TODO: Install pytest
	python -m pytest

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

.PHONY: install develop test clean test_fortran_compilation
