install:
	pip install .

develop:
	pip install meson-python numpy charset-normalizer # No installed from pyproject.toml in this case...
	pip install --no-build-isolation -e .

test: develop
	# TODO: use something like nox instead.
	# TODO: Install pytest and hypothesis?
	python -m pytest

clean:
	rm -f capytaine/green_functions/libs/*.so
	rm -rf build/
	rm -rf dist/
	rm -rf capytaine.egg-info/
	rm -rf docs/_build
	rm -rf .pytest_cache/
	rm -rf .hypothesis/
	rm -rf __pycache__ */__pycache__ */*/__pycache__ */*/*/__pycache__

.PHONY: install develop test clean
