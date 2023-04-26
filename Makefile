install:
	pip install .

develop:
	pip install meson-python numpy charset-normalizer # No installed from pyproject.toml in this case...
	pip install --no-build-isolation -e .

test: develop
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

full_archive: clean
	# See https://github.com/capytaine/capytaine/issues/221 for the motivation
	read -p "Version number:" version && tar caf /tmp/capytaine-full-$$version.tar.gz --transform "s,^./,capytaine-$$version/," --exclude '*.tar.gz' --exclude './.git/*' --exclude-vcs-ignores . && mv /tmp/capytaine-full-$$version.tar.gz ./

pypi: clean
	python setup.py sdist
	python -m twine upload dist/capytaine*.tar.gz

.PHONY: install develop test clean full_archive pypi
