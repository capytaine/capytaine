install:
	pip install .

develop:
	pip install -e .[develop]

test: develop
	python -m pytest

clean:
	rm -f capytaine/green_functions/*.so
	rm -rf build/
	rm -rf dist/
	rm -rf capytaine.egg-info/
	rm -rf .pytest_cache/
	rm -rf __pycache__ */__pycache__ */*/__pycache__

pypi: clean
	python setup.py sdist
	python -m twine upload dist/capytaine*.tar.gz

.PHONY: install develop test clean pypi
