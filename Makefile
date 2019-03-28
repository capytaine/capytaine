# Development Makefile

DIRECTORY=$(PWD)/capytaine/bem/

F90FILES = \
$(DIRECTORY)/NemohCore/constants.f90             \
$(DIRECTORY)/NemohCore/Green_Rankine.f90         \
$(DIRECTORY)/NemohCore/Green_wave.f90            \
$(DIRECTORY)/NemohCore/Initialize_Green_wave.f90 \
$(DIRECTORY)/NemohCore/old_Prony_decomposition.f90

SOFILE = \
$(DIRECTORY)/NemohCore.cpython-36m-x86_64-linux-gnu.so

$(SOFILE): $(F90FILES)
	python setup.py build_ext --inplace

update_fortran: $(SOFILE)

develop: $(SOFILE)
	pip install -e .

test: develop
	python -m pytest

clean:
	rm -f $(SOFILE)
	rm -rf build
	rm -rf capytaine.egg-info/
	rm -rf .pytest_cache/
	rm -rf __pycache__ */__pycache__ */*/__pycache__

pypi: clean
	python setup.py sdist
	python -m twine upload dist/capytaine*.tar.gz

.PHONY: update_fortran develop test clean pypi
