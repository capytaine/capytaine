# Development Makefile

DIRECTORY=$(PWD)/capytaine/green_functions

F90FILES = \
$(DIRECTORY)/Delhommeau_f90/constants.f90             \
$(DIRECTORY)/Delhommeau_f90/matrices.f90             \
$(DIRECTORY)/Delhommeau_f90/Green_Rankine.f90         \
$(DIRECTORY)/Delhommeau_f90/Green_wave.f90            \
$(DIRECTORY)/Delhommeau_f90/Initialize_Green_wave.f90 \
$(DIRECTORY)/Delhommeau_f90/old_Prony_decomposition.f90

SOFILE = \
$(DIRECTORY)/Delhommeau_f90.cpython-36m-x86_64-linux-gnu.so

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
