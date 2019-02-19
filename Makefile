# Development Makefile

DIRECTORY=capytaine/bem

F90FILES = \
$(DIRECTORY)/NemohCore/constants.f90             \
$(DIRECTORY)/NemohCore/Green_Rankine.f90         \
$(DIRECTORY)/NemohCore/Green_wave.f90            \
$(DIRECTORY)/NemohCore/Initialize_Green_wave.f90 \
$(DIRECTORY)/NemohCore/old_Prony_decomposition.f90

SOFILES = \
$(DIRECTORY)/NemohCore.cpython-36m-x86_64-linux-gnu.so

$(SOFILES): $(F90FILES)
	python setup.py build_ext --inplace

update_fortran: $(SOFILES)

develop: $(SOFILES)
	pip install -e .

test: develop
	python -m pytest

clean:
	rm -rf build
	rm -rf capytaine.egg-info/
	rm -f capytaine/**/*.so
	rm -rf __pycache__ */__pycache__ */*/__pycache__

.PHONY: update_fortran develop test clean
