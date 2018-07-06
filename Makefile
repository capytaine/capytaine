# Development Makefile

F90FILES = \
capytaine/NemohCore/Green_1.f90 \
capytaine/NemohCore/Green_2.f90 \
capytaine/NemohCore/Initialize_Green_2.f90 \
capytaine/NemohCore/old_Prony_decomposition.f90 \
capytaine/NemohCore/Wavenumber.f90

SOFILES = \
capytaine/_Green.cpython-36m-x86_64-linux-gnu.so \
capytaine/_Wavenumber.cpython-36m-x86_64-linux-gnu.so

$(SOFILES): $(F90FILES)
	python setup.py build_ext --inplace

update_fortran: $(SOFILES)

develop: $(SOFILES)
	pip install -e .

test: develop
	python -m pytest

clean:
	rm -rf build
	rm -f capytaine/*.so
	rm -rf */__pycache__

.PHONY: update_fortran develop test clean
