# Development Makefile

F90FILES = \
capytaine/NemohCore/Green_Rankine.f90         \
capytaine/NemohCore/Green_wave.f90            \
capytaine/NemohCore/Initialize_Green_wave.f90 \
capytaine/NemohCore/old_Prony_decomposition.f90

SOFILES = \
capytaine/NemohCore.cpython-36m-x86_64-linux-gnu.so

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
	rm -f capytaine/*.so
	rm -rf __pycache__ */__pycache__ */*/__pycache__

.PHONY: update_fortran develop test clean
