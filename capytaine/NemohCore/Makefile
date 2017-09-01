all:
	f2py -c -m _Green Green_1.f90 Initialize_Green_2.f90 Green_2.f90 Build_matrices.f90
	f2py -c -m _Wavenumber Wavenumber.f90

clean:
	rm -f *.o *.so *.mod
	rm -rf __pycache__
