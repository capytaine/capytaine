test_program: test_program.f90 Module_GreenFuncGlobal.o
	gfortran -fopenmp $^ -o $@

Module_GreenFuncGlobal.o: Module_GreenFuncGlobal.f90
	gfortran -c $< -o $@


