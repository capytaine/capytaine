### todo:
###   - deps between src files so that we can run make -j ?


### ==============================================
### dirs and flags

NAME = delhommeau
SRCDIR = src
OBJDIR = obj
LIBDIR = lib

COMP = gfortran
FLAGS = -Wall -fimplicit-none -ffree-line-length-0 -cpp -fopenmp
DEBUG = -O0 -g -Wextra -Wconversion -fbacktrace -fcheck=all
OPT = -O3 -march=native
PROF = -pg

### !! choose debug or optim
# FLAGS += $(DEBUG)
FLAGS += $(OPT)
#FLAGS += $(PROF)


### ==============================================
### source selection

SRC_S=src/float32.f90\
      src/constants.f90\
      src/Green_Rankine.f90\
      src/Delhommeau_integrals.f90\
      src/Green_wave.f90\
      src/matrices.f90\
      src/old_Prony_decomposition.f90

SRC_D=src/float64.f90\
      src/constants.f90\
      src/Green_Rankine.f90\
      src/Delhommeau_integrals.f90\
      src/Green_wave.f90\
      src/matrices.f90\
      src/old_Prony_decomposition.f90

### !! choose single or double comp.
# SRC=$(SRC_S)
SRC=$(SRC_D)

### objs
OBJ=$(SRC:.f90=.o)
OBJ:=$(addprefix $(OBJDIR)/,$(OBJ))
LIB_STATIC=$(LIBDIR)/lib$(NAME).a
LIB_DYNAMIC=$(LIBDIR)/lib$(NAME).so


### ==============================================
### example

EXAMPLES_SRC=examples/minimal/minimal_example.f90
EXAMPLES_BIN=$(EXAMPLES_SRC:.f90=.bin)


### ==============================================
### benchmarks

BENCH_SRC=benchmarks/openmp/benchmark_omp.f90\
	  benchmarks/tabulations/benchmark_tabulation.f90\
	  benchmarks/profiling/benchmark_profiling.f90
BENCH_BIN=$(BENCH_SRC:.f90=.bin)


### ==============================================
### targets

.PHONY: .init clean mrproper

$(OBJDIR)/%.o: %.f90
	$(COMP) $(FLAGS) -J$(OBJDIR) -c $< -o $@

%.bin: %.f90 $(OBJ)
	$(COMP) $(FLAGS) -J$(OBJDIR) -L$(LIBDIR) -l$(NAME) $^ -o $@

.init:
	@mkdir -p $(OBJDIR)/$(SRCDIR)
	@mkdir -p $(LIBDIR)

lib: .init $(LIB_STATIC) $(LIB_DYNAMIC)

$(LIB_STATIC): $(OBJ)
	ar r $@ $^

$(LIB_DYNAMIC): $(OBJ)
	$(COMP) -shared -fPIC $(FLAGS) -o $@ $^

example: lib $(EXAMPLES_BIN)

benchmark: lib $(BENCH_BIN)

run_min_ex: example
	@$(EXAMPLES_BIN)
	@echo ''
	@python3 examples/minimal/minimal_example.py

run_bench_tab: benchmark
	@benchmarks/tabulations/benchmark_tabulation.bin

run_bench_omp: benchmark
	@benchmarks/openmp/benchmark_omp.bin
	@mkdir -p results
	@mv benchmark_omp.csv results/benchmark_omp_$(shell git rev-parse --short HEAD)_$(shell git branch --show-current).csv
	@python3 benchmarks/openmp/read_output.py

# requires gprof2dot
run_bench_prof: benchmark
	@export OMP_NUM_THREADS=1 && benchmarks/profiling/benchmark_profiling.bin
	@gprof benchmarks/profiling/benchmark_profiling.bin gmon.out > prof.txt
	@gprof2dot < prof.txt > prof.dot && dot prof.dot -Tpng -o prof.png
	@rm gmon.out prof.txt prof.dot

clean:
	@find . -name '*~' -exec rm {} \;
	@rm -rf $(OBJ) $(OBJDIR) $(EXAMPLES_BIN) $(BENCH_BIN)

mrproper: clean
	@rm -rf $(LIBDIR)
