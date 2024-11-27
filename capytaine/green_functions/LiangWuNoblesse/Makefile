BUILD_DIR=build

run_test: $(BUILD_DIR)/test_program
	$(BUILD_DIR)/test_program

$(BUILD_DIR)/test_program: test_program.f90 $(BUILD_DIR)/LiangWuNoblesseWaveTerm.o
	mkdir -p $(BUILD_DIR)
	gfortran -fopenmp $^ -o $@ -J$(BUILD_DIR)

$(BUILD_DIR)/LiangWuNoblesseWaveTerm.o: LiangWuNoblesseWaveTerm.f90
	mkdir -p $(BUILD_DIR)
	gfortran -c $< -o $@ -J$(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)
.PHONY: run_test clean


