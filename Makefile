CC=nvcc
CXXFLAGS+=-std c++11

SPH: run-SPH.o particle-data-structures.o integrate.o calculate-field.o smoothing-kernels.o test-kernels.o
	$(CC) $(CXXFLAGS) integrate.o particle-data-structures.o calculate-field.o smoothing-kernels.o run-SPH.o test-kernels.o -o SPH

run-SPH.o: run-SPH.cu particle-data-structures.h simulation-parameters.h
	$(CC) $(CXXFLAGS) -c run-SPH.cu

particle-data-structures.o: particle-data-structures.cu initialize.h particle-data-structures.h simulation-parameters.h
	$(CC) $(CXXFLAGS) -c particle-data-structures.cu

integrate.o: integrate.cu integrate.h particle-data-structures.h simulation-parameters.h
	$(CC) $(CXXFLAGS) -c integrate.cu

calculate-field.o: calculate-field.cu calculate-field.h particle-data-structures.h simulation-parameters.h
	$(CC) $(CXXFLAGS) -c calculate-field.cu

smoothing-kernels.o: smoothing-kernels.cu smoothing-kernels.h particle-data-structures.h simulation-parameters.h
	$(CC) $(CXXFLAGS) -c smoothing-kernels.cu

test-kernels.o: test-kernels.h test-kernels.cu particle-data-structures.h simulation-parameters.h
	$(CC) $(CXXFLAGS) -c test-kernels.cu

clean:
	rm *.o SPH
