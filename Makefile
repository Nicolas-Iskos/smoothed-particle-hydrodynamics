SPH: run-SPH.o initialize.o integrate.o calculate-field.o smoothing-kernels.o
	nvcc integrate.o initialize.o calculate-field.o smoothing-kernels.o run-SPH.o -o SPH

run-SPH.o: run-SPH.cu particle-data-structures.h simulation-parameters.h
	nvcc -c run-SPH.cu

initialize.o: initialize.cu initialize.h particle-data-structures.h simulation-parameters.h
	nvcc -c initialize.cu

integrate.o: integrate.cu integrate.h particle-data-structures.h simulation-parameters.h
	nvcc -c integrate.cu

calculate-field.o: calculate-field.cu calculate-field.h particle-data-structures.h simulation-parameters.h
	nvcc -c calculate-field.cu

smoothing-kernels.o: smoothing-kernels.cu smoothing-kernels.h particle-data-structures.h simulation-parameters.h
	nvcc -c smoothing-kernels.cu

clean:
	rm *.o SPH
