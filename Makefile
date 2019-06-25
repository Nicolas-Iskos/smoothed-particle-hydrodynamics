SPH: run-SPH.o integrate.o calculate-field.o smoothing-kernels.o
	nvcc integrate.o calculate-field.o smoothing-kernels.o run-SPH.o -o SPH

run-SPH.o: run-SPH.cu run-SPH.h
	nvcc -c run-SPH.cu

integrate.o: integrate.cu integrate.h
	nvcc -c integrate.cu

calculate-field.o: calculate-field.cu calculate-field.h
	nvcc -c calculate-field.cu

smoothing-kernels.o: smoothing-kernels.cu smoothing-kernels.h
	nvcc -c smoothing-kernels.cu

clean:
	rm *.o SPH
