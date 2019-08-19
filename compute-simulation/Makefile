CC = nvcc

CFLAGS = -rdc=true -std c++14

OBJECTS = compute-simulation.o particle-data-structures.o integrate.o \
		  calculate-field.o smoothing-kernels.o


COM-DEPS = particle-data-structures.cu particle-data-structures.h



compute: $(OBJECTS)
	$(CC) $(OBJECTS) -o compute


%.o: %.cu
	$(CC) $(CFLAGS) -c $<

compute-simulation.o: compute-simulation.cu $(COM-DEPS)

particle-data-structures.o: $(COM-DEPS)

integrate.o: integrate.cu integrate.h $(COM-DEPS)

calculate-field.o: calculate-field.cu calculate-field.h $(COM-DEPS)

smoothing-kernels.o: smoothing-kernels.cu smoothing-kernels.h $(COM-DEPS)



clean:
	rm *.o compute
