CC = nvcc

CFLAGS = -rdc=true -std c++14

OBJECTS = run-SPH.o particle-data-structures.o integrate.o calculate-field.o \
		  smoothing-kernels.o test-functions.o

COM-DEPS = particle-data-structures.cu particle-data-structures.h \
		   simulation-parameters.h



SPH: $(OBJECTS)
	$(CC) $(OBJECTS) -o SPH

%.o: %.cu
	$(CC) $(CFLAGS) -c $<



run-SPH.o: run-SPH.cu $(COM-DEPS)

particle-data-structures.o: $(COM-DEPS)

integrate.o: integrate.cu integrate.h $(COM-DEPS)

calculate-field.o: calculate-field.cu calculate-field.h $(COM-DEPS)

smoothing-kernels.o: smoothing-kernels.cu smoothing-kernels.h $(COM-DEPS)

test-functions.o: test-functions.cu test-functions.h $(COM-DEPS)

clean:
	rm *.o SPH
