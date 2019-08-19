CC = g++

LDFLAGS += -lglut -lGL -lX11

COWFLAGS += -framework OpenGL -framework GLUT

OBJECTS = play-back-simulation.o

play-back: $(OBJECTS)
	$(CC) $(OBJECTS) $(COWFLAGS) -o play-back

%.o: %.cpp
	$(CC) -c $<

clean:
	rm *.o play-back
