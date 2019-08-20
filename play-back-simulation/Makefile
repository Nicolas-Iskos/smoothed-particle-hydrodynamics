CC = g++

LDFLAGS += -lglut -lGL -lX11

APPLEGLFLAGS += -framework OpenGL -framework GLUT
LINUXGLFLAGS += -lGL -lGLU -lglut

OBJECTS = play-back-simulation.o

# change APPLEGLFLAGS to LINUXGLFLAGS if compiling on linux
play-back: $(OBJECTS)
	$(CC) $(OBJECTS) $(APPLEGLFLAGS) -o play-back

%.o: %.cpp
	$(CC) -c $<

clean:
	rm *.o play-back
