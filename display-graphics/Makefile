CC = g++

LDFLAGS += -lglut -lGL -lX11

COWFLAGS += -framework OpenGL -framework GLUT

OBJECTS = display-graphics.o

display: $(OBJECTS)
	$(CC) $(OBJECTS) $(COWFLAGS) -o display

%.o: %.cpp
	$(CC) -c $<

clean:
	rm *.o display
