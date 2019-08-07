CC = g++

LDFLAGS += -lGL -lGLU -lglut

OBJECTS = display-graphics.o

display: $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o display

%.o: %.cpp
	$(CC) -c $<

clean:
	rm *.o display
