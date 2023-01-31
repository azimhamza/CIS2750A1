CC=gcc
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=test1.c test2.c test3.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=test1 test2 test3

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $@.o -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *o $(EXECUTABLE)

