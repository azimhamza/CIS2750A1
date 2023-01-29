CC=gcc
CFLAGS=-c -Wall

all: test1

test1: test1.o mol.o
	$(CC) test1.o mol.o -o test1

test1.o: test1.c
	$(CC) $(CFLAGS) test1.c

mol.o: mol.c
	$(CC) $(CFLAGS) mol.c

clean:
	rm -f *.o test1
