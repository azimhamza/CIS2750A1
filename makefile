
all: test1

test1: test1.o mol.o
	gcc -std=c99 -Wall -lm test1.o mol.o -o test1

test1.o: test1.c
	gcc -std=c99 -Wall -lm -c test1.c

mol.o: mol.c
	gcc -std=c99 -Wall -lm -c mol.c

clean:
	rm -f *.o test1