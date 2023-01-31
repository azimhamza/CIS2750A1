
all: test2

test2: test2.o mol.o
	gcc -std=c99 -Wall -lm test2.o mol.o -o test2

test2.o: test2.c
	gcc -std=c99 -Wall -lm -c test2.c

mol.o: mol.c
	gcc -std=c99 -Wall -lm -c mol.c

clean:
	rm -f *.o test2