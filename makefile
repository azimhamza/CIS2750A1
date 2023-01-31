
all: test3

test3: test3.o mol.o
	gcc -std=c99 -Wall -lm test1.o mol.o -o test3

test3.o: test3.c
	gcc -std=c99 -Wall -lm -c test3.c

mol.o: mol.c
	gcc -std=c99 -Wall -lm -c mol.c

clean:
	rm -f *.o test3