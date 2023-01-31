
all: test1 test2 test3

test1: test1.o
	gcc -std=c99 -Wall -lm test1.o -o test1

test1.o: test1.c
	gcc -std=c99 -Wall -lm -c test1.c

test2: test2.o
	gcc -std=c99 -Wall -lm test2.o -o test2

test2.o: test2.c
	gcc -std=c99 -Wall -lm -c test2.c

test3: test3.o mol.o
	gcc -std=c99 -Wall -lm test1.o mol.o -o test3

test3.o: test3.c
	gcc -std=c99 -Wall -lm -c test3.c

mol.o: mol.c
	gcc -std=c99 -Wall -lm -c mol.c

clean:
	rm -f *.o test1 test2 test3
