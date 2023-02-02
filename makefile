all: test1 test2 test3 testPart1.c

test1: test1.o mol.o
	gcc test1.o mol.o -o test1 -lm

test2: test2.o mol.o
	gcc test2.o mol.o -o test2 -lm

test3: test3.o mol.o
	gcc test3.o mol.o -o test3 -lm

testPart1: testPart1.o mol.o
	gcc testPart1.o mol.o -o testPart1 -lm

test1.o: test1.c mol.h
	gcc -c -Wall test1.c -o test1.o

test2.o: test2.c mol.h
	gcc -c -Wall test2.c -o test2.o

test3.o: test3.c mol.h
	gcc -c -Wall test3.c -o test3.o

testPart1.o : testPart1.c mol.h
	gcc -c -Wall testPart1.c -o testPart1.o

mol.o: mol.c mol.h
	gcc -c -Wall mol.c -o mol.o

clean:
	rm -f test1 test2 test3 test1.o test2.o test3.o testPart1.o mol.o
