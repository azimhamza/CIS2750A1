CC = clang
CFLAGS = -Wall -std=c99 -pedantic
INCLUDES = /usr/include/python3.7m
LIB = /usr/lib/python3.7/config-3.7m-x86_64-linux-gnu
SWIG = swig

all: libmol.so _molecule.so

T: MolDisplay.py
	python3 MolDisplay.py

MolTest: _molecule.so moltest.py
	python moltest.py

_molecule.so: molecule_wrap.o libmol.so
	$(CC) $(CFLAGS) -shared molecule_wrap.o -L. -lmol -L$(LIB) -lpython3.7m -o _molecule.so

test: test.o libmol.so 
	$(CC) test.o -L. -lmol -lm -o test 

test.o: test.c mol.h
	$(CC) -c $(CFLAGS) test.c 

mol.o: mol.c
	$(CC) -c $(CFLAGS) -fPIC mol.c

libmol.so: mol.o
	$(CC) mol.o -shared -o libmol.so

molecule_wrap.c: molecule.i
	$(SWIG) -python molecule.i

molecule_wrap.o: molecule_wrap.c
	$(CC) -c $(CFLAGS) molecule_wrap.c -I$(INCLUDES) -fPIC -o molecule_wrap.o

clean: 
	rm -f *.o *.so mol molecule.py molecule_wrap.c test

.PHONY: clean swig
