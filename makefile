CC = clang
CFLAGS = -Wall -std=c99 -pedantic
INCLUDES = /usr/include/python3.7m
LIB = /usr/lib/python3.7/config-3.7m-x86_64-linux-gnu

all: libmol.so _molecule.so

T: MolDisplay.py
	python3.7 MolDisplay.py

MolTest: _molecule.so moltest.py
	python3.7 moltest.py

swig: molecule.i
	swig -python molecule.i

molecule_wrap.o: swig
	$(CC) -c $(CFLAGS) -c molecule_wrap.c -I$(INCLUDES) -fPIC -o molecule_wrap.o

_molecule.so: molecule_wrap.o libmol.so
	$(CC) $(CFLAGS) -shared molecule_wrap.o -L. -lmol -L$(LIB) -lpython3.7m -o _molecule.so

test: test.o libmol.so
	$(CC) test.o -L. -lmol -lm -o test

test.o: test.c mol.h
	$(CC) -c $(CFLAGS) test.c

mol.o: mol.c
	$(CC) -c $(CFLAGS) -fPIC mol.c

libmol.so: mol.o
	$(CC) -shared -o libmol.so mol.o

run-server:
	python3.7 server.py $(PORT)

clean:
	rm -f *.o *.so mol molecule.py molecule_wrap.c test