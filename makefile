CC = clang
CFLAGS = -Wall -std=c99 -pedantic
INCLUDES = /usr/include/python3.7m
LIB = /usr/lib/python3.7/config-3.7m-x86_64-linux-gnu
PORT = 8080

all: T MolTest

T: molecule.py MolDisplay.py
	python3.7m MolDisplay.py

MolTest: molecule.py _molecule.so moltest.py
	python3.7m moltest.py

molecule.py: molecule.i mol.c mol.h
	swig -python -o molecule_wrap.c molecule.i
	$(CC) -c -fpic molecule_wrap.c -I$(INCLUDES)
	$(CC) -c -fpic mol.c -I$(INCLUDES)
	$(CC) -shared molecule_wrap.o mol.o -o _molecule.so -L$(LIB) -lpython3.7m

_molecule.so: molecule.py mol.o
	$(CC) -shared -o _molecule.so mol.o

test: test.o libmol.so
	$(CC) test.o -L. -lmol -lm -o test

test.o: test.c mol.h
	$(CC) -c $(CFLAGS) test.c

mol.o: mol.c
	$(CC) -c $(CFLAGS) -fPIC mol.c

libmol.so: mol.o
	$(CC) -shared -o libmol.so mol.o

# run-server:
# 	python3.11 server.py $(PORT)

sql:
	python3.7 molsql.py

clean:
	rm -f *.o *.so molecule.py molecule_wrap.c test
