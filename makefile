CC = gcc
CFLAGS = -Wall -std=c99 -pedantic
INCLUDES = /Library/Frameworks/Python.framework/Versions/3.11/include/python3.11
LIB = /Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/config-3.11-darwin
PORT = 8080

all: T MolTest

T: molecule.py MolDisplay.py
	python3.11 MolDisplay.py

MolTest: molecule.py _molecule.so moltest.py
	python3.11 moltest.py

molecule.py: molecule.i mol.c mol.h
	swig -python -o molecule_wrap.c molecule.i
	$(CC) -c -fpic molecule_wrap.c -I$(INCLUDES)
	$(CC) -c -fpic mol.c -I$(INCLUDES)
	$(CC) -shared molecule_wrap.o mol.o -o _molecule.so -L$(LIB) -lpython3.11

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
	python3.11 molsql.py

clean:
	rm -f *.o *.so molecule.py molecule_wrap.c test
