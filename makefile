CC = clang
CFLAGS = -Wall -std=c99 -pedantic
INCLUDES = /usr/include/python3.7m
LIB = /usr/lib/python3.7/config-3.7m-x86_64-linux-gnu

PORT = 8080

all: libmol.so _molecule.so

T: MolDisplay.py
	python3.7 MolDisplay.py

MolTest: _molecule.so moltest.py
	python3.7 moltest.py

molecule_wrap.c: molecule.i
	swig -python -o molecule_wrap.c molecule.i

molecule_wrap.o: molecule_wrap.c
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

sql:
	python3.7 molsql.py

clean:
	rm -f *.o *.so *.dylib molecule.py molecule_wrap.c test
