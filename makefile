CC = clang
CFLAGS = -Wall -std=c99 -pedantic
INCLUDES = /Library/Frameworks/Python.framework/Versions/3.11/include/python3.11
LIB = /Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/config-3.11-darwin

PORT = 8080

all: libmol.dylib _molecule.so

T: MolDisplay.py
	python3.7 MolDisplay.py

MolTest: _molecule.so moltest.py
	python3.7 moltest.py

molecule_wrap.c: molecule.i
	swig -python -o molecule_wrap.c molecule.i

molecule_wrap.o: molecule_wrap.c
	$(CC) -c $(CFLAGS) -c molecule_wrap.c -I$(INCLUDES) -fPIC -o molecule_wrap.o

_molecule.so: molecule_wrap.o libmol.dylib
	$(CC) $(CFLAGS) -shared molecule_wrap.o -L. -lmol -L$(LIB) -lpython3.11 -o _molecule.so

test: test.o libmol.dylib
	$(CC) test.o -L. -lmol -lm -o test

test.o: test.c mol.h
	$(CC) -c $(CFLAGS) test.c

mol.o: mol.c
	$(CC) -c $(CFLAGS) -fPIC mol.c

libmol.dylib: mol.o
	$(CC) -dynamiclib -o libmol.dylib mol.o

run-server:
	python3.11 server.py $(PORT)\

sql:
	python3.11 molsql.py

clean:
	rm -f *.o *.so *.dylib molecule.py molecule_wrap.c test
