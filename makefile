CC = clang
CFLAGS = -Wall -std=c99 -pedantic
INCLUDES = /usr/include/python3.7m
LIB = /usr/lib/python3.7/config-3.7m-x86_64-linux-gnu

all: libmol.so _molecule.so

swig: molecule.i
	swig -python molecule.i

_molecule.so: molecule_wrap.o libmol.so
	$(CC) $(CFLAGS) -shared molecule_wrap.o -L. -lmol -L$(LIB) -lpython3.7m -o _molecule.so

molecule_wrap.o: swig
	$(CC) -c $(CFLAGS) molecule_wrap.c -I$(INCLUDES) -fPIC -o molecule_wrap.o

libmol.so: mol.o
	$(CC) mol.o -shared -o libmol.so

mol.o: mol.c
	$(CC) -c $(CFLAGS) -fPIC mol.c

test: test.o libmol.so 
	$(CC) test.o -L. -lmol -lm -o test 

test.o: test.c mol.h
	$(CC) -c $(CFLAGS) test.c 

clean: 
	rm -f *.o *.so mol molecule.py molecule_wrap.c test

.PHONY: server
server:
	python3 server.py

.PHONY: run-server
run-server:
	@read -p "Enter port number: " port; python3 server.py localhost $$port &
	sleep 1
	open http://localhost:$$port
