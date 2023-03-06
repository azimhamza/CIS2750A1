CC = clang
CFLAGS = -Wall -std=c99 -pedantic
INCLUDES = /Library/Frameworks/Python.framework/Versions/3.11/include/python3.11 
INCLUDES = /usr/include/python3.7
LIB = /usr/lib/python3.7/config-3.7m-x86_64-linux-gnu
export LDFLAGS

all: _molecule.so

_molecule.so: molecule_wrap.o libmol.so
	$(CC) $(CFLAGS) -shared molecule_wrap.o -L. -lmol $(LDFLAGS) -lpython3.11 -dynamiclib -o _molecule.so

molecule_wrap.o: molecule_wrap.c
	$(CC) -c $(CFLAGS) molecule_wrap.c -I$(INCLUDES) -fPIC -o molecule_wrap.o

libmol.so: mol.o
	$(CC) mol.o -shared -o libmol.so

mol.o: mol.c
	$(CC) -c $(CFLAGS) -fPIC mol.c

clean:
	rm -f *.o *.so