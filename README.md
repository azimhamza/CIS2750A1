# CIS2750A1

**Introduction**

This repository contains the source code for a molecule manipulation library in C. The library provides functions to define atoms and bonds, create and modify molecules, perform 3-d affine transformations, and sort atoms and bonds.

**Dependencies**

**GCC compiler**
make utility
**Usage**
To build the library and test programs, run:
- make all
This will produce three test executables: test1, test2, and test3.

**To clean up the build files, run:**
make clean

**File Structure**
mol.c: source code for the molecule manipulation library
mol.h: header file for the molecule manipulation library
test1.c, test2.c, test3.c: source code for the test programs
Makefile: make file to build the library and test programs

**Library Functions**
- atomset(atom *atom, char element[3], double *x, double *y, double *z): **sets the element, x, y, and z coordinates of an atom**
- atomget(atom *atom, char element[3], double *x, double *y, double *z): **retrieves the element, x, y, and z coordinates of an atom**
- bondset(bond *bond, atom *a1, atom *a2, unsigned char epairs): **sets the atom pointers and electron pair count of a bond**
- bondget(bond *bond, atom **a1, atom **a2, unsigned char *epairs): **retrieves the atom pointers and electron pair count of a bond**
- molmalloc(unsigned short atom_max, unsigned short bond_max): **creates and returns a new molecule with specified maximum number of atoms and bonds**
- molcopy(molecule *src): **creates and returns a copy of a molecule**
- molfree(molecule *ptr): **frees a molecule and its associated memory**
- molappend_atom(molecule *molecule, atom *at): **appends an atom to a molecule**
- molappend_bond(molecule *molecule, bond *bo): **appends a bond to a molecule**
- molsort(molecule *molecule): **sorts the atoms and bonds of a molecule by z coordinate**
- compare_atom_z(const void *a, const void *b): **comparison function for sorting atoms by z coordinate**
- compare_bond_z(const void *a, const void *b): **comparison function for sorting bonds by z coordinate**
- xrotation(xform_matrix xform_matrix, unsigned short deg): **performs an x-axis rotation transformation**
- yrotation(xform_matrix xform_matrix, unsigned short deg): **performs a y-axis rotation transformation**
- zrotation(xform_matrix xform_matrix, unsigned short deg): **performs a z-axis rotation transformation**
- mol_xform(molecule *molecule, xform_matrix matrix): **applies a transformation matrix to a molecule.**
