#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


// define the structure of the object
typedef struct atom {
  char element[3];
  double x, y, z;
} atom;

// bond definiton 
typedef struct bond {
  atom *a1, *a2;
  unsigned char epairs;
} bond;

typedef struct molecule {
  unsigned short atom_max, atom_no;
  atom *atoms, **atom_ptrs;
  unsigned short bond_max, bond_no;
  bond *bonds, **bond_ptrs;
} molecule;

// 3-d affine transformation matrix
typedef double xform_matrix[3][3];


void atomset(atom *atom, char element[3], double *x, double *y, double *z);
void atomget(atom *atom, char element[3], double *x, double *y, double *z);
void bondset(bond *bond, atom *a1, atom *a2, unsigned char epairs);
void bondget(bond *bond, atom **a1, atom **a2, unsigned char *epairs);
molecule *molmalloc(unsigned short atom_max, unsigned short bond_max);
molecule *molcopy(molecule *src);
void molfree(molecule *ptr);
void molappend_atom(molecule *molecule, atom *atom);
void molappend_bond(molecule *molecule, bond *bond);
int compare_atom_z(const void *a, const void *b);
int compare_bond_z(const void *a, const void *b);
void molsort(molecule *molecule);
void xrotation(xform_matrix xform_matrix, unsigned short deg);
void yrotation(xform_matrix xform_matrix, unsigned short deg);
void zrotation(xform_matrix xform_matrix, unsigned short deg);

