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
typedef struct bond
{
unsigned short a1, a2; // index of atoms
unsigned char epairs;
atom *atoms;
double x1, x2, y1, y2, z, len, dx, dy;
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
void bondset( bond *bond, unsigned short *a1, unsigned short *a2, atom **atoms, unsigned char *epairs );
void bondget( bond *bond, unsigned short *a1, unsigned short *a2, atom **atoms, unsigned char *epairs );
molecule *molmalloc(unsigned short atom_max, unsigned short bond_max);
molecule *molcopy(molecule *src);
void molfree(molecule *ptr);
void molappend_atom(molecule *molecule, atom *at);
void molappend_bond(molecule *molecule, bond *bo);
void molsort(molecule *molecule);
int compare_atom(const void *a, const void *b);
int compare_bond(const void *a, const void *b);
void xrotation(xform_matrix xform_matrix, unsigned short deg);
void yrotation(xform_matrix xform_matrix, unsigned short deg);
void zrotation(xform_matrix xform_matrix, unsigned short deg);
void mol_xform(molecule *molecule, xform_matrix matrix);
void compute_coords (bond *bond);

