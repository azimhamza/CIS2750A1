#include <stdlib.h>
#include "mol.h"
#include <math.h>
#include <string.h>

#define PI 3.14159265
// setter for an atom
void atomset(atom *atom, char element[3], double *x, double *y, double *z)
{
    // setting x,y,z
    atom->x = *x;
    atom->y = *y;
    atom->z = *z;
    // copies the items into the addresses.
    strcpy(atom->element, element);
}

// getter for an atom
void atomget(atom *atom, char element[3], double *x, double *y, double *z)
{
    // assigning x,y,z
    *x = atom->x;
    *y = atom->y;
    *z = atom->z;
    // strncpy to copy all three element into the atom structure
    for (int i = 0; i < 3; i++)
    {
        element[i] = atom->element[i];
    }
}

// bond setter
void bondset(bond *bond, atom *a1, atom *a2, unsigned char epairs)
{
    // setting the bonds a1 and a2
    bond->a1 = a1;
    bond->a2 = a2;
    bond->epairs = epairs;
}

// bond getter
void bondget(bond *bond, atom **a1, atom **a2, unsigned char *epairs)
{
    // getting the bonds a1 and a2
    *a1 = bond->a1;
    *a2 = bond->a2;
    *epairs = bond->epairs;
}

molecule *molmalloc(unsigned short atom_max, unsigned short bond_max)
{
    molecule *m = (molecule *)malloc(sizeof(molecule));
    if (m == NULL)
    {
        return NULL;
    }

    m->atom_max = atom_max;
    m->atom_no = 0;
    m->atoms = (atom *)malloc(atom_max * sizeof(atom));
    if (m->atoms == NULL)
    {
        perror("null atoms in molmalloc <Error>");
        return NULL;
    }

    m->atom_ptrs = (atom **)malloc(atom_max * sizeof(atom *));
    if (m->atom_ptrs == NULL)
    {
        perror("null atom_ptrs in molmalloc <Error>");
        return NULL;
    }

    m->bond_max = bond_max;
    m->bond_no = 0;
    m->bonds = (bond *)malloc(bond_max * sizeof(bond));
    if (m->bonds == NULL)
    {
        perror("null bonds in molmalloc <Error>");
        return NULL;
    }

    m->bond_ptrs = (bond **)malloc(bond_max * sizeof(bond *));
    if (m->bond_ptrs == NULL)
    {
        perror("null bond_ptrs in molmalloc <Error>");
        return NULL;
    }

    return m;
}

molecule *molcopy(molecule *src)
{
    molecule *dst = (molecule *)malloc(sizeof(molecule));
    dst->atom_max = src->atom_max;
    dst->atom_no = src->atom_no;
    dst->bond_max = src->bond_max;
    dst->bond_no = src->bond_no;
    dst->atoms = (atom *)malloc(src->atom_max * sizeof(atom));
    dst->atom_ptrs = (atom **)malloc(src->atom_max * sizeof(atom *));
    dst->bonds = (bond *)malloc(src->bond_max * sizeof(bond));
    dst->bond_ptrs = (bond **)malloc(src->bond_max * sizeof(bond *));
    for (int i = 0; i < src->atom_no; i++)
    {
        molappend_atom(dst, &src->atoms[i]);
    }
    for (int i = 0; i < src->bond_no; i++)
    {
        molappend_bond(dst, &src->bonds[i]);
    }
    return dst;
}

// free memory
void molfree(molecule *ptr)
{
    /*
     * This function frees the memory associated with the molecule pointed to by ptr.
     * This includes the arrays atoms, atom_ptrs, bonds, and bond_ptrs.
     */

    // Error handler for when the pointer is null
    if (ptr == NULL)
    {
        printf("<Error> molecule struct is empty\n");
        return;
    }
    // Error handler for when the atoms and atom_ptrs arrays are null
    if (ptr->atoms == NULL || ptr->atom_ptrs == NULL)
    {
        printf("<Error> molecule struct is not initialized\n");
        return;
    }
    // Freeing the memory for the atoms array
    free(ptr->atoms);
    // Freeing the memory for the atom_ptrs array
    free(ptr->atom_ptrs);
    // Freeing the memory for the bond array
    free(ptr->bonds);
    // Freeing the memory for the bond_ptrs array
    free(ptr->bond_ptrs);
    // Freeing the memory for the molecule struct
    free(ptr);
}




/*void molappend_atom( molecule *molecule, atom *atom );
This function should copy the data pointed to by atom to the first “empty” atom in atoms in the
molecule pointed to by molecule, and set the first “empty” pointer in atom_ptrs to the same
atom in the atoms array incrementing the value of atom_no. If atom_no equals atom_max, then
atom_max must be incremented, and the capacity of the atoms, and atom_ptrs arrays
increased accordingly. If atom_max was 0, it should be incremented to 1, otherwise it should be
doubled. Increasing the capacity of atoms, and atom_ptrs should be done using realloc so
that a larger amount of memory is allocated and the existing data is copied to the new location.
IMPORTANT: After mallocing or reallocing enough memory for atom_ptrs, these pointers
should be made to point to the corresponding atoms in the new atoms array (not the old array
which may have been freed)*/

// adds it to the end
void molappend_atom( molecule *molecule, atom *at ) {
   int new_capacity;
   int i;

   // Check if the current atom_no equals atom_max
   if(molecule->atom_no == molecule->atom_max) {
      // If atom_max is 0, set it to 1
      if(molecule->atom_max == 0) {
         molecule->atom_max = 1;
      }
      // Else, double atom_max
      else {
         molecule->atom_max *= 2;
      }
      // Increase the capacity of atoms and atom_ptrs arrays
      new_capacity = molecule->atom_max * sizeof(atom);
      molecule->atoms = realloc(molecule->atoms, new_capacity);
      new_capacity = molecule->atom_max * sizeof(atom *);
      molecule->atom_ptrs = realloc(molecule->atom_ptrs, new_capacity);

      // Update the pointers in atom_ptrs to point to the corresponding atoms in the new atoms array
      for(i = 0; i < molecule->atom_no; i++) {
         molecule->atom_ptrs[i] = &(molecule->atoms[i]);
      }
   }
   // Copy the data pointed to by atom to the first "empty" atom in atoms in the molecule pointed to by molecule
   molecule->atoms[molecule->atom_no] = *at;
   molecule->atom_ptrs[molecule->atom_no] = &(molecule->atoms[molecule->atom_no]);
   // Increment the value of atom_no
   molecule->atom_no++;
}



void molappend_bond(molecule *molecule, bond *bo) {
  if (molecule->bond_no == molecule->bond_max) {
    if (molecule->bond_max == 0) {
      molecule->bond_max = 1;
    } else {
      molecule->bond_max *= 2;
    }
    molecule->bonds = realloc(molecule->bonds,
                              molecule->bond_max * sizeof(bond));
    molecule->bond_ptrs = realloc(molecule->bond_ptrs,
                                  molecule->bond_max * sizeof(bond*));
  }
  molecule->bonds[molecule->bond_no] = *bo;
  molecule->bond_ptrs[molecule->bond_no] = &(molecule->bonds[molecule->bond_no]);
  molecule->bond_no++;
}


void molsort(molecule *molecule)
{
    int (*compare)(const void *a, const void *b);
    compare = compare_atom_z;
    qsort(molecule->atom_ptrs, molecule->atom_no, sizeof(atom *), compare);

    compare = compare_bond_z;
    qsort(molecule->bond_ptrs, molecule->bond_no, sizeof(bond *), compare);
}

/* helper function for molsort to compare z values of atoms */
int compare_atom_z(const void *a, const void *b)
{
    atom *aa = *(atom **)a;
    atom *bb = *(atom **)b;
    double za = aa->z;
    double zb = bb->z;
    if (za < zb)
        return -1;
    else if (za > zb)
        return 1;
    return 0;
}

/* helper function for molsort to compare z values of bonds */
int compare_bond_z(const void *a, const void *b)
{
    bond *ab = *(bond **)a;
    bond *bb = *(bond **)b;
    double za = (ab->a1->z + ab->a2->z) / 2.0;
    double zb = (bb->a1->z + bb->a2->z) / 2.0;
    if (za < zb)
        return -1;
    else if (za > zb)
        return 1;
    return 0;
}

// Function to apply a rotation transformation in the X-axis to a given transformation matrix
// deg - angle of rotation in degrees
// xform_matrix - transformation matrix to be updated
void xrotation(xform_matrix xform_matrix, unsigned short deg)
{
    double rad = deg * PI / 180.0; // Convert degrees to radians

    xform_matrix[0][0] = 1;
    xform_matrix[0][1] = 0;
    xform_matrix[0][2] = 0;
    xform_matrix[1][0] = 0;
    xform_matrix[1][1] = cos(rad);
    xform_matrix[1][2] = -sin(rad);
    xform_matrix[2][0] = 0;
    xform_matrix[2][1] = sin(rad);
    xform_matrix[2][2] = cos(rad);
}

// This function calculates the transformation matrix for y-axis rotation
// The input is a 3x3 transformation matrix and an angle in degrees
// The output is the transformation matrix with y-axis rotation applied
void yrotation(xform_matrix xform_matrix, unsigned short deg)
{
    double rad = deg * PI / 180; // Convert degrees to radians

    xform_matrix[0][0] = cos(rad);
    xform_matrix[0][1] = 0;
    xform_matrix[0][2] = sin(rad);
    xform_matrix[1][0] = 0;
    xform_matrix[1][1] = 1;
    xform_matrix[1][2] = 0;
    xform_matrix[2][0] = -sin(rad);
    xform_matrix[2][1] = 0;
    xform_matrix[2][2] = cos(rad);
}

void zrotation(xform_matrix xform_matrix, unsigned short deg)
{
    // Convert the rotation angle from degrees to radians
    double rad = deg * PI / 180;

    // Initialize the transformation matrix to the identity matrix
    xform_matrix[0][0] = cos(rad);
    xform_matrix[0][1] = -sin(rad);
    xform_matrix[0][2] = 0;
    xform_matrix[1][0] = sin(rad);
    xform_matrix[1][1] = cos(rad);
    xform_matrix[1][2] = 0;
    xform_matrix[2][0] = 0;
    xform_matrix[2][1] = 0;
    xform_matrix[2][2] = 1;
}

// This helper function computes the transformation matrix for rotation about an arbitrary axis
// Input is the 3x3 transformation matrix and an angle in radians
// Output is the transformation matrix with rotation applied

// This function computes the transformation matrix for rotation about the z-axis
// Input is the 3x3 transformation matrix and an angle in degrees
// Output is the transformation matrix with z-axis rotation applied

// Convert the input angle from degrees to radians

// This function transforms an atom by applying the specified transformation matrix
// Input is a pointer to an atom and a 3x3 transformation matrix
// Output is the transformed atom
void mol_xform(molecule *molecule, xform_matrix matrix)
{
    for (int i = 0; i < molecule->atoms; i++)
    {
        atom *transform_atom = &molecule->atoms[i];
        double x = transform_atom->x;
        double y = transform_atom->y;
        double z = transform_atom->z;
        transform_atom->x = x * matrix[0][0] + y * matrix[0][1] + z * matrix[0][2];
        transform_atom->y = x * matrix[1][0] + y * matrix[1][1] + z * matrix[1][2];
        transform_atom->z = x * matrix[2][0] + y * matrix[2][1] + z * matrix[2][2];
    }
}
