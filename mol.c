#include <stdlib.h>
#include "mol.h"
#include <math.h>
#include <string.h>

// define PI as linux server has a bug with math.h and PI
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
    strcpy(element, atom->element);
    // assigning x,y,z
    *x = atom->x;
    *y = atom->y;
    *z = atom->z;
    // strncpy to copy all three element into the atom structure
}

// updated bond setter
void bondset(bond *bond, unsigned short *a1, unsigned short *a2, atom **atoms, unsigned char *epairs)
{
    // setting the bonds a1 and a2  inside atoms and epairs
    bond->a1 = atoms[*a1];
    bond->a2 = atoms[*a2];
    bond->atoms = *atoms;
    bond->epairs = epairs;
    compute_coords(bond);
}

// bond getter
void bondget(bond *bond, unsigned short *a1, unsigned short *a2, atom **atoms, unsigned char *epairs)
{
    // getting the bonds a1 and a2
    *a1 = bond->a1;
    *a2 = bond->a2;
    *epairs = bond->epairs;
    *atoms = bond->atoms;
}

molecule *molmalloc(unsigned short atom_max, unsigned short bond_max)
{
    // allocating memory for the molecule structure and checking for errors
    molecule *ptr = (molecule *)malloc(sizeof(molecule));
    if (ptr == NULL)
    {
        perror("null molecule in molmalloc <Error>");
        return NULL;
    }
    else
    {
        ptr->atom_max = atom_max;
        ptr->atom_no = 0;
    }

    if (ptr->atom_max == 0)
    {
        ptr->atoms = NULL;
        ptr->atom_ptrs = NULL;
    }
    else
    {
        ptr->atoms = (struct atom *)malloc(sizeof(struct atom) * atom_max); // allocate memory for atoms
        ptr->atom_ptrs = (struct atom **)malloc(sizeof(struct atom *) * atom_max);
    }

    // allocating memory for the bonds and bond pointers and checking for errors
    ptr->bond_max = bond_max;
    ptr->bonds = 0;

    if (ptr->bond_max == 0)
    {
        ptr->bond_ptrs = NULL;
        ptr->bonds = NULL;
    }
    else
    {
        ptr->bond_ptrs = (bond **)malloc(bond_max * sizeof(bond *));
        ptr->bonds = (struct bond *)malloc(sizeof(struct bond) * bond_max);

        if (ptr->bond_ptrs == NULL || ptr->bonds == NULL)
        {
            perror("memory allocation error in molmalloc <Error>");
            free(ptr);
            return NULL;
        }
    }

    // return the molecule
    return ptr;
}

// copy the molecule to a new molecule and return the new molecule
molecule *molcopy(molecule *src)
{

    molecule *holder;
    holder = molmalloc(src->atom_max, src->bond_max);

    int i = 0;
    while (i < src->atom_no)
    {
        atom *holder2 = src->atoms + i;
        molappend_atom(holder, holder2);
        i++;
    }

    while (i < src->bond_no)
    {
        bond *holder3 = src->bonds + i;
        molappend_bond(holder, holder3);
        i++;
    }
    holder->atom_no = src->atom_no;
    holder->bond_no = src->bond_no;
    return holder;
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
    if ( ptr->bonds == NULL || ptr->bond_ptrs == NULL)
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

// append a bond to the molecule structure pointed to by molecule and return the bond pointer
void molappend_atom(molecule *molecule, atom *at)
{
  if (molecule->atom_max == 0) // if atom_max is 0 make it one and malloc space
    {

        molecule->atom_max = 1;
        molecule->atoms = (struct atom *)malloc(sizeof(struct atom) * molecule->atom_max);
        molecule->atom_ptrs = (struct atom **)malloc(sizeof(struct atom) * molecule->atom_max);
    }
    else if (molecule->atom_max == molecule->atom_no) // if atom_max is equal to atom_no realloc space
    {

        molecule->atom_max = molecule->atom_max * 2;
        molecule->atoms = (struct atom *)realloc(molecule->atoms, sizeof(struct atom) * molecule->atom_max);
        molecule->atom_ptrs = (struct atom **)realloc(molecule->atom_ptrs, sizeof(struct atom *) * molecule->atom_max);
        for (int i = 0; i < molecule->atom_no; i++)
        { // ensure atom_ptrs are pointed to the correct new atom!
            molecule->atom_ptrs[i] = molecule->atoms + (i);
        }
    }

    molecule->atoms[molecule->atom_no] = *at;
    molecule->atom_ptrs[molecule->atom_no] = molecule->atoms + (molecule->atom_no); // addres of pointer given to doubel pointer
    molecule->atom_no++;
}

void molappend_bond(molecule *molecule, bond *bo)
{
  
     if (molecule->bond_max == 0) // if atom_bond is 0 make it one and malloc space
    {
        molecule->bond_max = 1;
        molecule->bonds = (struct bond *)malloc(sizeof(struct bond) * molecule->bond_max);
        molecule->bond_ptrs = (struct bond **)malloc(sizeof(struct bond *) * molecule->bond_max);
    }
    if (molecule->bond_max == molecule->bond_no) // if bond_max is equal to atom_no realloc space
    {
        molecule->bond_max = molecule->bond_max * 2;
        molecule->bonds = (struct bond *)realloc(molecule->bonds, sizeof(struct bond) * molecule->bond_max);
        molecule->bond_ptrs = (struct bond **)realloc(molecule->bond_ptrs, sizeof(struct bond *) * molecule->bond_max);
        for (int i = 0; i < molecule->bond_no; i++)
        { // ensure bond_ptrs are pointed to the correct new atom!
            molecule->bond_ptrs[i] = molecule->bonds + (i);
        }
    }

    molecule->bonds[molecule->bond_no] = *bo;
    molecule->bond_ptrs[molecule->bond_no] = molecule->bonds + (molecule->bond_no); // addres of pointer given to doubel pointer
    molecule->bond_no++;
}

// void print_molecule function prints the molecule structure pointed to by molecule to the file pointed to by fp.
int bond_comp(const void *a, const void *b)
{
    bond *bond1 = *(bond **)a;
    bond *bond2 = *(bond **)b;
    // compare z
    if (bond1->z < bond2->z)
        return -1;
    else if (bond1->z > bond2->z)
        return 1;
    else
        return 0;
}

// void computer_coords function computes the z, x1, y1, x2, y2, len, dx, and dy values of the bond and set
// them in the appropriate structure member variables.

void compute_coords(bond *bond)
{

    bond->x1 = bond->atoms[bond->a1].x;
    bond->x2 = bond->atoms[bond->a2].x;

    bond->y1 = bond->atoms[bond->a1].y;
    bond->y2 = bond->atoms[bond->a2].y;

    // computing z as the average distance between the two atoms
    bond->z = ((bond->atoms[bond->a2].z + bond->atoms[bond->a1].z) / 2); // average of 2 z values

    // computing the length of the bond
    bond->len = sqrt(pow(bond->x1 - bond->x2, 2) + pow(bond->y1 - bond->y2, 2));
    // computing dx and dy
    bond->dx = (bond->x2 - bond->x1) / bond->len;
    bond->dy = (bond->y2 - bond->y1) / bond->len;
}

// Comparison functions for qsort
int compare_atom(const void *a, const void *b)
{
    const atom *atom_a = *(struct atom **)a;
    const atom *atom_b = *(struct atom **)b;
    if (atom_a->z < atom_b->z)
    {
        return -1;
    }
    else if (atom_a->z > atom_b->z)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int compare_bond(const void *a, const void *b)
{
    const bond *bond_a = *(struct bond **)a;
    const bond *bond_b = *(struct bond **)b;
    if (bond_a->z < bond_b->z)
    {
        return -1;
    }
    else if (bond_a->z > bond_b->z)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void molsort(molecule *molecule)
{
    // Sort the atoms by z-coordinate
    qsort(molecule->atoms, molecule->atom_no, sizeof(atom), compare_atom);
    // Sort the bonds by z-coordinate
    qsort(molecule->bond_ptrs, molecule->bond_no, sizeof(bond *), compare_bond);
}

// Function to apply a rotation transformation in the X-axis to a given transformation matrix
// deg - angle of rotation in degrees
// xform_matrix - transformation matrix to be updated
void xrotation(xform_matrix xform_matrix, unsigned short deg)
{
    // Convert degrees to radians
    double rad = deg * PI / 180.0; // Convert degrees to radians

    // Initialize the transformation matrix to the identity matrix
    // Apply the rotation transformation to the transformation matrix

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
    // Convert the rotation angle from degrees to radians
    double rad = deg * PI / 180; // Convert degrees to radians
    // Initialize the transformation matrix to the identity matrix
    // Apply the rotation transformation to the transformation matrix
    // The rotation transformation is a 3x3 matrix

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

// This function computes the transformation matrix for rotation about the z-axis
// Input is the 3x3 transformation matrix and an angle in degrees
// Output is the transformation matrix with z-axis rotation applied

void zrotation(xform_matrix xform_matrix, unsigned short deg)
{
    // Convert the rotation angle from degrees to radians
    double rad = deg * PI / 180;

    // Initialize the transformation matrix to the identity matrix
    xform_matrix[0][0] = cos(rad);
    xform_matrix[0][1] = sin(rad) * -1;
    xform_matrix[0][2] = 0;

    xform_matrix[1][0] = sin(rad);
    xform_matrix[1][1] = cos(rad);
    xform_matrix[1][2] = 0;

    xform_matrix[2][0] = 0;
    xform_matrix[2][1] = 0;
    xform_matrix[2][2] = 1;
}

// Convert the input angle from degrees to radians

// This function transforms an atom by applying the specified transformation matrix
// Input is a pointer to an atom and a 3x3 transformation matrix
// Output is the transformed atom
void mol_xform(molecule *molecule, xform_matrix matrix)
{
    // define variables
    int i;
    double x, y, z;
    // get the number of atoms in the molecule
    int num_atoms = molecule->atom_no;

    // loop through the atoms in the molecule
    for (i = 0; i < num_atoms; i++)
    {
        // get the x, y, and z coordinates of the atom
        x = molecule->atoms[i].x;
        y = molecule->atoms[i].y;
        z = molecule->atoms[i].z;

        // apply the transformation matrix to the atom coordinates and update the atom coordinates
        molecule->atoms[i].x = ((x * matrix[0][0]) + (y * matrix[0][1]) + (z * matrix[0][2]));
        molecule->atoms[i].y = ((x * matrix[1][0]) + (y * matrix[1][1]) + (z * matrix[1][2]));
        molecule->atoms[i].z = ((x * matrix[2][0]) + (y * matrix[2][1]) + (z * matrix[2][2]));
    }
    // apply rotation matrix to bond coordinate
    for (i = 0; i < molecule->bond_no; i++)
    {
        bondset(molecule->bonds + i, &molecule->bonds[i].a1, &molecule->bonds[i].a2, &(molecule->atoms), &molecule->bonds[i].epairs);
    }
}