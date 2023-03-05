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


// updated bond setter
void bondset(bond *bond, unsigned short *a1, unsigned short *a2, atom **atoms, unsigned char *epairs)
{
    // setting the bonds a1 and a2  inside atoms and epairs
    bond->a1 = atoms[*a1];
    bond->a2 = atoms[*a2];
    bond->atoms = *atoms;
    bond->epairs = *epairs;
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
    
    // allocating memory for the atoms and atom pointers and checking for errors
    ptr->atom_max = atom_max;
    ptr->atom_no = 0;
    ptr->atoms = (atom *)malloc(atom_max * sizeof(atom));
    if (ptr->atoms == NULL)
    {
        perror("null atoms in molmalloc <Error>");
        free(ptr); // free the molecule memory before returning NULL
        return NULL;
    }
    ptr->atom_ptrs = (atom **)malloc(atom_max * sizeof(atom *));
    if (ptr->atom_ptrs == NULL)
    {
        perror("null atom_ptrs in molmalloc <Error>");
        free(ptr->atoms); // free the atoms memory before returning NULL
        free(ptr); // free the molecule memory before returning NULL
        return NULL;
    }

    // allocating memory for the bonds and bond pointers and checking for errors
    ptr->bond_max = bond_max;
    ptr->bond_no = 0;
    ptr->bonds = (bond *)malloc(bond_max * sizeof(bond));
    if (ptr->bonds == NULL)
    {
        perror("null bonds in molmalloc <Error>");
        free(ptr->atom_ptrs); // free the atom pointers memory before returning NULL
        free(ptr->atoms); // free the atoms memory before returning NULL
        free(ptr); // free the molecule memory before returning NULL
        return NULL;
    }
    ptr->bond_ptrs = (bond **)malloc(bond_max * sizeof(bond *));
    if (ptr->bond_ptrs == NULL)
    {
        perror("null bond_ptrs in molmalloc <Error>");
        free(ptr->bonds); // free the bonds memory before returning NULL
        free(ptr->atom_ptrs); // free the atom pointers memory before returning NULL
        free(ptr->atoms); // free the atoms memory before returning NULL
        free(ptr); // free the molecule memory before returning NULL
        return NULL;
    }

    // return the molecule
    return ptr;
}


// copy the molecule to a new molecule and return the new molecule
molecule *molcopy(molecule *src)
{
    molecule *new_mole = (molecule *)malloc(sizeof(molecule));
    new_mole->atom_max = src->atom_max;
    new_mole->atom_no = src->atom_no;
    new_mole->bond_max = src->bond_max;
    new_mole->bond_no = src->bond_no;
    new_mole->atoms = (atom *)malloc(src->atom_max * sizeof(atom));
    new_mole->atom_ptrs = (atom **)malloc(src->atom_max * sizeof(atom *));
    new_mole->bonds = (bond *)malloc(src->bond_max * sizeof(bond));
    new_mole->bond_ptrs = (bond **)malloc(src->bond_max * sizeof(bond *));
    for (int i = 0; i < src->atom_no; i++)
    {
        molappend_atom(new_mole, &src->atoms[i]);
    }
    for (int i = 0; i < src->bond_no; i++)
    {
        molappend_bond(new_mole, &src->bonds[i]);
    }
    return new_mole;
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

// append a bond to the molecule structure pointed to by molecule and return the bond pointer
void molappend_atom(molecule *molecule, atom *at)
{
    // variables for the new capacity and the index
    int new_capacity;
    int i;

    // Check if the current atom_no equals atom_max
    if (molecule->atom_no == molecule->atom_max)
    {
        // If atom_max is 0, set it to 1
        if (molecule->atom_max == 0)
        {
            // Set atom_max to 1
            molecule->atom_max = 1;
        }
        // Else, double atom_max
        else
        {
            molecule->atom_max *= 2;
        }
        // Increase the capacity of atoms and atom_ptrs arrays
        new_capacity = molecule->atom_max * sizeof(atom);
        // Reallocate memory for the atoms array
        molecule->atoms = realloc(molecule->atoms, new_capacity);
        // Reallocate memory for the atom_ptrs array
        new_capacity = molecule->atom_max * sizeof(atom *);
        // Reallocate memory for the atom_ptrs array
        molecule->atom_ptrs = realloc(molecule->atom_ptrs, new_capacity);

        // Update the pointers in atom_ptrs to point to the corresponding atoms in the new atoms array
        for (i = 0; i < molecule->atom_no; i++)
        {
            // Set the pointer in atom_ptrs to point to the corresponding atom in the atoms array
            molecule->atom_ptrs[i] = &(molecule->atoms[i]);
        }
    }
    // Copy the data pointed to by atom to the first "empty" atom in atoms in the molecule pointed to by molecule
    molecule->atoms[molecule->atom_no] = *at;
    molecule->atom_ptrs[molecule->atom_no] = &(molecule->atoms[molecule->atom_no]);
    // Increment the value of atom_no
    molecule->atom_no++;
}

void molappend_bond(molecule *molecule, bond *bo)
{
    // Check if the current bond_no equals bond_max
    if (molecule->bond_no == molecule->bond_max)
    {
        // If bond_max is 0, set it to 1
        if (molecule->bond_max == 0)
        {
            molecule->bond_max = 1;
            molecule->bonds = (struct bond *)malloc(sizeof(struct bond) * molecule->bond_max);
            molecule->bond_ptrs = (struct bond **)malloc(sizeof(struct bond *) * molecule->bond_max);
        }
        if (molecule->bond_max == molecule->bond_no)
        {
            molecule->bond_max * 2;
            molecule->bonds = (struct bond *)realloc(molecule->bonds, sizeof(struct bond) * molecule->bond_max);
            molecule->bond_ptrs = (struct bond **)realloc(molecule->bond_ptrs, sizeof(struct bond *) * molecule->bond_max);

            for (int i = 0; i < molecule->bond_no; i++)
            {
                molecule->bond_ptrs[i] = molecule->bonds + (i);
            }
        }
        // Copy the data pointed to by bond to the first "empty" bond in bonds in the molecule pointed to by molecule
        molecule->bonds[molecule->bond_no] = *bo;
        molecule->bond_ptrs[molecule->bond_no] = &(molecule->bonds[molecule->bond_no]);
        // Increment the value of bond_no
        molecule->bond_no++;
    }
}

// void print_molecule function prints the molecule structure pointed to by molecule to the file pointed to by fp.
int bond_comp( const void *a, const void *b){
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

    // call all the values from bond
    atom *a1 = bond->a1;
    atom *a2 = bond->a2;
    // storing x1 and x2
    bond->x1 = bond->atoms[bond->a1].x;
    bond->x2 = bond->atoms[bond->a2].x;
    // storing y1 and y2
    bond->y1 = bond->atoms[bond->a1].y;
    bond->y2 = bond->atoms[bond->a2].y;
    // computing z as the average distance between the two atoms
    bond->z = (bond->atoms[bond->a1].z + bond->atoms[bond->a2].z) / 2;
    // computing the length of the bond
    bond->z = sqrt(pow(bond->x2 - bond->x1, 2) + pow(bond->y2 - bond->y1, 2));
    // computing dx and dy
    bond->dx = (bond->x2 - bond->x1) / bond->len;
    bond->dy = (bond->y2 - bond->y1) / bond->len;
}

// sort the atoms and bonds in the molecule pointed to by molecule by z value
void molsort(molecule *molecule)
{
    // set the compare function to the compare_atom_z function
    int (*compare)(const void *a, const void *b);
    // set the compare function to the compare_atom_z function
    compare = compare_atom_z;
    // qsort the atoms and bonds in the molecule pointed to by molecule by z value
    qsort(molecule->atom_ptrs, molecule->atom_no, sizeof(atom *), compare);
    // set the compare function to the compare_bond_z function
    compare = compare_bond_z;
    // qsort the atoms and bonds in the molecule pointed to by molecule by z value
    qsort(molecule->bond_ptrs, molecule->bond_no, sizeof(bond *), compare);
}

/* helper function for molsort to compare z values of atoms */
int compare_atom_z(const void *a, const void *b)
{
    // cast the void pointers to atom pointers
    atom *aa = *(atom **)a;
    atom *bb = *(atom **)b;
    // get the z values of the atoms
    double za = aa->z;
    double zb = bb->z;
    // compare the z values if za < zb return -1, if za > zb return 1, else return 0
    if (za < zb)
        return -1;
    else if (za > zb)
        return 1;
    return 0;
}

/* helper function for molsort to compare z values of bonds */
int compare_bond_z(const void *a, const void *b)
{
    // cast the void pointers to bond pointers
    bond *ab = *(bond **)a;
    bond *bb = *(bond **)b;
    // get the z values of the bonds
    double za = ab->z;
    double zb = bb->z;
    // compare the z values if za < zb return -1, if za > zb return 1, else return 0
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
    xform_matrix[0][1] = -sin(rad);
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
   for ( i = 0; i <molecule->bond_no; i++){
        bondset(molecule->bonds + i, &molecule->bonds[i].a1, &molecule->bonds[i].a2, &(molecule->atoms), &molecule->bonds[i].epairs);
   }
}