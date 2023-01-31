#include <stdlib.h>
#include "mol.h"
#include <math.h>
#include <string.h>

// setter for an atom
void atomset(atom *atom, char element[3], double *x, double *y, double *z)
{
    // setting x,y,z
    atom->x = *x;
    atom->y = *y;
    atom->z = *z;
    // copies the items into the addresses.
    strncpy(atom->element, element, 3);
};

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
};

// bond setter
void bondset(bond *bond, atom *a1, atom *a2, unsigned char epairs)
{
    // setting the bonds a1 and a2
    bond->a1 = a1;
    bond->a2 = a2;
    bond->epairs = epairs;
};

// bond getter
void bondget(bond *bond, atom **a1, atom **a2, unsigned char *epairs)
{
    // getting the bonds a1 and a2
    *a1 = bond->a1;
    *a2 = bond->a2;
    *epairs = bond->epairs;
};

/* function to allocate memory for a molecule structure */
molecule *molmalloc(unsigned short atom_max, unsigned short bond_max)
{
    /* allocate memory for the molecule structure */
    struct molecule *mole = (molecule *)malloc(sizeof(molecule));

    /* check if memory allocation was successful */
    if (mole == NULL)
    {
        printf("<Memory Error> for Molecule Struct (allocating)\n"); // print error message
        exit(1);                                                     // exit the program
    }

    /* set the maximum number of atoms */
    mole->atom_max = atom_max;
    /* initialize the number of atoms to zero */
    mole->atom_no = 0;
    /* allocate memory for the atoms array */
    mole->atoms = (atom *)malloc(atom_max * sizeof(atom));

    /* check if memory allocation was successful */
    if (mole->atoms == NULL)
    {
        printf("<Memory Error> for Atom Array (allocating)\n"); // print error message
        exit(1);                                                // exit the program
    }
    /* allocate memory for the atom pointers array */
    mole->atom_ptrs = (atom **)malloc(atom_max * sizeof(atom *));
    /* check if memory allocation was successful */
    if (mole->atom_ptrs == NULL)
    {
        printf("<Memory Error> for Atom Pointer (allocating)\n"); // print error message
        exit(1);                                                  // exit the program
    }

    /* set the maximum number of bonds */
    mole->bond_max = bond_max;
    /* initialize the number of bonds to zero */
    mole->bond_no = 0;
    /* allocate memory for the bonds array */
    mole->bonds = (bond *)malloc(bond_max * sizeof(bond));
    /* check if memory allocation was successful */
    if (mole->bonds == NULL)
    {
        printf("<Memory Error> for Bond Array (allocating)\n"); // print error message
        exit(1);                                                // exit the program
    }
    /* allocate memory for the bond pointers array */
    mole->bond_ptrs = (bond **)malloc(bond_max * sizeof(bond *));
    /* check if memory allocation was successful */
    if (mole->bond_ptrs == NULL)
    {
        printf("<Memory Error> for Bond Pointer (allocating)\n"); // print error message
        exit(1);                                                  // exit the program
    }

    /* return the pointer to the allocated memory */
    return mole;
};

molecule *molcopy(molecule *src)
{
    // error checker if the source molecule is empty or NULL
    if (src == NULL)
    {
        printf("<Error> Source Molecule Struct is Empty\n");
        return NULL;
    }
    // error checker if atoms and atom_ptr are empty or NULL
    if (src->atoms == NULL || src->atom_ptrs == NULL || src->bonds == NULL || src->bond_ptrs == NULL)
    {
        printf("<Error> Source Molecule Struct is not Initialized Properly\n");
        return NULL;
    }
    // assigning the necessary space for the new molecule
    molecule *new_mole = molmalloc(src->atom_max, src->bond_max);
    // Null checker for new molecule
    if (new_mole == NULL)
    {
        printf("<Error> Failure to allocate memory for new molecule\n");
        return NULL;
    }

    // copying the pointers
    for (int i = 0; i < src->atom_no; i++)
    {
        new_mole->atom_ptrs[i] = &(new_mole->atoms[i]);
    }
    // copying the pointers
    for (int i = 0; i < src->bond_no; i++)
    {
        new_mole->bond_ptrs[i] = &(new_mole->bonds[i]);
    }
    return new_mole;
};

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

// adds it to the end
void molappend_atom(molecule *mol, atom *at)
{
    int i;
    // Loop through the atom pointers in the molecule
    for (i = 0; i < mol->atom_no; i++)
    {
        // If there is an empty slot, add the atom to it
        if (mol->atom_ptrs[i] == NULL)
        {
            mol->atoms[i] = *at;
            break;
        }
    }
    // Add the address of the newly added atom to the atom pointers in the molecule
    mol->atom_ptrs[mol->atom_no] = &mol->atoms[i];
    // Increment the atom number in the molecule
    mol->atom_no++;
    // If the number of atoms in the molecule equals the max number of atoms
    if (mol->atom_no == mol->atom_max)
    {
        // If the max number of atoms is 0, set it to 1
        if (mol->atom_max == 0)
        {
            mol->atom_max = 1;
        }
        // Otherwise, double the max number of atoms
        else
        {
            mol->atom_max *= 2;
        }
        // Reallocate memory for the atoms array in the molecule to be twice its current size
        mol->atoms = realloc(mol->atoms, sizeof(at) * mol->atom_max);
        // Reallocate memory for the atom pointers array in the molecule to be twice its current size

        mol->atom_ptrs = realloc(mol->atom_ptrs, sizeof(*at) * mol->atom_max);
        // Update the addresses of the atoms in the molecule

        for (i = 0; i < mol->atom_no; i++)
        {
            mol->atom_ptrs[i] = &mol->atoms[i];
        }
    }
};

// Function to append a bond to the end of a molecule
void molappend_bond(molecule *mol, bond *bonder)
{
    int i;
    // Loop to find the first NULL element in bond_ptrs array
    for (i = 0; i < mol->bond_no; i++)
    {
        if (mol->bond_ptrs[i] == NULL)
        {
            // Assign the bond to the first NULL element in bonds array
            mol->bonds[i] = *bonder;
            break;
        }
    }
    // Add the address of the bond to the bond_ptrs array
    mol->bond_ptrs[mol->bond_no] = &mol->bonds[i];
    // Increment the number of bonds in the molecule
    mol->bond_no++;

    // If the number of bonds equals the maximum number of bonds, reallocate memory
    if (mol->bond_no == mol->bond_max)
    {
        // If the maximum number of bonds is 0, set it to 1
        if (mol->bond_max == 0)
        {
            mol->bond_max = 1;
        }
        else
        {
            // Double the maximum number of bonds

            mol->bond_max *= 2;
        }
        // Reallocate memory for the bonds and bond_ptrs arrays

        mol->bonds = realloc(mol->bonds, sizeof(bonder) * mol->bond_max);
        mol->bond_ptrs = realloc(mol->bond_ptrs, sizeof(*bonder) * mol->bond_max);
        // Update the addresses of the bonds in the bond_ptrs array
        for (i = 0; i < mol->bond_no; i++)
        {
            mol->bond_ptrs[i] = &mol->bonds[i];
        }
    }
}

// cmp function compares two atom pointers based on their z value
int atom_cmp(const void *a, const void *b)
{
    // Cast the void pointers to atom pointers
    atom **atom1 = (atom **)a;
    // Cast the void pointers to atom pointers
    atom **atom2 = (atom **)b;
    // Compare the z values of the atoms
    return ((*atom1)->z - (*atom2)->z);
}

int bond_cmp(const void *a, const void *b)
{
    bond **bond_a = (bond **)a;
    bond **bond_b = (bond **)b;
    double z_a = ((*bond_a)->a1->z + (*bond_a)->a2->z) / 2;
    double z_b = ((*bond_b)->a1->z + (*bond_b)->a2->z) / 2;
    return (z_a - z_b);
}

void molsort(molecule *molecule)
{
    qsort(molecule->atom_ptrs, molecule->atom_no, sizeof(atom *), atom_cmp);
    qsort(molecule->bond_ptrs, molecule->bond_no, sizeof(bond *), bond_cmp);
}

// Function to apply a rotation transformation in the X-axis to a given transformation matrix
// deg - angle of rotation in degrees
// xform_matrix - transformation matrix to be updated
void xrotation(xform_matrix xform, unsigned short deg)
{
    // Convert the rotation angle from degrees to radians
    double rad = deg * M_PI / 180;

    // Initialize the transformation matrix to the identity matrix
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            xform[i][j] = (i == j) ? 1 : 0;
        }
    }

    // Update the transformation matrix to reflect a rotation around the x-axis
    xform[1][1] = cos(rad);
    xform[1][2] = -sin(rad);
    xform[2][1] = sin(rad);
    xform[2][2] = cos(rad);
}

// This function calculates the transformation matrix for y-axis rotation
// The input is a 3x3 transformation matrix and an angle in degrees
// The output is the transformation matrix with y-axis rotation applied
void yrotation(xform_matrix xform, unsigned short deg)
{
    double rad = deg * M_PI / 180; // Convert degrees to radians

    // initailize the transformation matrix to the identity matrix
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            xform[i][j] = 0;
        }
    }
    xform[0][0] = cos(rad);
    xform[0][2] = sin(rad);
    xform[2][0] = -sin(rad);
    xform[2][2] = cos(rad);
    xform[1][1] = 1;
}

void zrotation(xform_matrix xform, unsigned short deg)
{
    // Convert the rotation angle from degrees to radians
    double rad = deg * M_PI / 180;

    // Initialize the transformation matrix to the identity matrix
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            xform[i][j] = (i == j) ? 1 : 0;
        }
    }

    // Update the transformation matrix to reflect a rotation around the z-axis
    xform[0][0] = cos(rad);
    xform[0][1] = -sin(rad);
    xform[1][0] = sin(rad);
    xform[1][1] = cos(rad);
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
        atom *curr_atom = &molecule->atoms[i];
        double x = curr_atom->x;
        double y = curr_atom->y;
        double z = curr_atom->z;
        curr_atom->x = x * matrix[0][0] + y * matrix[0][1] + z * matrix[0][2];
        curr_atom->y = x * matrix[1][0] + y * matrix[1][1] + z * matrix[1][2];
        curr_atom->z = x * matrix[2][0] + y * matrix[2][1] + z * matrix[2][2];
    }
}
