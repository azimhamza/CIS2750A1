#include <stdlib.h>
#include "mol.h"
#include <math.h>
#include <string.h>

#define M_PI 3.14159265358979323846
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
    strncpy(element, atom->element, 3);
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

// creating space
molecule *molmalloc(unsigned short atom_max, unsigned short bond_max)
{
    molecule *mole = (molecule *)malloc(sizeof(molecule));
    // test case if its null
    if (mole == NULL)
    {
        printf("<Memory Error> for Molecule Struct (allocating)\n"); //
        exit(1);
    }
    // set atom_max
    mole->atom_max = atom_max;
    // initate array
    mole->atom_no = 0;
    // mallocing the atom, multiplying the total atom_max with the size of atom
    mole->atoms = (atom *)malloc(atom_max * sizeof(atom));

    // checks if the atoms are null
    if (mole->atoms == NULL)
    {
        printf("<Memory Error> for Atom Array (allocating)\n"); //
        exit(1);
    }
    // allocating pointer memory
    mole->atom_ptrs = (atom **)malloc(atom_max * sizeof(atom *));
    // null checker for atom_ptr
    if (mole->atom_ptrs == NULL)
    {
        printf("<Memory Error> for Atom Pointer (allocating)\n");
        exit(1);
    }

    // set atom max
    mole->bond_max = bond_max;
    // initate array
    mole->bond_no = 0;
    // mallocing the bond, multiplying the total bond_max with the size of bond
    mole->bonds = (bond *)malloc(bond_max * sizeof(bond));
    // checks if the bonds are null
    if (mole->bonds == NULL)
    {
        printf("<Memory Error> for Bond Array (allocating)\n");
        exit(1);
    }
    // allocating pointer memory
    mole->bond_ptrs = (bond **)malloc(bond_max * sizeof(bond *));
    // null checker for bond_ptr
    if (mole->bond_ptrs == NULL)
    {
        printf("<Memory Error> for Bond Pointer (allocating)\n");
        exit(1);
    }
    return mole;
};

molecule *molcopy(molecule *src)
{
    // error checker if the source moleculee is empty or NULL
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
    // assigning the new_mole atom_no and bond_no
    new_mole->atom_no = src->atom_no;
    new_mole->bond_no = src->bond_no;

    // copying all data for atoms, src->atom_no, src->atom_max, src->bond_no, src->bond_max

    memcpy(new_mole->atoms, src->atoms, src->atom_no * sizeof(atom));
    memcpy(new_mole->atom_ptrs, src->atom_ptrs, src->atom_no * sizeof(atom *));
    memcpy(new_mole->bonds, src->bonds, src->bond_no * sizeof(bond));
    memcpy(new_mole->bond_ptrs, src->bond_ptrs, src->bond_no * sizeof(bond *));

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
void molappend_atom(molecule *molecule, atom *atom)
{
    int i;
    // Loop through the atom pointers in the molecule
    for (i = 0; i < molecule->atom_no; i++)
    {
        // If there is an empty slot, add the atom to it
        if (molecule->atom_ptrs[i] == NULL)
        {
            molecule->atoms[i] = *atom;
            break;
        }
    }
    // Add the address of the newly added atom to the atom pointers in the molecule
    molecule->atom_ptrs[molecule->atom_no] = &molecule->atoms[i];
    // Increment the atom number in the molecule
    molecule->atom_no++;
    // If the number of atoms in the molecule equals the max number of atoms
    if (molecule->atom_no == molecule->atom_max)
    {
        // If the max number of atoms is 0, set it to 1
        if (molecule->atom_max == 0)
        {
            molecule->atom_max = 1;
        }
        // Otherwise, double the max number of atoms
        else
        {
            molecule->atom_max *= 2;
        }
        // Reallocate memory for the atoms array in the molecule to be twice its current size
        molecule->atoms = realloc(molecule->atoms, sizeof(atom) * molecule->atom_max);
        // Reallocate memory for the atom pointers array in the molecule to be twice its current size

        molecule->atom_ptrs = realloc(molecule->atom_ptrs, sizeof(*atom) * molecule->atom_max);
        // Update the addresses of the atoms in the molecule

        for (i = 0; i < molecule->atom_no; i++)
        {
            molecule->atom_ptrs[i] = &molecule->atoms[i];
        }
    }
};

// Function to append a bond to the end of a molecule
void molappend_bond(molecule *molecule, bond *bond)
{
    int i;
    // Loop to find the first NULL element in bond_ptrs array
    for (i = 0; i < molecule->bond_no; i++)
    {
        if (molecule->bond_ptrs[i] == NULL)
        {
            // Assign the bond to the first NULL element in bonds array
            molecule->bonds[i] = *bond;
            break;
        }
    }
    // Add the address of the bond to the bond_ptrs array
    molecule->bond_ptrs[molecule->bond_no] = &molecule->bonds[i];
    // Increment the number of bonds in the molecule
    molecule->bond_no++;

    // If the number of bonds equals the maximum number of bonds, reallocate memory
    if (molecule->bond_no == molecule->bond_max)
    {
        // If the maximum number of bonds is 0, set it to 1
        if (molecule->bond_max == 0)
        {
            molecule->bond_max = 1;
        }
        else
        {
            // Double the maximum number of bonds

            molecule->bond_max *= 2;
        }
        // Reallocate memory for the bonds and bond_ptrs arrays

        molecule->bonds = realloc(molecule->bonds, sizeof(bond) * molecule->bond_max);
        molecule->bond_ptrs = realloc(molecule->bond_ptrs, sizeof(*bond) * molecule->bond_max);
        // Update the addresses of the bonds in the bond_ptrs array
        for (i = 0; i < molecule->bond_no; i++)
        {
            molecule->bond_ptrs[i] = &molecule->bonds[i];
        }
    }
}

// cmp function compares two atom pointers based on their z value
int cmp(const void *a, const void *b)
{
    // calculate z value of first atom
    double a_z = (*(atom **)a)->z;
    // calculate z value of second atom
    double b_z = (*(atom **)b)->z;
    // return 1 if the z value of first atom is greater than second
    if (a_z > b_z)
        return 1;
    // return -1 if the z value of first atom is less than second
    else if (a_z < b_z)
        return -1;
    // return 0 if the z values of both atoms are equal
    else
        return 0;
}

// bondcmp function compares two bond pointers based on their z value
int bondcmp(const void *a, const void *b)
{
    // calculate z value of first bond by taking average of the z values of its two atoms
    double a_z = ((*(bond **)a)->a1->z + (*(bond **)a)->a2->z) / 2;
    // calculate z value of second bond by taking average of the z values of its two atoms
    double b_z = ((*(bond **)b)->a1->z + (*(bond **)b)->a2->z) / 2;
    // return 1 if the z value of first bond is greater than second
    if (a_z > b_z)
        return 1;
    // return -1 if the z value of first bond is less than second
    else if (a_z < b_z)
        return -1;
    // return 0 if the z values of both bonds are equal
    else
        return 0;
}

// sort atoms and bonds in a molecule based on the Z coordinate value
void molsort(molecule *molecule)
{
// sort atoms using qsort and the comparison function cmp
    qsort(molecule->atom_ptrs, molecule->atom_no, sizeof(atom *), cmp);
// sort bonds using qsort and the comparison function bondcmp
    qsort(molecule->bond_ptrs, molecule->bond_no, sizeof(bond *), bondcmp);
}

// Function to apply a rotation transformation in the X-axis to a given transformation matrix
// deg - angle of rotation in degrees
// xform_matrix - transformation matrix to be updated
void xrotation(xform_matrix xform, unsigned short deg) {
float rad = deg * M_PI / 180.0f;
xform[0][0] = 1;
xform[1][1] = cos(rad);
xform[1][2] = -sin(rad);
xform[2][1] = sin(rad);
xform[2][2] = cos(rad);
xform[3][3] = 1;
}

// This function calculates the transformation matrix for y-axis rotation
// The input is a 3x3 transformation matrix and an angle in degrees
// The output is the transformation matrix with y-axis rotation applied
void yrotation(xform_matrix xform, unsigned short deg) {
    float rad = deg * M_PI / 180.0f;
    xform[0][0] = cos(rad);
    xform[0][2] = sin(rad);
    xform[1][1] = 1;
    xform[2][0] = -sin(rad);
    xform[2][2] = cos(rad);
    xform[3][3] = 1;
}

void zrotation(xform_matrix xform, unsigned short deg) {
    float rad = deg * M_PI / 180.0f;
    xform[0][0] = cos(rad);
    xform[0][1] = -sin(rad);
    xform[1][0] = sin(rad);
    xform[1][1] = cos(rad);
    xform[2][2] = 1;
    xform[3][3] = 1;
}


// This helper function computes the transformation matrix for rotation about an arbitrary axis
// Input is the 3x3 transformation matrix and an angle in radians
// Output is the transformation matrix with rotation applied

void computeRotationMatrix(xform_matrix xform_matrix, double rad)
{
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

// This function computes the transformation matrix for rotation about the z-axis
// Input is the 3x3 transformation matrix and an angle in degrees
// Output is the transformation matrix with z-axis rotation applied

// Convert the input angle from degrees to radians



// This function transforms an atom by applying the specified transformation matrix
// Input is a pointer to an atom and a 3x3 transformation matrix
// Output is the transformed atom
void transform_atom(atom *a, xform_matrix matrix)
{
    double x = a->x;
    double y = a->y;
    double z = a->z;

    a->x = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z;
    a->y = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z;
    a->z = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z;
}

// This function transforms a molecule by applying the specified transformation matrix to each atom
// Input is a pointer to a molecule and a 3x3 transformation matrix
// Output is the transformed molecule
void mol_xform(molecule *molecule, xform_matrix matrix)
{
    for (int i = 0; i < molecule->atom_no; i++)
    {
        transform_atom(molecule->atom_ptrs[i], matrix);
    }
}
