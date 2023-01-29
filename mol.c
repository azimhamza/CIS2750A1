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
    for (i = 0; i < molecule->atom_no; i++)
    {
        if (molecule->atom_ptrs[i] == NULL)
        {
            molecule->atoms[i] = *atom;
            break;
        }
    }
    molecule->atom_ptrs[molecule->atom_no] = &molecule->atoms[i];
    molecule->atom_no++;

    if (molecule->atom_no == molecule->atom_max)
    {
        if (molecule->atom_max == 0)
        {
            molecule->atom_max = 1;
        }
        else
        {
            molecule->atom_max *= 2;
        }
        molecule->atoms = realloc(molecule->atoms, sizeof(atom) * molecule->atom_max);
        molecule->atom_ptrs = realloc(molecule->atom_ptrs, sizeof(*atom) * molecule->atom_max);
        for (i = 0; i < molecule->atom_no; i++)
        {
            molecule->atom_ptrs[i] = &molecule->atoms[i];
        }
    }
};

void molappend_bond(molecule *molecule, bond *bond)
{
    int i;
    for (i = 0; i < molecule->bond_no; i++)
    {
        if (molecule->bond_ptrs[i] == NULL)
        {
            molecule->bonds[i] = *bond;
            break;
        }
    }
    molecule->bond_ptrs[molecule->bond_no] = &molecule->bonds[i];
    molecule->bond_no++;

    if (molecule->bond_no == molecule->bond_max)
    {
        if (molecule->bond_max == 0)
        {
            molecule->bond_max = 1;
        }
        else
        {
            molecule->bond_max *= 2;
        }
        molecule->bonds = realloc(molecule->bonds, sizeof(bond) * molecule->bond_max);
        molecule->bond_ptrs = realloc(molecule->bond_ptrs, sizeof(*bond) * molecule->bond_max);
        for (i = 0; i < molecule->bond_no; i++)
        {
            molecule->bond_ptrs[i] = &molecule->bonds[i];
        }
    }
}

int cmp(const void *a, const void *b)
{
    double a_z = (*(atom **)a)->z;
    double b_z = (*(atom **)b)->z;
    if (a_z > b_z)
        return 1;
    else if (a_z < b_z)
        return -1;
    else
        return 0;
}

int bondcmp(const void *a, const void *b)
{
    double a_z = ((*(bond **)a)->a1->z + (*(bond **)a)->a2->z) / 2;
    double b_z = ((*(bond **)b)->a1->z + (*(bond **)b)->a2->z) / 2;
    if (a_z > b_z)
        return 1;
    else if (a_z < b_z)
        return -1;
    else
        return 0;
}

void molsort(molecule *molecule)
{

    qsort(molecule->atom_ptrs, molecule->atom_no, sizeof(atom *), cmp);

    qsort(molecule->bond_ptrs, molecule->bond_no, sizeof(bond *), bondcmp);
}

void xrotation(xform_matrix xform_matrix, unsigned short deg)
{
    // Convert deg to radians
    double rad = deg * (M_PI / 180);
    // Compute the transformation matrix
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

void yrotation(xform_matrix xform_matrix, unsigned short deg)
{
    // Convert deg to radians
    double rad = deg * (M_PI / 180);
    // Compute the transformation matrix
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

void zrotation(xform_matrix xform_matrix, unsigned short deg)
{
    // Convert deg to radians
    double rad = deg * (M_PI / 180);
    // Compute the transformation matrix
    computeRotationMatrix(xform_matrix, rad);
}

void transform_atom(atom *a, xform_matrix matrix)
{
    double x = a->x;
    double y = a->y;
    double z = a->z;

    a->x = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z;
    a->y = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z;
    a->z = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z;
}

void mol_xform(molecule *molecule, xform_matrix matrix)
{
    for (int i = 0; i < molecule->atom_no; i++)
    {
        transform_atom(molecule->atom_ptrs[i], matrix);
    }
}
