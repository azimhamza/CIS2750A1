#include "mol.h"

#ifndef STDIO_H
#include <stdio.h>
#endif // !STDIO_H

void display_mol(molecule *mol)
{
    int i;
    // display atom array params
    printf("atom_max: %hu, atom_no: %hu\n", mol->atom_max, mol->atom_no);
    // display all atoms
    for (i = 0; i < mol->atom_no; i++)
    {
        char name[3];
        double x, y, z;
        atom *atom_ptr;

        atom_ptr = mol->atoms + i;

        atomget(atom_ptr, name, &x, &y, &z);
        printf("atom[%d] %p %-2s: %lf %lf %lf\n",
               i, (void *)atom_ptr, name, x, y, z);
    }

    // display atom pointers and the atoms that they point to
    for (int i = 0; i < mol->atom_no; i++)
    {
        char name[3];
        double x, y, z;
        atom *atom_ptr;

        // get pointer from atom_ptrs array
        atom_ptr = mol->atom_ptrs[i];

        printf("atom_ptrs[%d]: %p", i, (void *)(atom_ptr));

        // get the info form the atom
        atomget(atom_ptr, name, &x, &y, &z);
        printf("  -> atom[%lu] %p %-2s: %lf %lf %lf\n",
               (atom_ptr - mol->atoms),
               (void *)atom_ptr, name, x, y, z);
    }

    // bond params
    printf("bond_max: %hu, bond_no: %hu\n", mol->bond_max, mol->bond_no);

    // display all bonds
    for (i = 0; i < mol->bond_no; i++)
    {
        bond *bond_ptr;
        unsigned short a1, a2;
        unsigned char epairs;
        atom *atoms;
        bond_ptr = mol->bonds + i;
        bondget(bond_ptr, &a1, &a2, &atoms, &epairs);
        // printf("bond[%d] %p: %p(index %hu)(%s) %p(index %hu)(%s) %hhu\n", i, (void *)bond_ptr, (void *)(bond_ptr->atoms + a1), a1, bond_ptr->atoms[a1].element, (void *)(bond_ptr->atoms + a2), a2, bond_ptr->atoms[a2].element, epairs);
        printf("bond[%d] (%p): z: %lf, len: %lf, dx: %lf, dy: %lf, atoms: %p\n",
               i, (void *)bond_ptr, bond_ptr->z, bond_ptr->len, bond_ptr->dx, bond_ptr->dy, (void *)atoms);
        printf("\ta1: atom[%hu](%s), x1: %lf, y1: %lf\n",
               a1, bond_ptr->atoms[a1].element, bond_ptr->x1, bond_ptr->y1);
        printf("\ta2: atom[%hu](%s), x2: %lf, y2: %lf\n",
               a2, bond_ptr->atoms[a2].element, bond_ptr->x2, bond_ptr->y2);
    }

    // display bond pointers & what they point to
    for (int i = 0; i < mol->bond_no; i++)
    {
        bond *bond_ptr;

        bond_ptr = mol->bond_ptrs[i];
        printf("bond_ptrs[%d]: %p", i, (void *)(mol->bond_ptrs[i]));

        // printf("  -> bond[%d] %p: %p(%s) %p(%s) %hhu\n",
        printf("  -> %p(%s) %p(%s) %hhu\n",
               // i, (void *)bond_ptr,
               //(void *)bond_ptr->a1,
               (void *)(bond_ptr->atoms + bond_ptr->a1),
               // bond_ptr->a1->element,
               bond_ptr->atoms[bond_ptr->a1].element,
               //(void *)bond_ptr->a2,
               (void *)(bond_ptr->atoms + bond_ptr->a2),
               // bond_ptr->a2->element,
               bond_ptr->atoms[bond_ptr->a2].element,
               bond_ptr->epairs);
    }
}

int main(int argc, char **argv)
{
    molecule *mol;
    atom a1, a2;
    bond single;
    unsigned char epair = 1;
    unsigned short locn0 = 0, locn1 = 1;
    double x = 1.0, y = 0.5, z = 1.5;

    /*
     Creating the molecule.  This tests:
     - atomset
     - bondset
     - molappend_atom
     - molappend bond
     It also indirectly tests:
     - compute_coords
    */
    atomset(&a1, "O", &x, &y, &z);
    x = 0.0;
    y = 1.0;
    z = 0.5;
    atomset(&a2, "H", &x, &y, &z);

    mol = molmalloc(1, 0);

    molappend_atom(mol, &a1);
    molappend_atom(mol, &a2);

    bondset(&single, &locn0, &locn1, &(mol->atoms), &epair);

    molappend_bond(mol, &single);

    display_mol(mol);
    /* End of Creating the molecule */

    /*
     Sorting the molecule. This tests:
     - molsort.
     It also indirectly tests:
     - bond_comp
    */
    molsort(mol);
    printf("\nsorted molecule:\n");
    display_mol(mol);
    /* End of Sorting the molecule */

    /*
     Copying the molecule. This tests:
     - molcopy
    */
    molecule *xrot;
    xrot = molcopy(mol);
    printf("\ncopied molecule (same molecule as x-rotated molecule below):\n");
    display_mol(xrot);
    /* End of Copying the molecule */

    /*
     Rotating the (copied) molecule on the x-axis. This tests:
     - xrotation
     - mol_xform
    */
    xform_matrix matrix;
    xrotation(matrix, 70);
    mol_xform(xrot, matrix); // invalid read
    printf("\nx-rotated molecule:\n");
    display_mol(xrot);
    /* End of Rotating the (copied) molecule on the x-axis */

    /*
     Rotating a copied molecule on the y-axis. This tests:
     - molcopy
     - yrotation
     - mol_xform
    */
    molecule *yrot;
    yrot = molcopy(mol);
    yrotation(matrix, 70);
    mol_xform(yrot, matrix);
    printf("\ny-rotated molecule:\n");
    display_mol(yrot);
    /* End of Rotating a copied molecule on the y-axis*/

    /*
     Rotating a copied molecule on the z-axis. This tests:
     - molcopy
     - zrotation
     - mol_xform
    */
    molecule *zrot;
    zrot = molcopy(mol);
    zrotation(matrix, 70);
    mol_xform(zrot, matrix);
    printf("\ny-rotated molecule:\n");
    display_mol(zrot);
    /* End of Rotating a copied molecule on the z-axis*/

    /* Free the molecules. This tests: molfree */
    molfree(mol);
    molfree(xrot);
    molfree(yrot);
    molfree(zrot);
    /* End of Free the molecules */

    return 0;
}
