# CIS2750A2

The following has been changed in the code, 

In the molsql.py, the add_bond, in the table a bond_ID was added and subsequently fixed. 
In the MolDisplay, the indexing was incorrect and was subsequently fixed

Program works correctly within the MACOS system, but for some reason the on the LINUX server
there is a segmentation fault and upon further inspection the following error occurred: 

Program received signal SIGSEGV, Segmentation fault.
__memmove_avx_unaligned_erms () at ../sysdeps/x86_64/multiarch/memmove-vec-unaligned-erms.S:369
369     ../sysdeps/x86_64/multiarch/memmove-vec-unaligned-erms.S: No such file or directory.>>

This error only occurs with the molsql.py file and could not be further assessed. 

To run the code, simply type make, and the code will generate the swig files