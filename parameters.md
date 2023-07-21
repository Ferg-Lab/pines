#### Supported parameters in PINES

##### Template
```
PINES ... 
LABEL=pines             
PRECISION=10000        
NLIST                 
ONLYCROSS               
REF_FILE=ref.pdb       
PINESATOMS=13          
ATOMTYPES=QC1,QC2,QP1,1N1,1N2,1N3,QP2,2N1,2N2,2N3,QP3,3N1,3N2,3N3  
SFACTOR=1.0
SORT=1.0
SWITCH1={RATIONAL R_0=0.56 MM=10 NN=5}
NL_CUTOFF=4.0,4.0
NL_STRIDE=1.0
NL_SKIN=0.1,0.1
NL_CONSTANT_SIZE = 1
WIRITEPIVESTRAJ 
... PINES
```

##### Notes
- LABEL:
Label for referencing this module elsewhere in the PLUMED file.
  
- PRECISION: For applying integer sorting algorithm. Precision is used to convert floats to integers. For example, a precision of 100000 implies a float of 0.9999 is 9999. Similarly, a float of 0.99999 to 9999. A minimum of PRECISION=100 is required.
  
- NLIST: Unused parameter.
  
- ONLYCROSS: Consider atoms pairs of different atom types in constructing PIV. Other supported value is ONLYDIRECT, which considers identical atom type pairs for constructing PIV. ONLYDIRECT is suitable for Argon system containing identical particles and is tested only with that system.
  
- REF_FILE: Input reference PDB file. It should contain the same order of atoms with same indices as the structure file used for performing the simulations. Atom names in PDB could be modified as needed to define the ATOMTYPES values set below for constructing the PIV.
  
- PINESATOMS: Total number of atom types. This is the same as the length of ATOMTYPES and is used for detecting errors in ATOMTYPES declaration.
  
- ATOMTYPES: All atomtypes of solute + solvent for defining the blocks in PIV. In the ONLYCROSS case, the order of the final PIV is based on the sequence of the atom types defined here. In brief, the order of PIV is solute1-solute2, solute1-solute3, ..., solute1-soluteN, solute1-OW [of size equal to NL_CONSTANT_SIZE defined below], solute1-H1 [of size equal to NL_CONSTANT_SIZE defined below], solute2-solute3, solute2-solute4, ..., solute2-soluteN, solute2-OW [of size equal to NL_CONSTANT_SIZE defined below], solute2-H1 [of size equal to NL_CONSTANT_SIZE defined below], ..., soluteN-OW [of size equal to NL_CONSTANT_SIZE defined below], soluteN-H1 [of size equal to NL_CONSTANT_SIZE defined below]. The oxygen atoms are extracted based on NL_CUTOFF and only the number of oxygen atoms equal to NL_CONSTANT_SIZE are retained. Similarly, hydrogen atoms are extracted based on atom indices that share the residue numbers with the oxygen atom indices and only NL_CONSTANT_SIZE hydrogen atoms are retained. Consequently, the total number of elements in PIV is equal to (PINESATOMS-3)x(PINESATOMS-3)x0.5 + (PINESATOMS-3)xNL_CONSTANT_SIZE + 2x(PINESATOMS-3)xNL_CONSTANT_SIZE. The number of blocks is (PINESATOMS-3)x(PINESATOMS-3)x0.5 + (PINESATOMS-3) + (PINESATOMS-3)
  
- SFACTOR: Should be the same size as of nblocks = (PINESATOMS-3)x(PINESATOMS-3)x0.5 + (PINESATOMS-3) + (PINESATOMS-3)
  
- SORT: Should be the same size as of nblocks = (PINESATOMS-3)x(PINESATOMS-3)x0.5 + (PINESATOMS-3) + (PINESATOMS-3) 
  
- SWITCH1: Should be the same size as of nblocks = (PINESATOMS-3)x(PINESATOMS-3)x0.5 + (PINESATOMS-3) + (PINESATOMS-3)
  
- NL_CUTOFF: Neighbor list cutoff for each element in PIV. Same size as nelements. Use a small value but update the neighbor list frequently, which can be expensive. Alternatively, use a high value and update less frequently to save computational time.
  
- NL_STRIDE: Frequency to update neighbor list. For post-processing, set it to 1 to read each frame.
  
- NL_SKIN: Additional buffer cutoff to nl_cutoff. Neighborlist reconstruction is triggered if atoms move out of a radius equal to nlcutoff+skin.
  
- NL_CONSTANT_SIZE: Number of oxygen and hydrogen atoms to retain in each solute-solvent block of PIV.
  
- WIRITEPIVESTRAJ: Write PIV in each frame. Alternatively, use `print STRIDE=1 ARG=pines.* FILE=colvar.out`. Depending on switching function parameters, some of the values could be ~0.9999. If the value exceeds the precision, then these values will be stored as NAN after converting to integer. In that case colvar reports these values as zeros but WIRITEPIVESTRAJ could be used for reporting the exact values.

##### Miscellaneous

- Use unique atom names in the input reference PDB file to let `PINES` recognize the atom types of solute. Only `OW`, `HW1` and `HW2` are currently recognized as atom types in PINES for solvent. The input reference PDB file is primarily used for fetching the indices of atoms in each element of PIV. Therefore, it can be a separate file than that is used for performing the simulations.
- Further information on the parameters can be found in [src/pines/PINES.cpp](src/pines/PINES.cpp) or see below snapshot.

![pines](https://github.com/Ferg-Lab/pines/assets/38693318/c197d449-7a22-4111-be78-911f77170a91)

