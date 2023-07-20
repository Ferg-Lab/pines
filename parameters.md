#### Supported parameters in PINES

```
PINES ... 
LABEL=pines             # Label for referencing this module elsewhere in the PLUMED file.
PRECISION=10000         # For applying integer sorting algorithm. Precision is used to convert floats to integers. For example, a precision of 100000 implies a float of 0.9999 is 9999. Similarly, a float of 0.99999 to 9999. A minimum of PRECISION=100 is required.
NLIST                   # Unused parameter.
ONLYCROSS               # Consider atoms pairs of different atom types in constructing PIV. Other supported value is ONLYDIRECT, which considers identical atom type pairs for constructing PIV. ONLYDIRECT is suitable for Argon system containing identical particles and is tested only with that system.
REF_FILE=ref.pdb        # Input reference PDB file. It should contain the same order of atoms with same indices as the structure file used for performing the simulations. Atom names in PDB could be modified as needed to define the ATOMTYPES values set below for constructing the PIV.
PINESATOMS=13           # Total number of atom types. This is the same as the length of ATOMTYPES and is used for detecting errors in ATOMTYPES declaration.
ATOMTYPES=QC1,QC2,QP1,1N1,1N2,1N3,QP2,2N1,2N2,2N3,QP3,3N1,3N2,3N3  # All atomtypes of solute + solvent for defining the blocks in PIV. In the ONLYCROSS case, the order of the final PIV is based on the sequence of the atom types defined here. In brief, the order of PIV is solute1-solute2, solute1-solute3, ..., solute1-soluteN, solute1-OW [of size equal to NL_CONSTANT_SIZE defined below], solute1-H1 [of size equal to NL_CONSTANT_SIZE defined below], solute2-solute3, solute2-solute4, ..., solute2-soluteN, solute2-OW [of size equal to NL_CONSTANT_SIZE defined below], solute2-H1 [of size equal to NL_CONSTANT_SIZE defined below], ..., soluteN-OW [of size equal to NL_CONSTANT_SIZE defined below], soluteN-H1 [of size equal to NL_CONSTANT_SIZE defined below]. The oxygen atoms are extracted based on NL_CUTOFF and only the number of oxygen atoms equal to NL_CONSTANT_SIZE are retained. Similarly, hydrogen atoms are extracted based on atom indices that share the residue numbers with the oxygen atom indices and only NL_CONSTANT_SIZE hydrogen atoms are retained. Consequently, the total number of elements in PIV is equal to (PINESATOMS-3)x(PINESATOMS-3)x0.5 + (PINESATOMS-3)xNL_CONSTANT_SIZE + (PINESATOMS-3)xNL_CONSTANT_SIZE.
SFACTOR=1.0,1.0,...,1.0 # Should be the same size as of nelements = (PINESATOMS-3)x(PINESATOMS-3)x0.5 + (PINESATOMS-3)xNL_CONSTANT_SIZE + (PINESATOMS-3)xNL_CONSTANT_SIZE 
SORT=1.0,1.0,...,1.0 # Should be the same size as of nelements = (PINESATOMS-3)x(PINESATOMS-3)x0.5 + (PINESATOMS-3)xNL_CONSTANT_SIZE + (PINESATOMS-3)xNL_CONSTANT_SIZE 
SWITCH1={RATIONAL R_0=0.56 MM=10 NN=5},...,{RATIONAL R_0=0.56 MM=10 NN=5} # Should be the same size as of nelements = (PINESATOMS-3)x(PINESATOMS-3)x0.5 + (PINESATOMS-3)xNL_CONSTANT_SIZE + (PINESATOMS-3)xNL_CONSTANT_SIZE
NL_CUTOFF=4.0 #  Neighbor list cutoff for each element in PIV. Use a small value but update the neighbor list frequently, which can be expensive. Alternatively, use a high value and update less frequently to save computational time.
. small will atleast nl_Constant_size atoms. Code has a buffer number of atoms closest to a given atom.. it is hardcoded to 10. we can make it a variable.
NL_STRIDE=1.0 #frequency to update neighbor list.
NL_SKIN=0.1 # Additional buffer cutoff to nl_cutoff. Neighborlist reconstruction is triggered if atoms move out of a radius equal to nlcutoff+skin.
NL_CONSTANT_SIZE = 10 # Number of oxygen and hydrogen atoms to retain in each solute-solvent block of PIV.
WIRITEPIVESTRAJ # Write PIV in each frame. Alternatively, use `print STRIDE=1 ARG=pines.* FILE=colvar.out`. Depending on switching function parameters, some of the values could be ~0.9999. If the value exceeds the precision, then these values will be stored as NAN after converting to integer. In that case colvar reports these values as zeros but WIRITEPIVESTRAJ could be used for reporting the exact values.
... PINES
```
