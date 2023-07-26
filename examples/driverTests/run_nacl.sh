#plumed driver --plumed plumed_NaCl.dat --mf_pdb ../SimulationFiles/NaCl/md.pdb >& out_nacl.out
plumed driver --natoms 100000  --parse-only --kt 2.49 --plumed plumed_NaCl.dat >& out_nacl.out
