#plumed driver --plumed plumed_C45.dat --mf_pdb ../SimulationFiles/C45/md.pdb >& out_c45.out
plumed driver --natoms 100000  --parse-only --kt 2.49 --plumed plumed_C45.dat >& out_c45.out
