#plumed driver --plumed plumed_Ar13.dat --mf_pdb ../SimulationFiles/Ar13/md.pdb >& out_argon.out
plumed driver --natoms 100000  --parse-only --kt 1 --plumed plumed_Ar13.dat >& out_argon.out
