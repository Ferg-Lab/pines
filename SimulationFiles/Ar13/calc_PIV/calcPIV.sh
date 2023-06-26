#!/bin/bash

plumed driver --plumed plumed_calcPIV.dat --mf_xtc biased.xtc >& out.out

sed -i '/#END OF FRAME/d' PINES_representation_traj.dat

python3 genNPYfile.py
