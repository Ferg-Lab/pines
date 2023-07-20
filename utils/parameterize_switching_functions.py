#module load python/anaconda-2021.05
#source activate torch-ml

import torch
import torch.nn as nn
import numpy as np
import math
import re, shutil, tempfile
import sys
import os
import mdtraj as md
from sklearn.decomposition import PCA
import subprocess
from scipy.signal import savgol_filter

#inputPIV = np.load("nSFP_NaClpiv61.npy")
if len(sys.argv) < 2:
        print("COMMAND LINE ARGS:\n")
        print("1. Path to xtc file\n")
        print("2. Path to gro file\n")
        print("Remember to update atom list in script\n")
        sys.exit("Too few command line args")
dcd_file = sys.argv[2]
pdb_file = sys.argv[1]

#filepath = '/work2/03273/tg825722/shared-folder-siva/doe/hydrolysis/mlcolvar'
trj = md.load_dcd(dcd_file, top=pdb_file)


#chignolin alpha carbons
#atom_list = ["index 5", "index 26", "index 47", "index 67", "index 73", "index 88", "index 102", "index 109", "index 123", "index 147", "name OW"]

#chignolin hydrogen bonds
#atom_list = ["index ", "name OW"]


#atom_list = ["index 20", "index 17", "index 16", "index 54", "index 52", "index 49", "index 214", "index 212", "index 209", "index 180", "index 177", "index 176"]
#DNA hairpin
#atom_list = ["index 16", "index 17", "index 20", "index 49", "index 52", "index 54", "index 79", "index 80", "index 83", "index 112", "index 115", "index 144", "index 147", "index 176", "index 177", "index 180", "index 209", "index 212", "index 214", "name OW"]

# hydrolysis
atom_name_list = ['O1', 'O2', r'P$_{\alpha}$',  r'P$_{\beta}$', r'P$_{\gamma}$', r'O$_{\beta_{1}}$', r'O$_{\beta_{2}}$', r'O$_{\beta_{3}}$', r'O$_{\gamma_{1}}$', r'O$_{\gamma_{2}}$', r'O$_{\gamma_{3}}$']
atom_list = [3893, 3894, 13505, 13509, 13513, 13510, 13511, 13512, 13514, 13515, 13516]
atom_list = ['index ' + str(x-1) for x in atom_list]

def find_first_and_last_peak(rdf):
    rdf_i=0
    for i in range(len(rdf)):
        if i == 0:
            first_peak_index = i 
        elif i != 0 and rdf[i] >= rdf[i - 1]:
            first_peak_index = i
        else:
            if rdf_i >=10:
                break
            else:
                rdf_i+=1

    # Find the first peak from the right
    reverse_rdf = rdf[::-1]
    rdf_i=0
    for i in range(len(reverse_rdf)):
        if i == 0:
            last_peak_index = i 
        elif i != 0 and reverse_rdf[i] >= reverse_rdf[i - 1]:
            last_peak_index = i
        else:
            if rdf_i >=10:
                break
            else:
                rdf_i+=1

    return first_peak_index, len(rdf)-last_peak_index



interactions = []
for i in range(len(atom_list)-1):
    for j in range(i+1,len(atom_list)):
        interactions.append(atom_list[i]+" - "+atom_list[j])

starts = []
ends = []
r0s = []
for i in range(len(atom_list)-1):
    for j in range(i+1,len(atom_list)):
        atom_pairs = trj.top.select_pairs(atom_list[i], atom_list[j])
        dist = md.compute_distances(trj, atom_pairs)
        bins, counts = md.compute_rdf(trj, atom_pairs, (dist.min(), dist.max()), n_bins=100)
        #bins, counts = md.compute_rdf(trj, atom_pairs, (0.1, 2.5), n_bins=1000)
        
        smoothed_data = savgol_filter(counts, window_length=10, polyorder=2, mode='nearest') 
        pk1, pk2 = find_first_and_last_peak(smoothed_data) 


        r0 = (bins[pk1] + bins[pk2])/2.0
        starts.append(bins[pk1])
        ends.append(bins[pk2])
        r0s.append(round(r0,3))

def rat_switch(r, r0, d0, n, m):
    s = (1 - ((r-d0)/r0)**n)/(1 - ((r-d0)/r0)**m)
    return s

vec_switch = np.vectorize(rat_switch)

n_vals = np.zeros(len(interactions))
for n in range(1,100):
    m = 2*n

    closest_val = vec_switch(starts, r0s, 0, n, m) #input peak 1 start as r
    pk2_end = vec_switch(ends, r0s, 0, n, m) #input peak 2 end as r
    for i in range(len(starts)):
        if closest_val[i] > 0.9 and pk2_end[i] < 0.1:
            if n_vals[i] == 0:
                n_vals[i] = n

final_params = []
for i in range(len(interactions)):
    final_params.append([r0s[i], n_vals[i], interactions[i]])

r0s = np.asarray(r0s)
n_vals = np.asarray(n_vals)
string = ''

for i in range(len(r0s)):
    string += "SWITCH%d={RATIONAL R_0=%.3f MM=%d NN=%d}\n" %(i+1, r0s[i], 2*n_vals[i], n_vals[i])
    #print("SWITCH%d={RATIONAL R_0=%.3f MM=%d NN=%d}" %(i+1, r0s[i], 2*n_vals[i], n_vals[i]))

#print(string)

f = open("switching_function_parameters.dat", "w")
f.write(string)
#f.write("\n\n\n")
#for i in interactions:
#    f.write(str(i) + "\n")
f.close()
