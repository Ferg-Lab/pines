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
import pickle

#inputPIV = np.load("nSFP_NaClpiv61.npy")
if len(sys.argv) < 2:
        print("COMMAND LINE ARGS:\n")
        print("1. Path to xtc file\n")
        print("2. Path to gro file\n")
        print("Remember to update atom list in script\n")
        sys.exit("Too few command line args")
xtc_file = sys.argv[2]
pdb_file = sys.argv[1]
computeRDF = False
NL_CONSTANT_SIZE = 5

#filepath = '/work2/03273/tg825722/shared-folder-siva/doe/hydrolysis/mlcolvar'
#trj = md.load_dcd(dcd_file, top=pdb_file)
if computeRDF:
    trj = md.load_xtc(xtc_file, top=pdb_file)

# hydrolysis
#atom_name_list = ['O1', 'O2', r'P$_{\alpha}$',  r'P$_{\beta}$', r'P$_{\gamma}$', r'O$_{\beta_{1}}$', r'O$_{\beta_{2}}$', r'O$_{\beta_{3}}$', r'O$_{\gamma_{1}}$', r'O$_{\gamma_{2}}$', r'O$_{\gamma_{3}}$']
atom_list = [3893, 3894, 13505, 13509, 13513, 13510, 13511, 13512, 13514, 13515, 13516]
select_atom_list = ['index ' + str(x-1) for x in atom_list]

def find_first_and_last_peak(rdf, nbins_avg=6):
    rdf_i=0
    first_peak_index = 0
    for i in range(len(rdf)-3*nbins_avg):
        # compute average over N bins
        rdf_mean_Nbins = np.polyfit(np.arange(len(rdf[i:i+nbins_avg])), rdf[i:i+nbins_avg], 1)[0]
        # compute average over the next N bins
        rdf_mean_N_plus_bins = np.polyfit(np.arange(len(rdf[i+nbins_avg:i+2*nbins_avg])), rdf[i+nbins_avg:i+2*nbins_avg], 1)[0] 
        # compute average over the next next N bins
        rdf_mean_N_plusplus_bins = np.polyfit(np.arange(len(rdf[i+2*nbins_avg:i+3*nbins_avg])), rdf[i+2*nbins_avg:i+3*nbins_avg], 1)[0] 
        if rdf_mean_Nbins > 0:
            continue
        elif rdf_mean_N_plus_bins < 0:
            #print(f'rdf_mean_N_plus_bins:{rdf_mean_N_plus_bins}')
            first_peak_index = i
            if rdf_mean_N_plusplus_bins < 0: 
            # if no change in N continous bins
            #if rdf_i >=5:
                break
            #else:
            #rdf_i+=1

    # Find the second peak
    #second_peak_rdf = rdf[first_peak_index:]
    rdf_i=0
    second_peak_index=first_peak_index+1
    for i in range(first_peak_index+1, len(rdf)-3*nbins_avg):
        # compute average over N bins
        second_rdf_mean_Nbins = np.polyfit(np.arange(len(rdf[i:i+nbins_avg])), rdf[i:i+nbins_avg], 1)[0]
        # compute average over the next N bins
        second_rdf_mean_N_plus_bins = np.polyfit(np.arange(len(rdf[i+nbins_avg:i+2*nbins_avg])), rdf[i+nbins_avg:i+2*nbins_avg], 1)[0] 
        # compute average over the next next N bins
        second_rdf_mean_N_plusplus_bins = np.polyfit(np.arange(len(rdf[i+2*nbins_avg:i+3*nbins_avg])), rdf[i+2*nbins_avg:i+3*nbins_avg], 1)[0] 
        if second_rdf_mean_Nbins < 0:
            continue
        else:
            #print(f'reached at i:{i}, second_rdf_mean_Nbins:{second_rdf_mean_Nbins}, second_rdf_mean_N_plus_bins:{second_rdf_mean_N_plus_bins}, second_rdf_mean_N_plusplus_bins:{second_rdf_mean_N_plusplus_bins}')
            second_peak_index = i
            if second_rdf_mean_N_plus_bins < 0:
                break



        #if second_rdf_mean_Nbins > 0: # and second_rdf_mean_N_plusplus_bins >= second_rdf_mean_N_plus_bins:
        #    second_peak_index = i+1
        #elif second_rdf_mean_N_plus_bins < 0:
        #    print(f'reached at i:{i}, second_rdf_mean_Nbins:{second_rdf_mean_Nbins}, second_rdf_mean_N_plus_bins:{second_rdf_mean_N_plus_bins}')
        #    break

    #if second_peak_index == -1:
    #   second_peak_index = len(rdf)-1 
    #print(first_peak_index, second_peak_index)
    return first_peak_index, second_peak_index

#interactions = []
#for i in range(len(atom_list)-1):
#    for j in range(i+1,len(atom_list)):
#        interactions.append(atom_list[i]+" - "+atom_list[j])

starts = []
ends = []
starts_OW = []
ends_OW = []
starts_HW = []
ends_HW = []
r0s = []
OW_r0s = []
HW_r0s = []
count_ij = 0
count_i = 0
#for i in range(len(atom_list)-1):
#    for j in range(i+1,len(atom_list)):

for i in range(len(atom_list)):
    
    for j in range(len(atom_list)):

        if i <= j:
            continue

        if computeRDF:
            atom_pairs = trj.top.select_pairs(select_atom_list[i], select_atom_list[j])
            dist = md.compute_distances(trj, atom_pairs)
            bins, counts = md.compute_rdf(trj, atom_pairs, (dist.min(), dist.max()), n_bins=100)
            starts.append(dist.min())
            ends.append(dist.max())  
        else:
            with open('rdf_solute_solute.pkl', 'rb') as handle:
                data = pickle.load(handle)
                bins = data['rdf_solute_solute_pair'][count_ij]['bins']
                counts = data['rdf_solute_solute_pair'][count_ij]['counts'] 
                starts.append(data['starts'][count_ij])
                ends.append(data['ends'][count_ij])
        
        smoothed_data = savgol_filter(counts, window_length=20, polyorder=2, mode='nearest') 
        pk1, pk2 = find_first_and_last_peak(smoothed_data)

        #print(f'i:{i}, j:{j}, count_ij:{count_ij}')
        r0 = np.mean(bins[pk1:pk2]) # bins[np.argmin(np.abs(counts - np.mean(bins[pk1:pk2])))] #/2.0
        r0s.append(round(r0,3))

        count_ij += 1

for i in range(len(atom_list)):

    if computeRDF:
        atom_pairs = trj.topology.select_pairs(select_atom_list[i], 'name OT')
        dist_OW = md.compute_distances(trj, atom_pairs)
        bins_OW, counts_OW = md.compute_rdf(trj, atom_pairs, (0.05, 1.5), n_bins=100)
        starts_OW.append(0.05)
        ends_OW.append(2.0)  
    else:
        with open('rdf_solute_OW_2.5.pkl', 'rb') as handle:
            data = pickle.load(handle)
            bins_OW = data['rdf_solute_OW_pair'][count_i]['bins']
            counts_OW = data['rdf_solute_OW_pair'][count_i]['counts'] 
            starts_OW.append(data['starts_OW'][count_i])
            ends_OW.append(data['ends_OW'][count_i])

    smoothed_data = savgol_filter(counts_OW, window_length=20, polyorder=2, mode='nearest') 
    pk1, pk2 = find_first_and_last_peak(smoothed_data) 

    OW_r0 = np.mean(bins_OW[pk1:pk2])
    OW_r0s.append(round(OW_r0,3))

for i in range(len(atom_list)):
    
    if computeRDF:
        atom_pairs = trj.topology.select_pairs(select_atom_list[i], 'name HT')
        dist_HW = md.compute_distances(trj, atom_pairs)
        bins_HW, counts_HW = md.compute_rdf(trj, atom_pairs, (0.05, 1.5), n_bins=100)
        starts_HW.append(0.05)
        ends_HW.append(2.0)  
    else:
        with open('rdf_solute_HW_2.5.pkl', 'rb') as handle:
            data = pickle.load(handle)
            bins_HW = data['rdf_solute_HW_pair'][count_i]['bins']
            counts_HW = data['rdf_solute_HW_pair'][count_i]['counts'] 
            starts_HW.append(data['starts_HW'][count_i])
            ends_HW.append(data['ends_HW'][count_i])

    smoothed_data = savgol_filter(counts_HW, window_length=20, polyorder=2, mode='nearest') 
    pk1, pk2 = find_first_and_last_peak(smoothed_data) 

    HW_r0 = np.mean(bins_HW[pk1:pk2])
    HW_r0s.append(round(HW_r0,3))

    count_i += 1
         
def rat_switch(r, r0, d0, n, m):
    s = (1 - ((r-d0)/r0)**n)/(1 - ((r-d0)/r0)**m)
    return s

#vec_switch = np.vectorize(rat_switch)
#n_vals = np.zeros(len(interactions))
n_vals = np.ones(len(r0s))*6
for r0 in range(len(r0s)):  
    for n in range(1,100):
        m = 2*n

    #closest_val = vec_switch(starts, r0s, 0, n, m) #input peak 1 start as r
    #pk2_end = vec_switch(ends, r0s, 0, n, m) #input peak 2 end as r

        r = np.arange(starts[r0],ends[r0],0.05)
        s = rat_switch(r, r0s[r0], 0, n, m) 
    #for i in range(len(starts)):
        #if closest_val[i] > 0.9 and pk2_end[i] < 0.1:
        #    if n_vals[i] == 0:
        #        n_vals[i] = n
        if s.max() > 0.9 and s.min() < 0.1:
            n_vals[r0] = n
            break
    #print(n_vals[r0])

OW_n_vals = np.ones(len(OW_r0s))*6
for r0 in range(len(OW_r0s)):
    for n in range(1,100):
        m = 2*n
        r = np.arange(starts_OW[r0],ends_OW[r0],0.05)
        s = rat_switch(r, OW_r0s[r0], 0, n, m) 
        if s.max() > 0.9 and s.min() < 0.1:
            OW_n_vals[r0] = n
            break

HW_n_vals = np.ones(len(HW_r0s))*6
for r0 in range(len(HW_r0s)):
    for n in range(1,100):
        m = 2*n
        r = np.arange(starts_HW[r0],ends_HW[r0],0.05)
        s = rat_switch(r, HW_r0s[r0], 0, n, m) 
        if s.max() > 0.9 and s.min() < 0.1:
            HW_n_vals[r0] = n
            break
       
        

r0s = np.asarray(r0s)
n_vals = np.asarray(n_vals)
OW_r0s = np.asarray(OW_r0s)
OW_n_vals = np.asarray(OW_n_vals)
HW_r0s = np.asarray(HW_r0s)
HW_n_vals = np.asarray(HW_n_vals)

with open('final_switching_function_parameters.pkl', 'wb') as handle:
    pickle.dump({'r0s':r0s, 'n':n_vals, 'OW_r0s':OW_r0s, 'OW_n':OW_n_vals, 'HW_r0s': HW_r0s, 'HW_n': HW_n_vals}, handle)
