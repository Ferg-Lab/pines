#module load python/anaconda-2021.05
#source activate torch-ml

import torch
import torch.nn as nn
import numpy as np
import math
import re, shutil, tempfile
import sys
import os
from sklearn.decomposition import PCA
import subprocess

def sed_inplace(filename, pattern, repl):
    '''
    Perform the pure-Python equivalent of in-place `sed` substitution: e.g.,
    `sed -i -e 's/'${pattern}'/'${repl}' "${filename}"`.
    '''
    # For efficiency, precompile the passed regular expression.
    pattern_compiled = re.compile(pattern)

    # For portability, NamedTemporaryFile() defaults to mode "w+b" (i.e., binary
    # writing with updating). This is usually a good thing. In this case,
    # however, binary writing imposes non-trivial encoding constraints trivially
    # resolved by switching to text writing. Let's do that.
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_file:
        with open(filename) as src_file:
            for line in src_file:
                tmp_file.write(pattern_compiled.sub(repl, line))

    # Overwrite the original file with the munged temporary file in a
    # manner preserving file attributes (e.g., permissions).
    shutil.copystat(filename, tmp_file.name)
    shutil.move(tmp_file.name, filename)


#inputPIV = np.load("nSFP_NaClpiv61.npy")
if len(sys.argv) < 3:
        print("COMMAND LINE ARGS:\n")
        print("1. Path to checkpoint file\n")
        print("2. Grid mins: One per PC and comma separated (i.e. -3.5,-3.5)\n")
        print("3. Grid maxes: One per PC and comma separated (i.e. 3.5,3.5)\n")
        sys.exit("Too few command line args")
grid_mins = sys.argv[2]
grid_maxes = sys.argv[3]
latent = np.load("LS.npy")
eps = 0.00001
chkpnt_file_path = sys.argv[1]
cwd = os.getcwd()
shutil.copy(cwd+"/plumed_template.dat", cwd+"/plumed_parameterized.dat")
plumed_file = "plumed_parameterized.dat"
model = torch.load(chkpnt_file_path, map_location=torch.device('cpu'))

weights0= model['state_dict']['encoder.2.weight']
np_weights0 = weights0.numpy()
np_weights0_flt = np_weights0.flatten()
np_weights0_arrstr = np.char.mod('%f', np_weights0_flt)
np_weights0_str = ",".join(np_weights0_arrstr)
np_weights0_str = "GAMMAS0=" + np_weights0_str
sed_inplace(plumed_file, 'GAMMAS0=.*', np_weights0_str)

bias0= model['state_dict']['encoder.2.bias']
np_bias0 = bias0.numpy()
np_bias0_flt = np_bias0.flatten()
np_bias0_arrstr = np.char.mod('%f', np_bias0_flt)
np_bias0_str = ",".join(np_bias0_arrstr)
np_bias0_str = "BETAS0=" + np_bias0_str
sed_inplace(plumed_file, 'BETAS0=.*', np_bias0_str)


weights2= model['state_dict']['encoder.2.running_mean']
np_weights2 = weights2.numpy()
np_weights2_flt = np_weights2.flatten()
np_weights2_arrstr = np.char.mod('%f', np_weights2_flt)
np_weights2_str = ",".join(np_weights2_arrstr)
np_weights2_str = "EXPECTATIONS0=" + np_weights2_str
sed_inplace(plumed_file, 'EXPECTATIONS0=.*', np_weights2_str)

bias2= model['state_dict']['encoder.2.running_var']
np_bias2 = bias2.numpy()
np_bias2 = np.sqrt(np_bias2 + eps)
np_bias2_flt = np_bias2.flatten()
np_bias2_arrstr = np.char.mod('%f', np_bias2_flt)
np_bias2_str = ",".join(np_bias2_arrstr)
np_bias2_str = "VARIANCES0=" + np_bias2_str
sed_inplace(plumed_file, 'VARIANCES0=.*', np_bias2_str)

weights0= model['state_dict']['encoder.5.weight']
np_weights0 = weights0.numpy()
np_weights0_flt = np_weights0.flatten()
np_weights0_arrstr = np.char.mod('%f', np_weights0_flt)
np_weights0_str = ",".join(np_weights0_arrstr)
np_weights0_str = "GAMMAS1=" + np_weights0_str
sed_inplace(plumed_file, 'GAMMAS1=.*', np_weights0_str)

bias0= model['state_dict']['encoder.5.bias']
np_bias0 = bias0.numpy()
np_bias0_flt = np_bias0.flatten()
np_bias0_arrstr = np.char.mod('%f', np_bias0_flt)
np_bias0_str = ",".join(np_bias0_arrstr)
np_bias0_str = "BETAS1=" + np_bias0_str
sed_inplace(plumed_file, 'BETAS1=.*', np_bias0_str)

weights2= model['state_dict']['encoder.5.running_mean']
np_weights2 = weights2.numpy()
np_weights2_flt = np_weights2.flatten()
np_weights2_arrstr = np.char.mod('%f', np_weights2_flt)
np_weights2_str = ",".join(np_weights2_arrstr)
np_weights2_str = "EXPECTATIONS1=" + np_weights2_str
sed_inplace(plumed_file, 'EXPECTATIONS1=.*', np_weights2_str)

bias2= model['state_dict']['encoder.5.running_var']
np_bias2 = bias2.numpy()
np_bias2 = np.sqrt(np_bias2 + eps)
np_bias2_flt = np_bias2.flatten()
np_bias2_arrstr = np.char.mod('%f', np_bias2_flt)
np_bias2_str = ",".join(np_bias2_arrstr)
np_bias2_str = "VARIANCES1=" + np_bias2_str
sed_inplace(plumed_file, 'VARIANCES1=.*', np_bias2_str)

weights0= model['state_dict']['encoder.8.weight']
np_weights0 = weights0.numpy()
np_weights0_flt = np_weights0.flatten()
np_weights0_arrstr = np.char.mod('%f', np_weights0_flt)
np_weights0_str = ",".join(np_weights0_arrstr)
np_weights0_str = "GAMMAS2=" + np_weights0_str
sed_inplace(plumed_file, 'GAMMAS2=.*', np_weights0_str)

bias0= model['state_dict']['encoder.8.bias']
np_bias0 = bias0.numpy()
np_bias0_flt = np_bias0.flatten()
np_bias0_arrstr = np.char.mod('%f', np_bias0_flt)
np_bias0_str = ",".join(np_bias0_arrstr)
np_bias0_str = "BETAS2=" + np_bias0_str
sed_inplace(plumed_file, 'BETAS2=.*', np_bias0_str)


weights2= model['state_dict']['encoder.8.running_mean']
np_weights2 = weights2.numpy()
np_weights2_flt = np_weights2.flatten()
np_weights2_arrstr = np.char.mod('%f', np_weights2_flt)
np_weights2_str = ",".join(np_weights2_arrstr)
np_weights2_str = "EXPECTATIONS2=" + np_weights2_str
sed_inplace(plumed_file, 'EXPECTATIONS2=.*', np_weights2_str)

bias2= model['state_dict']['encoder.8.running_var']
np_bias2 = bias2.numpy()
np_bias2 = np.sqrt(np_bias2 + eps)
np_bias2_flt = np_bias2.flatten()
np_bias2_arrstr = np.char.mod('%f', np_bias2_flt)
np_bias2_str = ",".join(np_bias2_arrstr)
np_bias2_str = "VARIANCES2=" + np_bias2_str
sed_inplace(plumed_file, 'VARIANCES2=.*', np_bias2_str)

weights0= model['state_dict']['encoder.11.weight']
np_weights0 = weights0.numpy()
np_weights0_flt = np_weights0.flatten()
np_weights0_arrstr = np.char.mod('%f', np_weights0_flt)
np_weights0_str = ",".join(np_weights0_arrstr)
np_weights0_str = "GAMMAS3=" + np_weights0_str
sed_inplace(plumed_file, 'GAMMAS3=.*', np_weights0_str)

bias0= model['state_dict']['encoder.11.bias']
np_bias0 = bias0.numpy()
np_bias0_flt = np_bias0.flatten()
np_bias0_arrstr = np.char.mod('%f', np_bias0_flt)
np_bias0_str = ",".join(np_bias0_arrstr)
np_bias0_str = "BETAS3=" + np_bias0_str
sed_inplace(plumed_file, 'BETAS3=.*', np_bias0_str)


weights2= model['state_dict']['encoder.11.running_mean']
np_weights2 = weights2.numpy()
np_weights2_flt = np_weights2.flatten()
np_weights2_arrstr = np.char.mod('%f', np_weights2_flt)
np_weights2_str = ",".join(np_weights2_arrstr)
np_weights2_str = "EXPECTATIONS3=" + np_weights2_str
sed_inplace(plumed_file, 'EXPECTATIONS3=.*', np_weights2_str)

bias2= model['state_dict']['encoder.11.running_var']
np_bias2 = bias2.numpy()
np_bias2 = np.sqrt(np_bias2 + eps)
np_bias2_flt = np_bias2.flatten()
np_bias2_arrstr = np.char.mod('%f', np_bias2_flt)
np_bias2_str = ",".join(np_bias2_arrstr)
np_bias2_str = "VARIANCES3=" + np_bias2_str
sed_inplace(plumed_file, 'VARIANCES3=.*', np_bias2_str)


weights0= model['state_dict']['encoder.0.weight']
np_weights0 = weights0.numpy()
np_weights0_flt = np_weights0.flatten()
np_weights0_arrstr = np.char.mod('%f', np_weights0_flt)
np_weights0_str = ",".join(np_weights0_arrstr)
np_weights0_str = "WEIGHTS0=" + np_weights0_str
sed_inplace(plumed_file, 'WEIGHTS0=.*', np_weights0_str)

bias0= model['state_dict']['encoder.0.bias']
np_bias0 = bias0.numpy()
np_bias0_flt = np_bias0.flatten()
np_bias0_arrstr = np.char.mod('%f', np_bias0_flt)
np_bias0_str = ",".join(np_bias0_arrstr)
np_bias0_str = "BIASES0=" + np_bias0_str
sed_inplace(plumed_file, 'BIASES0=.*', np_bias0_str)


weights2= model['state_dict']['encoder.3.weight']
np_weights2 = weights2.numpy()
np_weights2_flt = np_weights2.flatten()
np_weights2_arrstr = np.char.mod('%f', np_weights2_flt)
np_weights2_str = ",".join(np_weights2_arrstr)
np_weights2_str = "WEIGHTS1=" + np_weights2_str
sed_inplace(plumed_file, 'WEIGHTS1=.*', np_weights2_str)

bias2= model['state_dict']['encoder.3.bias']
np_bias2 = bias2.numpy()
np_bias2_flt = np_bias2.flatten()
np_bias2_arrstr = np.char.mod('%f', np_bias2_flt)
np_bias2_str = ",".join(np_bias2_arrstr)
np_bias2_str = "BIASES1=" + np_bias2_str
sed_inplace(plumed_file, 'BIASES1=.*', np_bias2_str)

weights0= model['state_dict']['encoder.6.weight']
np_weights0 = weights0.numpy()
np_weights0_flt = np_weights0.flatten()
np_weights0_arrstr = np.char.mod('%f', np_weights0_flt)
np_weights0_str = ",".join(np_weights0_arrstr)
np_weights0_str = "WEIGHTS2=" + np_weights0_str
sed_inplace(plumed_file, 'WEIGHTS2=.*', np_weights0_str)

bias0= model['state_dict']['encoder.6.bias']
np_bias0 = bias0.numpy()
np_bias0_flt = np_bias0.flatten()
np_bias0_arrstr = np.char.mod('%f', np_bias0_flt)
np_bias0_str = ",".join(np_bias0_arrstr)
np_bias0_str = "BIASES2=" + np_bias0_str
sed_inplace(plumed_file, 'BIASES2=.*', np_bias0_str)


weights2= model['state_dict']['encoder.9.weight']
np_weights2 = weights2.numpy()
np_weights2_flt = np_weights2.flatten()
np_weights2_arrstr = np.char.mod('%f', np_weights2_flt)
np_weights2_str = ",".join(np_weights2_arrstr)
np_weights2_str = "WEIGHTS3=" + np_weights2_str
sed_inplace(plumed_file, 'WEIGHTS3=.*', np_weights2_str)

bias2= model['state_dict']['encoder.9.bias']
np_bias2 = bias2.numpy()
np_bias2_flt = np_bias2.flatten()
np_bias2_arrstr = np.char.mod('%f', np_bias2_flt)
np_bias2_str = ",".join(np_bias2_arrstr)
np_bias2_str = "BIASES3=" + np_bias2_str
sed_inplace(plumed_file, 'BIASES3=.*', np_bias2_str)

sed_inplace(plumed_file, 'GRID_MIN=.*', 'GRID_MIN='+grid_mins+' GRID_MAX='+grid_maxes+' GRID_BIN=5000,5000,5000')
sys.exit("\n\t\tSUCCESS!\n\nplumed_parameterized.dat file created!\n")
