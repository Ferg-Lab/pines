import numpy as np
import pandas as pd
PIVs = pd.read_csv("PINES_representation_traj.dat", delim_whitespace="tab", header=None)
np.save("NaCl_PIVs.npy", PIVs)
