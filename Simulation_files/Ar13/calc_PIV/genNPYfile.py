import numpy as np
import pandas as pd
PIVs = pd.read_csv("PIV_representation_traj.dat", delim_whitespace="tab", header=None)
np.save("Ar13_PIVs.npy", PIVs)
