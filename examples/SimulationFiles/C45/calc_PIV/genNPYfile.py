import numpy as np
import pandas as pd
C45_PIVs = pd.read_csv("PINES_representation_traj.dat", delim_whitespace="tab", header=None)
np.save("C45_PIVs.npy", C45_PIVs)
