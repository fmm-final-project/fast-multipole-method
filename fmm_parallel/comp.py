import pandas as pd
import numpy as np

#df1 = pd.read_csv("force_direct.csv", sep=",", header=None)
#df2 = pd.read_csv("force_tree2fmm_parallel.csv", sep=",", header=None)
df1 = np.fromfile("force_direct.bin", dtype=np.float64)
N = df1.shape[0] // 3
df1 = df1.reshape(N, 3)
df2 = np.fromfile("force_tree2fmm_parallel.bin", dtype=np.float64).reshape(N, 3)

error = ((df2 - df1) / df1)**2
error = np.sqrt(error.sum(axis=0)) / N
for i in range(3):
    print("{:.10e}".format(error[i]))
