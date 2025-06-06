import pandas as pd
import numpy as np

df1 = pd.read_csv("force_direct.csv", sep=",", header=None)
df2 = pd.read_csv("force_tree2fmm_L2.csv", sep=",", header=None)

error = ((df2 - df1) / df1)**2
error = np.sqrt(error.sum()) / (len(df1[0]))
for i in range(3):
    print("{:.10e}".format(error[i]))
