import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

outfolder = 'tree_output/'
df = pd.read_csv(outfolder + "summary_results.csv")
name = "Tree"

time_unit = 1000

# ----------- Plot 1: Time vs nThreads -----------
plt.figure(figsize=(8, 6))
plt.plot(df["nThreads"][1:], df["Building time"][1:]/time_unit, marker='o', label="Building time")
plt.plot(df["nThreads"][1:], df["Expansion time"][1:]/time_unit, marker='s', label="Expansion time")
# plt.plot(df["nThreads"][1:], df["Dual Tree Walk time"][1:]/time_unit, marker='^', label="Dual Tree Walk time")
# plt.plot(df["nThreads"][1:], df["Local Passing Down time"][1:]/time_unit, marker='d', label="Local Passing Down time")
plt.plot(df["nThreads"][1:], df["Total execution time"][1:]/time_unit, marker='x', label="Total execution time")
plt.plot(df["nThreads"][1:], df["Total execution time"][0]/time_unit/df["nThreads"][1:], '--', label="Ideal execution Time")


plt.xlabel("nThreads")
plt.ylabel("Time (s)")
plt.title(name + " Execution Time vs nThreads")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(outfolder + "time_vs_nThreads.png")
plt.show()

# ----------- Plot 2: SpeedUp vs nThreads -----------
plt.figure(figsize=(8, 6))
one_core_time = df["Total execution time"][0]
plt.plot(df["nThreads"][1:], one_core_time/df["Total execution time"][1:], marker='o', label="SpeedUp")
plt.plot(df["nThreads"][1:], df["nThreads"][1:], '--', label='Ideal (Linear SpeedUp)')
plt.xlabel("nThreads")
plt.ylabel("SpeedUp")
plt.title(name + " SpeedUp vs nThreads")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(outfolder + "SpeedUp_vs_nThreads.png")
plt.show()
