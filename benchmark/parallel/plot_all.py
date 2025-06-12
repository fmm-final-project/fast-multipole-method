import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


outfolder = ['naive_output/', 'tree_output/', 'fmm_output/']
name = ["Naive Gravity", "Tree", "FMM"]
marker = ["o", "s", "^"]

time_unit = 1000

# ----------- Plot 1: Time vs nThreads -----------
plt.figure(figsize=(8, 6))
for i in range(3):
    df = pd.read_csv(outfolder[i] + "summary_results.csv")
    plt.plot(df["nThreads"][1:], df["Total execution time"][1:]/time_unit, marker=marker[i], label=name[i]+" Total execution time")

plt.xlabel("nThreads")
plt.ylabel("Time (s)")
plt.yscale("log")
plt.title("Execution Time vs nThreads")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("ALL_time_vs_nThreads.png")
plt.show()

# ----------- Plot 2: SpeedUp vs nThreads -----------
plt.figure(figsize=(8, 6))
for i in range(3):
    df = pd.read_csv(outfolder[i] + "summary_results.csv")
    one_core_time = df["Total execution time"][0]
    plt.plot(df["nThreads"][1:], one_core_time/df["Total execution time"][1:], marker=marker[i], label=name[i]+" SpeedUp")
plt.plot(df["nThreads"][1:], df["nThreads"][1:], '--', label='Ideal (Linear SpeedUp)')
plt.xlabel("nThreads")
plt.ylabel("SpeedUp")
plt.title("SpeedUp vs nThreads")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("ALL_SpeedUp_vs_nThreads.png")
plt.show()
