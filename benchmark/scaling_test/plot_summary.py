import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

outfolder = ''


plt.figure(figsize=(8, 6))
N = np.linspace(1e1, 1e7, int(1e5))
plt.loglog(N, 0.01 * N**1, '--', label=r"Time $\propto N$", linewidth=1)
plt.loglog(N, 0.01 * N*np.log(N), '--', label=r"Time $\propto N\cdot log(N)$", linewidth=1)
N = np.linspace(1e1, 1e5, int(1e5))
plt.loglog(N, 0.000005 * N**2, '--', label=r"Time $\propto N^2$", linewidth=1)

# ----------- Plot 1: FMM Time vs n -----------
df = pd.read_csv(outfolder + "fmm_summary_results.csv")
plt.loglog(10**df["n"], df["Total execution time"], marker='x', label="FMM Total execution time")


# ----------- Plot 2: Tree Time vs n -----------
df = pd.read_csv(outfolder + "tree_summary_results.csv")
plt.loglog(10**df["n"], df["Total execution time"], marker='^', label="Tree Total execution time")

# ----------- Plot 3: Tree Time vs n -----------
df = pd.read_csv(outfolder + "naive_summary_results.csv")
plt.loglog(10**df["n"], df["Total execution time"], marker='^', label="Naive Gravity Total execution time")

plt.xlabel("N")
plt.ylabel("Time (ms)")
plt.title("Execution Time vs N")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(outfolder + "time_vs_N.png")
plt.show()
