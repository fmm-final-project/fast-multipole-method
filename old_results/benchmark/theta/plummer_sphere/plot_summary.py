import pandas as pd
import matplotlib.pyplot as plt

outfolder = 'tree_output/'
df = pd.read_csv(outfolder + "summary_results_tree.csv")

# ----------- Plot 1: Error vs Theta -----------
plt.figure(figsize=(8, 6))
plt.plot(df["theta"], df["Ex_error"], marker='o', label="Ex_error")
plt.plot(df["theta"], df["Ey_error"], marker='s', label="Ey_error")
plt.plot(df["theta"], df["Ez_error"], marker='^', label="Ez_error")
plt.xlabel("Theta")
plt.ylabel("Relative Error")
plt.title("Force Error vs Theta")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(outfolder + "error_vs_theta.png")
plt.show()

# ----------- Plot 2: Time vs Theta -----------
plt.figure(figsize=(8, 6))
plt.plot(df["theta"], df["Building time"], marker='o', label="Building time")
plt.plot(df["theta"], df["Expansion time"], marker='s', label="Expansion time")
# plt.plot(df["theta"], df["Dual Tree Walk time"], marker='^', label="Dual Tree Walk time")
# plt.plot(df["theta"], df["Local Passing Down time"], marker='d', label="Local Passing Down time")
plt.plot(df["theta"], df["Total execution time"], marker='x', label="Total execution time")
plt.xlabel("Theta")
plt.ylabel("Time (ms)")
plt.title("Execution Time vs Theta")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(outfolder + "time_vs_theta.png")
plt.show()
