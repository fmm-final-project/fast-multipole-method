import pandas as pd
import matplotlib.pyplot as plt

outfolder = 'tree_output/'
df = pd.read_csv(outfolder + "summary_results.csv")
name = "Tree"
time_unit = 1000

# ----------- Plot 1: Error vs Theta -----------
plt.figure(figsize=(8, 6))
plt.plot(df["theta"], df["Ex_error"], marker='o', label="Ex_error")
plt.plot(df["theta"], df["Ey_error"], marker='s', label="Ey_error")
plt.plot(df["theta"], df["Ez_error"], marker='^', label="Ez_error")
plt.xlabel("Theta")
plt.ylabel("Relative Error")
plt.title(name+" Force Error vs Theta")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(outfolder + "error_vs_theta.png")
plt.show()

# ----------- Plot 2: Time vs Theta -----------
plt.figure(figsize=(8, 6))
plt.plot(df["theta"], df["Building time"]/time_unit, marker='o', label="Building time")
plt.plot(df["theta"], df["Expansion time"]/time_unit, marker='s', label="Expansion time")
# plt.plot(df["theta"], df["Dual Tree Walk time"]/time_unit, marker='^', label="Dual Tree Walk time")
# plt.plot(df["theta"], df["Local Passing Down time"]/time_unit, marker='d', label="Local Passing Down time")
plt.plot(df["theta"], df["Total execution time"]/time_unit, marker='x', label="Total execution time")
plt.xlabel("Theta")
plt.ylabel("Time (s)")
plt.title(name+" Execution Time vs Theta")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(outfolder + "time_vs_theta.png")
plt.show()
