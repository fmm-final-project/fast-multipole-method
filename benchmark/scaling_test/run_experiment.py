import os
import subprocess
import pandas as pd
import numpy as np
import re

exe_name = "./naive_gravity"
outfolder = ''

outfile = outfolder + "force.csv"

def generate_input(datafile):
    with open("tree.in", "w") as f:
        f.write(f"""# G
1.0
# THETA
0.3
# MAX_PARTICLES_PER_CELL
10
# datafile
{datafile}
# outfile
{outfile}
""")

def run_fmm_and_capture_output():
    result = subprocess.run([exe_name], capture_output=True, text=True)
    return result.stdout

def parse_time(log_text):
    times = {}
    for key in ["Building time", "Expansion time", "Dual Tree Walk time", "Local Passing Down time", "Total execution time"]:
        match = re.search(rf"{key}:\s+([\d.]+) ms", log_text)
        if match:
            times[key] = float(match.group(1))
    return times

summary = []

for n in range(1, 1 + 6):
    generate_input(f'uniform_cube_3d_1e{n}.bin')
    print(f"Running n = {n}")
    log = run_fmm_and_capture_output()
    times = parse_time(log)

    if os.path.exists(outfile):
        os.remove(outfile)

    summary.append({
        "n": n,
        **times
    })

df_summary = pd.DataFrame(summary)
df_summary.to_csv(outfolder + "naive_summary_results.csv", index=False)
print(df_summary)