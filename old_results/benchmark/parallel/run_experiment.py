import os
import subprocess
import pandas as pd
import numpy as np
import re

exe_name = "./tree2fmm_parallel"
outfolder = 'fmm_output/'
datafile = "datafile/plummer_sphere_3d_1e6.bin"
Max_nThreads = 24

outfile = outfolder + "force.csv"

def generate_input(nThreads):
    with open("tree.in", "w") as f:
        f.write(f"""# G
1.0
# THETA
0.2
# MAX_PARTICLES_PER_CELL
10
# datafile
{datafile}
# outfile
{outfile}
# NUM_OF_THREADS
{nThreads}
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

for nThreads in range(1, 1 + Max_nThreads):
    generate_input(nThreads)
    print(f"Running nThreads = {nThreads}")
    log = run_fmm_and_capture_output()
    times = parse_time(log)

    if os.path.exists(outfile):
        os.remove(outfile)

    summary.append({
        "nThreads": nThreads,
        **times
    })

df_summary = pd.DataFrame(summary)
df_summary.to_csv(outfolder + "/summary_results_fmm.csv", index=False)
print(df_summary)