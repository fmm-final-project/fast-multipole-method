import os
import subprocess
import pandas as pd
import numpy as np
import re

thetas = np.linspace(0.2, 1.0, 17)
exe_name = "./tree"
outfolder = 'tree_output/'
datafile = "datafile/uniform_cube_3d_1e6.bin"
direct_result = "direct_result/force_direct_parallel.csv"

def generate_input(theta, output_file):
    with open("tree.in", "w") as f:
        f.write(f"""# G
1.0
# THETA
{theta}
# MAX_PARTICLES_PER_CELL
10
# datafile
{datafile}
# outfile
{output_file}
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

def compute_error(reference_csv, test_csv):
    df1 = pd.read_csv(reference_csv, sep=",", header=None)
    df2 = pd.read_csv(test_csv, sep=",", header=None)
    error = ((df2 - df1) / df1) ** 2
    error = np.sqrt(error.sum()) / len(df1[0])
    return error.values  # [Ex, Ey, Ez]

summary = []

for theta in thetas:
    outfile = outfolder + f"force_theta_{theta:.5f}.csv"
    generate_input(theta, outfile)
    print(f"Running theta = {theta}")
    log = run_fmm_and_capture_output()
    times = parse_time(log)
    errors = compute_error(direct_result, outfile)

    if os.path.exists(outfile):
        os.remove(outfile)

    summary.append({
        "theta": theta,
        "Ex_error": errors[0],
        "Ey_error": errors[1],
        "Ez_error": errors[2],
        **times
    })

df_summary = pd.DataFrame(summary)
df_summary.to_csv(outfolder + "/summary_results.csv", index=False)
print(df_summary)
