# Performance under varying `theta` parameters

This repository is used to benchmark the performance and accuracy of our FMM and tree implementation under varying `theta` parameters. Each test case is organized in its own folder.

Each test folder contains the following files:

- `run_experiment.py`: Automates the process of running the FMM or tree simulation for different `theta` values, collecting results such as execution time and force errors, and saving them to `summary_results.csv`.
- `plot_summary.py`: Reads `summary_results.csv` and generates plots for:
  - Force error (`Ex`, `Ey`, `Ez`) vs. `theta`
  - Runtime breakdown (build, expansion, dual tree walk, etc.) vs. `theta`
