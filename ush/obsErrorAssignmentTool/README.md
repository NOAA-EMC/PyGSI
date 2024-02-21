## 2/20/2024
## Azadeh Gholoubi
# Python Tool for Observation Error Assignment

## Overview
This tool automates the estimation of observation error parameters. It consists of a shell script, `run_error_model_estimate.sh`, which invokes the `ObsErrorAssignment.py` Python script to perform the following tasks:

1. Compute mean cloud amount and first guess standard deviation.
2. Generate four plots:
   - Standard deviation of FG departures in each channel as a function of mean cloud amount: (clw_guess + clw_obs)/2
   - Number of observations in each defined bin of mean cloud amount
   - Histogram of un-normalised and also normalised FG departures in each channel
   - Histogram of Errors 
3. Estimate error parameters (`error_cld`, `error_clr`, `clw_clr`, and `clw_cld`) and create CSV files for each channel.

## Usage
Before running the script, the user needs to specify certain parameters in `run_error_model_estimate.sh`:

- `config_path`: Directory for input config .nc files
- `output_path`: Path to output files
- `sensor`: Sensor name
- `Channels`: Channel numbers (enclosed in quotes, separated by commas, e.g., "1,2,3,5")
- `bin_size`: Size of bin for plotting (default=0.05)
- `bindir`: Output folder name
- `qc_flag`: QC flag for filtering (default=0)
- `PyGSI`: Path to the PyGSI branch on Hera or Orion to load the proper environment

To view all the required inputs for this tool, use the following command:

```bash
python ObsErrorAssignment.py -h



