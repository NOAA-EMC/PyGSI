### November 2024
### Azadeh Gholoubi
# Spatial and Temporal Analysis Tool for Satellite Observation Data

## Overview
**Purpose**:  This tool performs spatial and temporal analysis for satellite observation data, allowing users to create customizable grids, filter data by time and region, and generate statistical and summary plots.

### Key Functionalities:
- Grid-based Data Summaries: Creates spatial grids for data aggregation.
- Data Filtering: Processes data across specified time frames and geographical regions.
- Visualization: Generates evaluation plots for different data attributes and regions

The `SpatialTemporalStats` class is central to this tool, with methods for creating grids, reading observational data, filtering, plotting, and producing summary statistics.

### Important Methods of the SpatialTemporalStats Class
- `generate_grid(resolution=1)`: Generates a spatial grid with specified resolution (default: 1x1 degree).
- `read_obs_values()`:  Reads and filters observational data from NetCDF files, performs spatial averaging, and returns averaged values.
- `plot_obs()`: Plots observation data on a map, with options for different regions and grid sizes.
- `list_variable_names(file_path)`: Lists variable names available in a specified NetCDF file.
- `make_summary_plots()`: Generates scatter plots for counts, means, and standard deviations of observational data.

## Requirements
User need to load EVA environment when working on Hera, use the following commands:
```
cd GDASApp/modulefiles/
module load EVA/hera

```

## Usage
To get started, run the following command to see all available options and argument formats:
```python
python SpatialTemporalStats.py -h
```
This command will display detailed information on how to input your settings. Key parameters include:

- input: Path to input data files
- output: Path for saving the results
- sensor: Satellite sensor name (e.g., "atms_n20")
- var: Variable to analyze (e.g., "Obs_Minus_Forecast_adjusted")
- ch: Channel number for the analysis (e.g., 1)
- grid: Grid resolution for spatial analysis (choices: 0.5, 1, 2; default: 1)
- region: Region code for map plot:
1 – Global
2 – Polar region (+60° latitude and above)
3 – Northern mid-latitudes (20° to 60° latitude)
4 – Tropics (-20° to 20° latitude)
5 – Southern mid-latitudes (-60° to -20° latitude)
6 – Southern polar region (below -60° latitude)
- sdate / -edate: Start and end dates for the time period (e.g., "2023-01-27" to "2023-01-28")
These parameters allow you to customize the spatial and temporal analysis to suit specific data and regions.

   

 ## Notes
 Ensure that the `obs_files_path` and `output_path` variables are correctly set to the paths of observational files and output directory, respectively.
 Adjust method parameters and plotting settings as needed for your specific use case.
 Make sure to define the `filter_by_variable` method as needed for filtering observational data based on variable values.

## Example Usage

```python
python SpatialTemporalStats.py -input /PATH/TO/INPUT/DIAG/FILES -output ./Results -sensor "atms_n20" -var "Obs_Minus_Forecast_adjusted" -ch 1 -grid 2 -region 1 -sdate "2023-01-27" -edate "2023-01-28"
```

## Example output plots using different settings
```python
-sensor "atms_n20" -var "Obs_Minus_Forecast_adjusted" -ch 1 -grid 2 -region 1 -sdate "2023-01-27" -edate "2023-01-28"
```
![atms_n20_ch1_Obs_Minus_Forecast_adjusted_Average_region_1](https://github.com/user-attachments/assets/e0ddcf64-8ce1-4175-b646-71d1d38ec3d4)
![atms_n20_ch1_Obs_Minus_Forecast_adjusted_Count_region_1](https://github.com/user-attachments/assets/a33dd6c4-bfb0-4ae9-a46d-02086f7dc960)
![atms_n20_ch1_Obs_Minus_Forecast_adjusted_RMS_region_1](https://github.com/user-attachments/assets/f9b34e74-7511-464d-a27d-82f08cfa5c6b)


Example plot for filtering out the locations where the land fraction is less than 0.9
```python
 -filter_by_vars Land_Fraction,lt,0.9
```
![atms_n20_ch1_Obs_Minus_Forecast_adjusted_Average_region_1](https://github.com/user-attachments/assets/bc6b7215-9d26-41c8-b51d-0f51d42238c3)

Example of the summary plots:
![atms_n20_Obs_Minus_Forecast_adjusted_mean_std](https://github.com/user-attachments/assets/99b09315-1faa-4fd1-9c26-e7b591dba2fc)
![atms_n20_Obs_Minus_Forecast_adjusted_sumamryCounts](https://github.com/user-attachments/assets/449cd174-f50d-4521-ab9f-e0d4b6f5ad9b)







