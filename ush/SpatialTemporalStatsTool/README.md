### April 2024
### Azadeh Gholoubi
# Python Tool to do time/space (2D) evaluation

## Overview
This tool provides functionalities for processing and analyzing data over time and space.

The `SpatialTemporalStats` class is designed to perform spatial and temporal averaging of observational data stored in NetCDF files. It includes features for generating grids, reading observational values, filtering data, plotting observations, and creating summary plots.

### Important Methods
- `generate_grid(resolution=1)`
This method generates a grid for spatial averaging based on the specified resolution.

- `read_obs_values(obs_files_path, sensor, var_name, channel_no, start_date, end_date,
 filter_by_vars, QC_filter)`
This method reads observational values from NetCDF files, filters them based on various criteria, performs spatial averaging, and returns the averaged values.

- `plot_obs(selected_var_gdf, var_name, region, resolution, output_path)`
This method plots observational data on a map, showing different regions and their corresponding data values.

- `list_variable_names(file_path)`
This method lists variable names from a NetCDF file.

- `make_summary_plots(obs_files_path, sensor, var_name, start_date, end_date, QC_filter, output_path)`
This method generates summary plots of observational data, including scatter plots of counts, means, and standard deviations.

## Requirements
User need to load EVA environment when working on Hera, use the following commands:
```
cd GDASApp/modulefiles/
module load EVA/hera

```

## Usage
Before running the script, the user needs to specify certain parameters in `user_Analysis.py`:

- `input_path`: Directory for input .nc files
- `output_path`: Path to output plots
- `sensor`: Sensor name
- `Channels`: Channel number (e.g., 1, 2, 3, 5)
- `var_name`: variable name 
- `start_date, end_date`: Start and End date of the input files for evaluations
- `region`: Insert a number to select Global or Regional ouput plots (1: global (default), 2: polar region, 3: mid-latitudes region, 4: tropics region, 5: southern mid-latitudes region, 6: southern polar region)
- `resolution`: Resolution for grid generation (1: 1X1 degree(default), 2:2X2 degree, 3:3X3 degree)
- `filter_by_vars`: Filter by variable to generate plots based on surface type (land, water, snow, seaice) or can be an empty list for no filtering.
  
 ## Notes
 Ensure that the `obs_files_path` and `output_path` variables are correctly set to the paths of observational files and output directory, respectively.
 Adjust method parameters and plotting settings as needed for your specific use case.
 Make sure to define the `filter_by_variable` method as needed for filtering observational data based on variable values.

To run the tool:

```
python user_Analysis.py

```

## Example Usage
To run the `SpatialTemporalStats` tool with your specific settings, you can use the following example code:
![image](https://github.com/NOAA-EMC/PyGSI/assets/51101867/4379cb6e-e1a7-4167-8859-ae881f2c61c1)

## Example output plots
`var_name = "Obs_Minus_Forecast_adjusted"`

`region = 1`

`resolution = 2`

calling `my_tool.plot_obs()` method will produce three plots for ave,count, rms as shown below:
![atms_n20_ch1_Obs_Minus_Forecast_adjusted_Average_region_1](https://github.com/NOAA-EMC/PyGSI/assets/51101867/b838ae92-3303-45ca-b7ba-35b11c01213c)
![atms_n20_ch1_Obs_Minus_Forecast_adjusted_Count_region_1](https://github.com/NOAA-EMC/PyGSI/assets/51101867/113ef427-9771-462a-b543-f36166ed978e)
![atms_n20_ch1_Obs_Minus_Forecast_adjusted_RMS_region_1](https://github.com/NOAA-EMC/PyGSI/assets/51101867/ed4bc44c-6364-451b-811e-b2c8a0ce5d2a)

calling `my_tool.make_summary_plots()` method will generate two summary plots:
![atms_n20_Obs_Minus_Forecast_adjusted_mean_std](https://github.com/NOAA-EMC/PyGSI/assets/51101867/28cc26f4-c024-4713-82e1-b9a7ed5f5d1b)
![atms_n20_Obs_Minus_Forecast_adjusted_sumamryCounts](https://github.com/NOAA-EMC/PyGSI/assets/51101867/fd835f41-5b9c-4a14-be85-4c74d49571f6)





