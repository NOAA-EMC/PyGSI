### April 2024
### Azadeh Gholoubi
# Python Tool for Time/Space (2D) Evaluation

## Overview
This tool provides functionalities for processing and analyzing data over time and space.

The `SpatialTemporalStats` class is designed to perform spatial and temporal statistics of data stored in NetCDF files. It includes features for generating grids, reading observational values, filtering data, plotting observations, and creating summary plots based on user settings.

### Important Methods of the SpatialTemporalStats Class
- `generate_grid(resolution=1)`: Generates a grid for spatial averaging based on the specified resolution. (default resolution is 1X1)
- `read_obs_values()`: Reads observational values from NetCDF files, filters them based on various criteria, performs spatial averaging, and returns the averaged values.
- `plot_obs()`: Plots observational data on a map, showing different regions and their corresponding data values.
- `list_variable_names(file_path)`: Lists variable names from a NetCDF file.
- `make_summary_plots()`: Generates summary plots of observational data, including scatter plots of counts, means, and standard deviations.

## Requirements
User need to load EVA environment when working on Hera, use the following commands:
```
cd GDASApp/modulefiles/
module load EVA/hera

```

## Usage
`user_Analysis.py` contains the `SpatialTemporalStats` class, which encapsulates the functionalities of the tool. Here's how to use it:

1. Import the `SpatialTemporalStats` class:

   ```python
   from SpatialTemporalStats import SpatialTemporalStats
2. Create an instance of the SpatialTemporalStats class:

   ```python
   my_tool = SpatialTemporalStats()

3. Specify the parameters based on the type of plots that you want:
   
  - `input_path`: Directory for input .nc files
  - `output_path`: Path to output plots
  - `sensor`: Sensor name
  - `channel_no`: Channel number (e.g., 1, 2, 3, 5)
  - `var_name`: variable name 
  - `start_date, end_date`: Start and End date of the input files for evaluations
  - `region`: Insert a number to select Global or Regional ouput plots (1: global (default), 2: polar region, 3: mid-latitudes region, 4: tropics region, 5: southern mid-latitudes region, 6: southern polar region)
  - `resolution`: Resolution for grid generation (1: 1X1 degree(default), 2:2X2 degree, 3:3X3 degree)
  - `filter_by_vars`: Filter by variable to generate plots based on surface type  (land, water, snow, seaice) or can be an empty list for no filtering.

4. Call `read_obs_values` to Read observational values and perform analysis:
   
```python
o_minus_f_gdf = my_tool.read_obs_values(
    input_path,
    sensor,
    var_name,
    channel_no,
    start_date,
    end_date,
    filter_by_vars,
    QC_filter)
```   
5. Call `plot_obs` to plot evaluation plots based on your setting for grid size, channel, region, surface type, and filtering values:

```python
my_tool.plot_obs(o_minus_f_gdf, var_name, region, resolution, output_path)
```
6. Call `make_summary_plots` to generate summary plots:

```python
summary_results = my_tool.make_summary_plots(
    input_path, sensor, var_name, start_date, end_date, QC_filter, output_path
)
```
 ## Notes
 Ensure that the `obs_files_path` and `output_path` variables are correctly set to the paths of observational files and output directory, respectively.
 Adjust method parameters and plotting settings as needed for your specific use case.
 Make sure to define the `filter_by_variable` method as needed for filtering observational data based on variable values.

To run the tool:

```
python user_Analysis.py

```

## Example Usage
Here's a sample script demonstrating how to use the`SpatialTemporalStats` tool:
![image](https://github.com/NOAA-EMC/PyGSI/assets/51101867/4379cb6e-e1a7-4167-8859-ae881f2c61c1)

## Example output plots using different settings
```python
var_name = "Obs_Minus_Forecast_adjusted"
region = 1
resolution = 2
filter_by_vars=[]
```
Calling `read_obs_values` and then `my_tool.plot_obs()` method will produce three plots for ave,count, rms as shown below:
![atms_n20_ch1_Obs_Minus_Forecast_adjusted_Average_region_1](https://github.com/NOAA-EMC/PyGSI/assets/51101867/b838ae92-3303-45ca-b7ba-35b11c01213c)
![atms_n20_ch1_Obs_Minus_Forecast_adjusted_Count_region_1](https://github.com/NOAA-EMC/PyGSI/assets/51101867/113ef427-9771-462a-b543-f36166ed978e)
![atms_n20_ch1_Obs_Minus_Forecast_adjusted_RMS_region_1](https://github.com/NOAA-EMC/PyGSI/assets/51101867/ed4bc44c-6364-451b-811e-b2c8a0ce5d2a)

Example plot for filtering out the locations where the land fraction is less than 0.9
```python
filter_by_vars = [("Land_Fraction", "lt", 0.9),]
```
![atms_n20_ch1_Obs_Minus_Forecast_adjusted_Average_region_1](https://github.com/NOAA-EMC/PyGSI/assets/51101867/978e2677-4a7b-45b3-a2e2-67674bf0803e)

Calling read_obs_values and then my_tool.make_summary_plots() method will generate two summary plots:
![atms_n20_Obs_Minus_Forecast_adjusted_mean_std](https://github.com/NOAA-EMC/PyGSI/assets/51101867/28cc26f4-c024-4713-82e1-b9a7ed5f5d1b)
![atms_n20_Obs_Minus_Forecast_adjusted_sumamryCounts](https://github.com/NOAA-EMC/PyGSI/assets/51101867/fd835f41-5b9c-4a14-be85-4c74d49571f6)





