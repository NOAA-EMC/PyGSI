### April 2024
### Azadeh Gholoubi
# Python Tool to do time/space (2D) evaluation

## Overview
This tool provides functionalities for processing and analyzing data over time and space.
It includes features for generating grids, reading observational values, filtering data, plotting observations, and creating summary plots.

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
  
To run the tool:

```
python user_Analysis.py

```
