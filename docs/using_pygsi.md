# Using EMCPy

## Using `diags.py`

[`diags.py`](https://github.com/NOAA-EMC/PyGSI/blob/develop/src/pyGSI/diags.py) is the main driver script for extracting data from GSI diagnostic netCDF4 files. It performs all of the heavy lifting to read, filter, and perform quality control of the data so the user does not have to do this themselves. Given the data type and a diagnostic file, `diags.py` will return the user a [Pandas](https://pandas.pydata.org/docs/index.html) dataframe they can use to perform their analysis.

### Getting Started

`diags.py` has three classes, `Conventional`, `Radiance`, and `Ozone`. Each class handles specific input netCDF4 files based on the data type.

### Conventional

If the user is interested in using conventional data, the following code would return a dataframe including all data from the inputted netCDF4 file:

```
from PyGSI.diags import Conventional

file = "path/to/conventional_gsi_file.nc4"

diag = Conventional(diag_file)
df = diag.get_data()
```


**`get_data()`**:

A method within the `Conventional()` class. Given input parameters, returns indexed dataframe from a conventional diagnostic file.

```
Args:
    obsid        : (list of ints; default=None) observation type ID number;
                   default=None
    subtype      : (list of ints; default=None) observation measurement ID
                   subtype number, default=None
    station_id   : (list of str; default=None) station id, default=None
    analysis_use : (bool; default=False) if True, will return
                   three sets of data:
                   assimilated (analysis_use_flag=1, qc<7),
                   rejected (analysis_use_flag=-1, qc<8),
                   monitored (analysis_use_flag=-1, qc>7)
    lvls         : (list type; default=None) List of pressure or height
                   levels i.e. [250,500,750,1000]. List must be arranged
                   low to high. For pressure, will return a dictionary of data
                   with data greater than the low bound, and less than or equal
                   to the high bound. For height, will return a dictionary of 
                   data with data greater than or equal to the low bound, and
                   less than the high bound.
    lvl_type     : (str; default='pressure') lvls definition as
                   'pressure' or 'height'.
Returns:
    indexed_df   : requested indexed dataframe
```

**`list_obsids()`**:

Prints all of the unique observation IDs in the diagnostic file.

```
diag.list_obsids()
```

**`list_stationids()`**:

Prints all the unique station IDs in the diagnostic file.

```
diag.list_stationids()
```

### Radiance

If the user is interested in using radiance data, the following code would return a dataframe including all data from the inputted netCDF4 file:

```
from PyGSI.diags import Radiance

file = "path/to/radiance_gsi_file.nc4"

diag = Radiance(diag_file)
df = diag.get_data()
```

**`get_data()`**:

A method within the `Radiance()` class. Given input parameters, returns indexed dataframe from a radiance diagnostic file.

```
Args:
    channel           : (list of ints; default=None) observation channel number
    qcflag            : (list of ints; default=None) qc flag number
    analysis_use      : (bool; default=False) if True, will return three sets of data:
                        assimilated (QC_Flag=0, inv_observation_error!=0),
                        rejected (QC_Flag!=0), monitored (use_flag!=1)
    separate_channels : (bool; default=False) if True, returns dict of separate
                        data by specified channels
    separate_qc       : (bool; default=False) if True, returns dict of separate
                        data by specified qc flag
    use_flag          : (bool; default=False) if True, will only return where
                        use_flag==1
    errcheck          : (bool; default=True) when True and qcflag==0, will
                        toss out obs where inverse obs error is zero (i.e.
                        not assimilated in GSI)
Returns:
    indexed_df        : requested indexed dataframe
```

**`list_channels()`**:

Prints all of the unique channels in the diagnostic file.

```
diag.list_channels()
```

**`list_qcflags()`**:

Prints all the unique qcflags in the diagnostic file.

```
diag.list_qcflags()
```

### Ozone

If the user is interested in using ozone data, the following code would return a dataframe including all data from the inputted netCDF4 file:

```
from PyGSI.diags import Ozone

file = "path/to/ozone_gsi_file.nc4"

diag = ozone(diag_file)
df_dict = diag.get_data()
```

**`get_data()`**:

A method within the `Ozone()` class. Given input parameters, returns a dictionary of indexed dataframes from an ozone diagnostic file.

```
Args:
    analysis_use : (bool; default=False) if True, will return
                   two sets of data:
                   assimilated (analysis_use_flag=1),
                   monitored (analysis_use_flag=-1)
    errcheck     : (bool; default=True) when True, will toss out
                   obs where inverse obs error is zero (i.e.
                   not assimilated in GSI)
Returns:
    df_dict      : (dict) requested indexed dataframes within a
                   dictionary
```

**`list_pressures()`**:

Prints all the unique pressure levels in the diagnostic file.

```
diag.list_pressures()
```
